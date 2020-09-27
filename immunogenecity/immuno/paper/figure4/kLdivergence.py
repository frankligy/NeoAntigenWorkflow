import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import collections
import scipy.stats as sc

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42



def peptide_and_hla(peptide,hla_type,hla, dic_inventory):
    length = len(peptide)
    if length == 10:
        peptide = peptide
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]

    dic = paratope_dic(hla)
    try:
        seq = dic[hla_type]
    except KeyError:
        hla_type = rescue_unknown_hla(hla_type, dic_inventory)
        seq = dic[hla_type]

    total = peptide + seq   # concatenate them together
    return total

def construct_aa_array(df):
    slots = np.empty([df.shape[0],56],dtype='<U1')
    for i in range(df.shape[0]):
        peptide = df.iloc[i]['peptide']
        hla_type = df.iloc[i]['HLA']
        seq = peptide_and_hla(peptide,hla_type,hla,dic_inventory)
        seq = seq.replace('X','-')
        slots[i,:] = list(seq)
    return slots



def array_to_prob(array):

    length = len(array)
    counter = collections.Counter(array)
    return {k:v/length for k,v in counter.items()}

def dic_to_list(dic):
    # let's set the rule: ARNDCQEGHILKMFPSTWYV-
    lis = []
    rule = 'ARNDCQEGHILKMFPSTWYV-'
    for i in rule:
        try:
            lis.append(dic[i])
        except KeyError:
            lis.append(0.00001)
    return np.array(lis)




def convert_to_prob(array):
    slots = np.empty([21,56],dtype=np.float32)
    for j in range(array.shape[1]):
        dic = array_to_prob(array[:,j])
        lis = dic_to_list(dic)
        slots[:,j] = lis
    return slots


def KL(neg_prob,pos_prob):
    slots = np.empty(56,dtype=np.float32)
    for j in range(neg_prob.shape[1]):
        neg = neg_prob[:,j]
        pos = pos_prob[:,j]
        slots[j] = sc.entropy(pos,neg)
    return slots


def array_to_fasta(array,name):
    with open(name,'w') as f:
        for row in array:
            seq = ''.join(list(row))
            seq.replace('-','X')
            f.write('>header\n')
            f.write('{}\n'.format(seq))

def array_to_fasta_peptide(array,name):
    with open(name,'w') as f:
        for row in array:
            seq = ''.join(list(row)[0:10])
            seq.replace('-','X')
            f.write('>header\n')
            f.write('{}\n'.format(seq))

def dict_inventory(inventory):
    dicA, dicB, dicC = {}, {}, {}
    dic = {'A': dicA, 'B': dicB, 'C': dicC}

    for hla in inventory:
        type_ = hla[4]  # A,B,C
        first2 = hla[6:8]  # 01
        last2 = hla[8:]  # 01
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)

    return dic


if __name__ == '__main__':
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    immuno_training = pd.read_csv('data/shuffle_training_test.txt',sep='\t')
    immuno_val = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')
    immuno = pd.concat([immuno_training,immuno_val])
    neg,pos = immuno.groupby(by=['immunogenecity'])
    neg,pos = neg[1],pos[1]

    neg_array = construct_aa_array(neg)
    pos_array = construct_aa_array(pos)
    neg_dis = convert_to_prob(neg_array)
    pos_dis = convert_to_prob(pos_array)

    # compute KL divergence
    realKL = KL(neg_dis,pos_dis)
    # draw weblogo
    array_to_fasta(neg_array,'real_neg.fasta')
    array_to_fasta(pos_array,'real_pos.fasta')

    array_to_fasta_peptide(pos_array,'paper/figure4/real_pos_pep.fasta')



    # what the model learned, we use binding data training set
    ResLikeCNN = model()
    ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')
    binding = pd.read_csv('data/transfer_training.txt',sep='\t')

    dataset = construct(binding, hla, dic_inventory)   # [ (10,21,1),(46,21,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    result_kl_blosum = ResLikeCNN.predict(x=[input1, input2])

    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    ResLikeCNN_index = model_aaindex()
    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')
    dataset = construct_aaindex(binding,hla,dic_inventory,after_pca)
    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    result_kl_aaindex = ResLikeCNN_index.predict(x=[input1,input2])
    result_kl_ensembl = np.mean(np.concatenate([result_kl_aaindex,result_kl_blosum],axis=1),axis=1)
    result_df = pd.DataFrame({'aaindex':result_kl_aaindex[:,0],'blosum':result_kl_blosum[:,0],'ensemble':result_kl_ensembl})
    result_df.to_csv('paper/figure4/binding_non_overlapping_data.txt',sep='\t',index=None)

    result_df = pd.read_csv('paper/figure4/binding_non_overlapping_data.txt',sep='\t')
    result_kl_ensembl = result_df['ensemble']

    hard = [1 if i > 0.5 else 0 for i in result_kl_ensembl]
    binding['immunogenecity'] = hard
    neg_r,pos_r = binding.groupby(by=['immunogenecity'])
    neg_r,pos_r = neg_r[1],pos_r[1]

    neg_array_r = construct_aa_array(neg_r)
    pos_array_r = construct_aa_array(pos_r)
    neg_dis_r = convert_to_prob(neg_array_r)
    pos_dis_r = convert_to_prob(pos_array_r)

    # compute KL divergence
    learnedKL = KL(neg_dis_r,pos_dis_r)
    # draw weblogo
    array_to_fasta(neg_array_r,'learn_neg.fasta')
    array_to_fasta(pos_array_r,'learn_pos.fasta')

    array_to_fasta_peptide(pos_array_r,'paper/figure4/learn_pos_pep.fasta')

    # compare
    plt.bar(np.arange(56)+1,realKL,color='blue',alpha=0.3)
    plt.bar(np.arange(56)+1, learnedKL, color='orange',alpha=0.3)

    # plot KL plot
    plt.bar(np.arange(10)+1,realKL[0:10],color='blue',alpha=0.3)
    plt.xticks(np.arange(10)+1,labels=['P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'])
    plt.ylabel('Kullback-Leibler divergence')
    plt.xlabel('Peptide Positions')
    plt.savefig('paper/figure4/realKL.pdf')

    plt.bar(np.arange(10)+1,learnedKL[0:10],color='orange',alpha=0.3)
    plt.xticks(np.arange(10)+1,labels=['P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'])
    plt.ylabel('Kullback-Leibler divergence')
    plt.xlabel('Peptide Positions')
    plt.savefig('paper/figure4/learnedKL.pdf')

    '''
    Go to weblogo threeplus web server:
    size: large
    y_axis sclae =2.0
    y_ticks space = 0.5
    color: chemistryAA
    '''


