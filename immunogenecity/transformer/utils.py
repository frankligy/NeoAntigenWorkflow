#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:25:58 2020

@author: ligk2e
"""
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SubsMat import MatrixInfo
import json
import numpy as np
from collections import Counter
import pandas as pd
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split,Subset


def json2fsa(hla):
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope/{}.json'.format(hla),'r') as f:
        data = json.load(f)
        
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/{}.fsa'.format(hla),'w') as f:        
        for item in data:
            key = list(item.keys())[0]
            value = list(item.values())[0]
            f.write('>{}\n'.format(key))
            f.write('{}\n'.format(value))

def multiple_json2fsa():
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope/inventory_new.txt','r') as f1:
        for line in f1:
            line = line.rstrip('\n')
            json2fsa(line)
'''          
run clustal-omega:
    download binary: http://www.clustal.org/omega/
    chmod 777 ./clustal-omega-1.2.3-macosx
    ./clustal-omega-1.2.3-macosx -i "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/HLA-A*0101.fsa" -o "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/HLA-A*0101.aligned.fasta" --auto -v            
   
    
run multiple sequentially:
    cat /Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope/inventory_new.txt | while read line; do
    ./clustal-omega-1.2.3-macosx -i "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/${line}.fsa" -o "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/${line}.aligned.fasta" --auto -v; done         
    
only hla that has more than 1 paratope will be processed in clustal-omega   
'''

def single_paratope(hla):
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/{}.fsa'.format(hla),'r') as f:
        seq = f.readlines()[1].rstrip('\n')
    return hla,seq


def matrix2concensus(mat):
    final = ''
    for j in range(mat.shape[1]):
        most = Counter(mat[:,j]).most_common(1)[0][0]  # if most_common(2): [('A', 3), ('C', 1)]
        if most == '-':
            most = Counter(mat[:,j]).most_common(2)[1][0]
        final += most
    return final
    
    

def msa_paratope(hla):
    alignment = AlignIO.read(open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/{}.aligned.fasta'.format(hla)),"fasta")
    msa = []
    for record in alignment:
        msa.append(list(record.seq))   # another part is record.id
    mat = np.array(msa)
    final = matrix2concensus(mat)
    return hla,final
    
    
        
    

def hla_paratope():

            
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/inventory_single.txt','r') as single:
        singles = single.readlines()  # each one will contain '\n'
        singles = [item.rstrip('\n') for item in singles]
    
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/inventory_msa.txt','r') as multiple:
        multiples = multiple.readlines()
        multiples = [item.rstrip('\n') for item in multiples]
        
    
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable.txt','w') as f:
        for item in singles:
            hla,seq = single_paratope(item)   
            f.write('{0}\t{1}\n'.format(hla,seq))
        for item in multiples:
            hla,seq = msa_paratope(item)
            f.write('{0}\t{1}\n'.format(hla,seq))
            
########################################################################################################            

def clean_series(series):  # give a pandas series
    
    if series.dtype == object:  # pandas will store str as object since string has variable length, you can use astype('|S')
        clean = []
        for item in series:
            item = item.lstrip(' ')   # remove leading whitespace
            item = item.rstrip(' ')   # remove trailing whitespace
            item = item.replace(' ','')  # replace all whitespace in the middle
            clean.append(item)
    else:
        clean = series
        
    
        
    return pd.Series(clean)


def clean_data_frame(data):  # give a pandas dataFrame
    
    peptide_clean = clean_series(data['peptide'])
    hla_clean = clean_series(data['HLA'])
    immunogenecity_clean = clean_series(data['immunogenecity'])
    
    data_clean = pd.concat([peptide_clean,hla_clean,immunogenecity_clean],axis=1)
    data_clean.columns = ['peptide','HLA','immunogenecity']
    
    return data_clean

def convert_hla(hla):
    cond = True
    hla = hla.replace(':','')
    if len(hla) < 9: cond = False   # HLA-A3
    elif len(hla) == 9:   # HLA-A3002
        f = hla[0:5]  # HLA-A
        e = hla[5:]   # 3002
        hla = f+'*'+e
    return hla,cond

def convert_hla_series(df):
    new = []
    col = []
    for i in df['HLA']:
        hla,cond = convert_hla(i)
        col.append(cond)
        if cond == True: new.append(hla)
    df = df.loc[pd.Series(col)]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    df['HLA'] = new
    return df
        
def test_no_space(series):
    for i in series:
        if ' ' in i:
            print('damn')

   
'''
a = pd.read_excel('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data.xlsx')
a1 = clean_data_frame(a)
test_no_space(a1.iloc[:,0])
test_no_space(a1.iloc[:,1])
a1.iloc[:,2].dtype

a2 = convert_hla_series(a1)
a2.to_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data.txt',sep='\t',index=None)

then use:
    
{ cat data.txt | head -n 1; cat data.txt | tail -n +2 | sort -u -k1,2; } > data_new.txt, only 32669 training data left

ori = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data_new.txt',sep='\t')
hla = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable.txt',sep='\t',header=None,names=['hla','paratope'])
inventory = hla['hla']
dic_inventory = dict_inventory(inventory)


'''
def dict_inventory(inventory):
    dicA,dicB,dicC = {},{},{}
    dic = {'A':dicA,'B':dicB,'C':dicC}
    
    for hla in inventory:
        type_ = hla[4]  # A,B,C
        first2 = hla[6:8] # 01
        last2 = hla[8:]  # 01
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)
            
    return dic


def rescue_unknown_hla(hla,dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2)-int(i)) for i in small_category]
        optimal = min(zip(small_category,distance),key=lambda x:x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2)-int(i)) for i in small_category]   
        optimal = min(zip(small_category,distance),key=lambda x:x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + str(big_category[optimal][0])

class dataset(Dataset):
    # the output would be ([seq_len,21],[batch]),(),()
    def __init__(self,ori,hla,dic_inventory):
        self.ori = ori
        self.hla = hla
        self.dic_inventory = dic_inventory
        
        self.paratope_dic()
        self.middle =  self.convert()
        self.new = self.padding()
        #self.new = self.padding_onehot()
        
    def __len__(self):
        return len(self.new)
    
    def __getitem__(self,idx):
        return self.new[idx]
    
    
    def padding(self):
        len_values = [tup[0].shape[0] for tup in self.middle]
        #max_length = max(len_values)
        max_length = 50
        
        # padding
        bucket = []
        for item in self.middle:

            length = item[0].shape[0]
            gap = max_length - length
            if gap % 2 == 0:  # even number
                gapped_left, gapped_right = gap // 2, gap //2  # will be an int
            else:  # odd number
                if np.random.uniform() < 0.5:  # randomly decide which side will have one more padded value
                    gapped_left = gap // 2
                    gapped_right = gap - gapped_left
                else:
                    gapped_right = gap // 2
                    gapped_left = gap - gapped_right
                    
            padding_left = torch.empty([gapped_left,20]).fill_(-1.0)
            padding_right = torch.empty([gapped_right,20]).fill_(-1.0)
            final = torch.cat([padding_left,item[0],padding_right],dim=0)
            bucket.append((final,item[1])) 

        
        self.max_length = max_length
        
        return bucket
    
    def padding_onehot(self):
        len_values = [tup[0].shape[0] for tup in self.middle]
        max_length = max(len_values)
        #max_length = 48
        
        # padding
        bucket = []
        for item in self.middle:

            length = item[0].shape[0]
            gap = max_length - length
            if gap % 2 == 0:  # even number
                gapped_left, gapped_right = gap // 2, gap //2  # will be an int
            else:  # odd number
                if np.random.uniform() < 0.5:  # randomly decide which side will have one more padded value
                    gapped_left = gap // 2
                    gapped_right = gap - gapped_left
                else:
                    gapped_right = gap // 2
                    gapped_left = gap - gapped_right
                    
            padding_left = torch.empty([gapped_left,20]).fill_(0.05)
            padding_right = torch.empty([gapped_right,20]).fill_(0.05)
            final = torch.cat([padding_left,item[0],padding_right],dim=0)
            bucket.append((final,item[1])) 

        
        self.max_length = max_length
        
        return bucket



    def paratope_dic(self):
        df = self.hla
        self.dic = {}
        for i in range(df.shape[0]):
            hla = df['hla'].iloc[i]
            paratope = df['paratope'].iloc[i]
            self.dic[hla] = paratope
    
    @staticmethod
    def onehot_classic(peptide):
        amino = 'ARNDCQEGHILKMFPSTWYV'
        encoded = torch.empty([len(peptide),20])
        onehot = torch.eye(20)
        for i in range(len(peptide)):
            encoded[i,:] = onehot[:,amino.index(peptide[i])]

        return encoded

    @staticmethod
    def onehot_adapt(peptide):
        amino = 'ARNDCQEGHILKMFPSTWYV'
        encoded = torch.empty([len(peptide),20])
        onehot = torch.eye(20)
        mask = torch.eye(20)
        onehot = onehot.masked_fill(mask == 1, 0.9)
        onehot = onehot.masked_fill(mask == 0, 0.005)
        for i in range(len(peptide)):
            encoded[i,:] = onehot[:,amino.index(peptide[i])]
        return encoded




    @staticmethod
    def blosum50(peptide):
        amino = 'ARNDCQEGHILKMFPSTWYV'
        dic = MatrixInfo.blosum50
        matrix = np.zeros([20,20])
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                try:
                    matrix[i,j] = dic[(amino[i],amino[j])] 
                except KeyError:
                    matrix[i,j] = dic[(amino[j],amino[i])]
                    
        encoded = torch.empty([len(peptide),20])       # (seq_len,20)       
        for i in range(len(peptide)):

            encoded[i,:] = torch.from_numpy(matrix[:,amino.index(peptide[i])])
                
        return encoded
    
    def convert(self):
        lis = []
        df = self.ori
        for i in range(df.shape[0]):
            #print(i)
            peptide = df['peptide'].iloc[i]
            hla_type = df['HLA'].iloc[i]
            immuno = df['immunogenecity'].iloc[i]
            try:
                cat = self.dic[hla_type] + peptide
            except KeyError:
                hla_type = rescue_unknown_hla(hla_type, self.dic_inventory)
                cat = self.dic[hla_type] + peptide
            cat = cat.upper()
            if 'X' in cat: continue
            X = dataset.blosum50(cat).float()   # 2-d tensor
            #X = dataset.onehot_classic(cat).float()
            y = torch.tensor(immuno).long()  # 0-d tensor
            lis.append((X,y))
        return lis
            
            
            

def balancedBinaryLoader(dataset,batch_size):
    
    dic = {'0':[],'1':[]}
    for i in range(len(dataset)):
        X = dataset[i][0]
        y = dataset[i][1]
        if y == 1: dic['1'].append(i)
        elif y == 0: dic['0'].append(i)
        
    

    sample_size = batch_size // 2  # will be an int, make sure batch_size is an even number
    
    negative = Subset(dataset,dic['0']) # dataset.Subset object
    positive = Subset(dataset,dic['1'])
    # print(len(positive),type(positive)) 
    
    negative_loader = DataLoader(negative,batch_size=sample_size,shuffle=True,drop_last=True)
    positive_loader = DataLoader(positive,batch_size=sample_size,shuffle=True,drop_last=True)
    



    neg_chunks_X = []
    neg_chunks_y = []
    for idx,(X,y) in enumerate(negative_loader):
        neg_chunks_X.append(X)
        neg_chunks_y.append(y)

    
    pos_chunks_X = []
    pos_chunks_y = []
    for idx,(X,y) in enumerate(positive_loader):
        pos_chunks_X.append(X)
        pos_chunks_y.append(y)

    
    pos_chunks_X_cycle = pos_chunks_X * 10
    pos_chunks_y_cycle = pos_chunks_y * 10


    chunks_X_list = []
    chunks_y_list = []    
    for i in range(len(neg_chunks_X)):
        chunks_X = torch.cat([neg_chunks_X[i],pos_chunks_X_cycle[i]],dim=0)
        chunks_y = torch.cat([neg_chunks_y[i],pos_chunks_y_cycle[i]],dim=0)
        chunks_X_list.append(chunks_X)
        chunks_y_list.append(chunks_y)
        
    
        

    
    loader = list(zip(chunks_X_list,chunks_y_list)) # zip can only be iterated once
    return loader             

    
##############################################################################################
'''
all_assay = pd.read_csv('/Users/ligk2e/Desktop/immunogenecity/iedb/all_assay.csv',skiprows=0,header=1)
all_assay_extract = all_assay[['Description','Qualitative Measure','Allele Name']]

all_assay_extract = replace_immunogenecity_df(all_assay_extract)
all_assay_extract = replace_hla_allele_df(all_assay_extract)
all_assay_extract = purge_peptide_df(all_assay_extract)
all_assay_extract.columns = ['peptide','immunogenecity','HLA']
all_assay_extract = all_assay_extract[['peptide','HLA','immunogenecity']]
all_assay_extract.to_csv('/Users/ligk2e/Desktop/immunogenecity/iedb/all_assay.txt',sep='\t',index=None)

{ cat all_assay.txt | head -n 1; cat all_assay.txt | tail -n +2 | sort -u -k1,2; } > all_assay_new.txt
cat all_assay_new.txt | awk 'length($1) < 15 {print $0}' > all_assay_new_filster15.txt

on linux system:
cat all_assay_new_filster15.txt | tail -n +2 | shuf > shuffle_all.txt
{ echo -e "peptide\tHLA\timmunogenecity"; cat shuffle_all.txt | head -n 30000; } > shuffle_training.txt
{ echo -e "peptide\tHLA\timmunogenecity"; cat shuffle_all.txt | tail -n 3444; } > shuffle_testing.txt

'''

def purge_peptide(entry):
    cond = True
    if '+' in entry:
        cond = False
    return entry,cond
        
def purge_peptide_df(df):

    col = []
    conds = []


    for i in range(df.shape[0]):
        entry = df.iloc[i]['Description']
        entry,cond = purge_peptide(entry)
        col.append(entry)
        conds.append(cond)
    df.update(pd.DataFrame({'Description':col}))
    df = df.loc[pd.Series(conds)]  
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    return df


def replace_immunogenecity(entry):
    if entry == 'Positive' or entry == 'Positive-High' or entry == 'Positive-Intermediate' or entry == 'Positive-Low':
        entry = 1
    else:
        entry = 0
    return entry

            
            
def replace_immunogenecity_df(df):
    col = []
    for i in range(df.shape[0]):
        entry = df.iloc[i]['Qualitative Measure']
        entry = replace_immunogenecity(entry)
        col.append(entry)
    df.update(pd.DataFrame({'Qualitative Measure':col}))
    return df
                 
            
def replace_hla_allele(entry):
    cond = True
    entry = entry.replace(':','')
    if len(entry) < 9:
        a = entry[0:5]
        b = entry[5:]
        if len(b) == 1: b = '0' + b
        if 'w' in b: b = b.replace('w','0')
        entry = a + '*' + b + '01'
    if ' ' in entry: cond = False
    return entry,cond

def replace_hla_allele_df(df):
    col = []
    conds = []
    for i in range(df.shape[0]):
        entry = df.iloc[i]['Allele Name']
        entry,cond = replace_hla_allele(entry)
        col.append(entry)
        conds.append(cond)
    df.update(pd.DataFrame({'Allele Name':col}))
    df = df.loc[pd.Series(conds)]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))

    return df


######################################################################################
# preprocess AA index
'''
cplit aaindex3 '/\/\//+1' {45}    # 47 indices in total, repeat 45 more times, means do it 46 times in total, the remaing one will automatically be written into a new file xx46
​a=0; for i in xx*; do cat $i | tail -n 21 | head -n 20 > index$a.txt; a=$[$a+1]; done
for i in xx*; do grep -e "^D" $i; done
# index42,index43,index45 are non-symmetrical, discard them


# get [210 40] matrix
result = convert_all_index_before_pca()
execute:
    normalize_before_pca()
    scale_before_pca()
    pca_get_components()
    pca_apply_reduction()   # [210,25] matrix
    
'''


def count_character(string,tar):
    counter = 0
    for char in string:
        if char == tar:
            counter += 1
    return counter
               

def impute_triangle(int_path):
    with open(int_path,'r') as f, open('{0}.tmp.txt'.format(int_path),'w') as f1:
        data = f.readlines()
        for row in data:  # each item in the data array corrspond to a row in lower triangle matrix
            row = row.lstrip(' ')   # remove leading whitespace
            row = re.sub(' +',' ',row)  # compress multiple whitespace to one
            num_space = count_character(row,' ')   # how many whitespace
            num_newline = count_character(row,'\n')  # how many newline
            num = num_space + num_newline   # how many items
            diff = 20 - num
            row = row.rstrip('\n') + ' '  # will be '-3.4 -4.3 '
            row += '0.0 '*diff
            row = row.rstrip(' ') + '\n'  # will be '-3.4 -4.3 0.0 0.0 0.0\n'
            f1.write(row)

    index = np.loadtxt('{0}.tmp.txt'.format(int_path))
        
    return index
        
    
    
def extract_pair_metircs(index):
    amino = 'ARNDCQEGHILKMFPSTWYV'
    frame = np.empty([210,1])
    counter = -1
    for i in range(index.shape[0]):
        for j in range(index.shape[1]):
            entry = index[i,j]
            if not entry==0:
                counter += 1
                frame[counter] = entry
    return frame

def convert_all_index_before_pca():
    array = []
    for i in range(47):
        if not i in [42,43,45]:
            if not i in [2,22,23,25]:  # Glu is not available, NA, just discard for now
                index=impute_triangle('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/AAindex3/index{0}.txt'.format(i))
                frame=extract_pair_metircs(index)
                array.append(frame)
    result = np.concatenate(array,axis=1)
    return result   # should be [210,40]


def normalize_before_pca(result):
    # first figure out what columns have extreme value and needs to be normalizced beforehand
    troubles = [0,1,7,9,11,14,15,28,29,30]
    for i in troubles:
        subject = result[:,i].reshape(1,-1)   # Normalizer have to be applied to a matrix, and perform by row, each subject [1,210]
        t1 = Normalizer()
        new = t1.fit_transform(subject).reshape(-1)   # [210,] 1d
        result[:,i] = new
    return result
    
def scale_before_pca(result):
    t1 = StandardScaler()
    new = t1.fit_transform(result)
    return new
        
    

def pca_get_components(result):
    pca= PCA()
    pca.fit(result)
    result = pca.explained_variance_ratio_
    sum_ = 0
    for index,var in enumerate(result):
        sum_ += var
        if sum_ > 0.95:
            return index    # 25 components
        
def pca_apply_reduction(result):
    pca = PCA(n_components=25)   # or strictly speaking ,should be 26, since python is 0-index
    new = pca.fit_transform(result)
    return new
    
def wrapper_preprocess():
    result = convert_all_index_before_pca()
    result = normalize_before_pca(result)
    result = scale_before_pca(result)
    result = pca_apply_reduction(result)
    return result
    


def sum_to_itself(a):
    if a == 0:
        return 0
    elif a == 1:
        return 1
    elif a == 2:
        return 1+2
    else:
        return sum_to_itself(a-1) + a  
    
class CNN_dataset(Dataset):
    def __init__(self,ori,hla,dic_inventory,index):   # index is [210,25] matrix
        self.ori= ori
        self.hla = hla 
        self.dic_inventory = dic_inventory
        self.index = index
        
        self.paratope_dic()   # will get self.dic in order to get hla sequence in self.convert function
        self.new = self.convert()   # self.new [  (tensor(25,pep_max,hla_max),0-d tensor as y), (), ()     ]

    
    def __len__(self):
        return len(self.new)
    
    
    def __getitem__(self,index):
        return self.new[index]





    def paratope_dic(self):
        df = self.hla
        self.dic = {}
        for i in range(df.shape[0]):
            hla = df['hla'].iloc[i]
            paratope = df['paratope'].iloc[i]
            self.dic[hla] = paratope
    



    def get_index_value(self,tup): # input will be a tuple (a,b), a, b will be the amino acid one letter character like ('D','Q')
        amino = 'ARNDCQEGHILKMFPSTWYV'
        first = amino.index(tup[0].upper())
        second = amino.index(tup[1].upper())
        if first < second: 
            first,second = second,first  # first will be row index, second will be column index in [20*20] matrix
        row_index = sum_to_itself(first) + second  # the row index in [210,25] matrix
        values = self.index[row_index,:]
        return values   # 1d numpy array
        
            
        
    def convert(self):
        df = self.ori
        peptides = df['peptide']
        self.pep_max = max([len(pep) for pep in peptides])
        self.hla_max = 46   # hardcode it
        self.components = 25  # hardcode it

        
        result = torch.empty([25,self.pep_max,self.hla_max])
        
        '''
        clear them up:
        we initialize a empty 3-d array, the same shape as our desired output,
        then at each 2d plane, given [i,j] we should have a vector of length 25,
        the question is, how to find the corresponding row in [210,15] matrix?
        
        let's think about this in an example, say we wanna explore the row index of D -> Q
        amino = 'ARNDCQEGHILKMFPSTWYV'
        amino.index('D') will return 3, amino.index('Q') will return 5, considering the way we generate and order
        the 210 aa-pair, the way we get (D->Q) value should be when row is Q and column is D, extropolating that means
        always from larger index(Q,5) to smaller index(D,3), it will be clear when draw a 20*20 lower triangle and figure
        out how the 210 pair is genrated one by one. 
        
        then still focus on D -> Q example, how many pair been pulled out before Q to D, row-wise, there are 5 aa before
        Q because the index of Q is 5, so there are 1 + 2 + 3 + 4 + 5 pair already, then comes to Q's row, before arrive at
        D, there are 3 aa in front of that, because index of D is 3, so (1 + 2 + 3 + 4 + 5) + 3, there are 18 items before
        (D -> Q), given python is 0-index, so (D -> Q) is the 18th row in index_result [210,25] matrix
        
        Let's generalize the rule we observed, we can write the function
        
        
        '''
        final = []
        for row in range(df.shape[0]):
            print(row)
            peptide = df['peptide'].iloc[row]
            hla_type = df['HLA'].iloc[row]
            try:
                hla_seq = self.dic[hla_type]
            except KeyError:
                hla_type = rescue_unknown_hla(hla_type, self.dic_inventory)
                hla_seq = self.dic[hla_type]
                
            immuno = df['immunogenecity'].iloc[row]
            y = torch.tensor(immuno).long()
            
            
            '''
            Understanding the padding process below:
            say for a 8-mer peptide, we know self.pep_max = 14
            gap_left = 3, gap_right =3
            gap_left_indices: 0,1,2  these three indices will be padded instead of encoding with real value
            gap_right_indices: 11,12,13, so range(14), take from 14-3=11, we get 11,12,13
            
            '''
            pep_len = len(peptide)
            diff_len = self.pep_max - pep_len
            gap_left = diff_len // 2 
            gap_right = diff_len - gap_left
            gap_left_indices = list(range(gap_left))
            gap_right_indices = list(range(self.pep_max))[self.pep_max-gap_right:]
            
            
            
            
            
            for i in range(result.shape[1]):
                for j in range(result.shape[2]):
                    hla_index = j
                    hla_aa = hla_seq[hla_index]
                    if i in gap_left_indices or i in gap_right_indices:
                        result[:,i,j] = torch.empty(25).fill_(0.0).float()

                    elif hla_aa == '-':
                        result[:,i,j] = torch.empty(25).fill_(0.0).float()
                    else:
                        real_peptide_index = i - gap_left  # say the i=4, [0,1,2] is gapped, 4-3=1, the second aa in peptide
                        real_peptide = peptide[real_peptide_index]
                        if real_peptide == 'X':
                            result[:,i,j] = torch.empty(25).fill_(0.0).float()
                        else:
                            


                            try:
                                result[:,i,j] = torch.from_numpy(self.get_index_value((real_peptide,hla_aa))).float()
                            except: print(real_peptide,hla_aa); raise Exception
                            

            final.append((result,y))
            
        return final
            
            
######################################################################################
def pytorch_training(modelObj,training_dataset,optimizer,criterion,batch_size,num_epochs,outdir):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = modelObj().to(device)
    training_loader = DataLoader(training_dataset,batch_size=batch_size,shuffle=True)
 
    num_epochs = num_epochs
    for epoch in range(num_epochs):
        loss_list = []
        acc_list = []

        for i in training_loader:
            X = i[0].to(device)
            y = i[1].to(device)
            optimizer.zero_grad()
            
            y_pred = model(X)
            print(y_pred)

            loss = criterion(y_pred,y)
            loss.backward()
            optimizer.step()
            loss_list.append(loss.item())
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)
            print(predictions)
            print(y)
    
            num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it
    
            num_samples  += predictions.size(0)
    
            acc_list.append(float(num_correct)/float(num_samples)*100)
            
        loss,acc = sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    
        print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(epoch+1,num_epochs,loss,acc))

    torch.save(model.state_dict(),outdir)
    
                    
            
            
            
            
            
            



















        
    
    