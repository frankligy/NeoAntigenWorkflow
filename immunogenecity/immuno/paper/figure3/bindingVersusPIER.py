import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
from utils import *
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import matplotlib.pyplot as plt
import numpy as np


class ResBlock(layers.Layer):
    def __init__(self,in_channel,pool_size):
        super(ResBlock,self).__init__()
        intermediate_channel = in_channel
        out_channel = in_channel * 2
        self.conv1 = layers.Conv2D(filters=intermediate_channel,kernel_size=(1,1),strides=(1,1),padding='same')
        self.bn1 = layers.BatchNormalization()
        self.conv2 = layers.Conv2D(filters=intermediate_channel,kernel_size=(3,1),strides=(1,1),padding='same')
        self.bn2 = layers.BatchNormalization()
        self.conv3 = layers.Conv2D(filters=out_channel,kernel_size=(1,1),strides=(1,1),padding='same')
        self.bn3 = layers.BatchNormalization()
        self.identity = layers.Conv2D(filters=out_channel,kernel_size=(1,1),strides=(1,1))
        self.maxpool = layers.MaxPool2D(pool_size=pool_size,strides=pool_size)



    def call(self,x):
        out = keras.activations.relu(self.bn1(self.conv1(x)))   # (8,1,16)
        out = keras.activations.relu(self.bn2(self.conv2(out)))  # (8,1,16)
        out = keras.activations.relu(self.bn3(self.conv3(out)))   # (8,1,32)
        identity_map = self.identity(x)   # change original input (8,1,16)  --> (8,1,32)
        out = out + identity_map    # (8,1,32)
        out = self.maxpool(out)    # (4,1,32)

        return out


class CNN_peptide(layers.Layer):
    def __init__(self):
        super(CNN_peptide,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(3,21),strides=(1,1))
        self.block1 = ResBlock(16,(2,1))
        self.block2 = ResBlock(32,(2,1))
        self.block3 = ResBlock(64,(2,1))

    def call(self,x):    # (10,21,1)
        out = self.conv(x)   # (8,1,16)
        out = self.block1(out)   # (4,1,32)
        out = self.block2(out)   # (2,1,64)
        out = self.block3(out)   # (1,1,128)
        return out


class CNN_MHC(layers.Layer):
    def __init__(self):
        super(CNN_MHC,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(15,21),strides=(1,1)) # (32,1,16)
        self.block1 = ResBlock(16, (2, 1))    # (16,1,32)
        self.block2 = ResBlock(32, (2, 1))    # (8,1,64)
        self.block3 = ResBlock(64, (2, 1))    # (4,1,128)
        self.conv_add = layers.Conv2D(filters=128,kernel_size=(4,1),strides=(1,1))
        self.bn = layers.BatchNormalization()


    def call(self, x):
        out = self.conv(x)
        out = self.block1(out)
        out = self.block2(out)
        out = self.block3(out)
        out = keras.activations.relu(self.bn(self.conv_add(out)))   # (1,1,128)
        return out


class model(keras.Model):
    def __init__(self):
        super(model,self).__init__()
        self.br_pep = CNN_peptide()
        self.br_mhc = CNN_MHC()
        self.flatten = layers.Flatten()
        self.fc1 = layers.Dense(128,activation='relu')
        self.fc2 = layers.Dense(1,activation='sigmoid')

    def call(self,input):
        x1,x2 = input[0],input[1]  # x1: (10,21,1)    x2: (46,21,1)
        out1 = self.flatten(self.br_pep(x1))
        out2 = self.flatten(self.br_mhc(x2))
        out = layers.concatenate([out1,out2])
        out = self.fc1(out)
        out = self.fc2(out)
        return out

    def model(self):
        x1 = keras.Input(shape=(10,21,1))
        x2 = keras.Input(shape=(46,21,1))
        return keras.Model(inputs=[x1,x2],outputs=self.call([x1,x2]))



class CNN_peptide_aaindex(layers.Layer):
    def __init__(self):
        super(CNN_peptide_aaindex,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(3,12),strides=(1,1))
        self.block1 = ResBlock(16,(2,1))
        self.block2 = ResBlock(32,(2,1))
        self.block3 = ResBlock(64,(2,1))

    def call(self,x):    # (10,21,1)
        out = self.conv(x)   # (8,1,16)
        out = self.block1(out)   # (4,1,32)
        out = self.block2(out)   # (2,1,64)
        out = self.block3(out)   # (1,1,128)
        return out


class CNN_MHC_aaindex(layers.Layer):
    def __init__(self):
        super(CNN_MHC_aaindex,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(15,12),strides=(1,1)) # (32,1,16)
        self.block1 = ResBlock(16, (2, 1))    # (16,1,32)
        self.block2 = ResBlock(32, (2, 1))    # (8,1,64)
        self.block3 = ResBlock(64, (2, 1))    # (4,1,128)
        self.conv_add = layers.Conv2D(filters=128,kernel_size=(4,1),strides=(1,1))
        self.bn = layers.BatchNormalization()


    def call(self, x):
        out = self.conv(x)
        out = self.block1(out)
        out = self.block2(out)
        out = self.block3(out)
        out = keras.activations.relu(self.bn(self.conv_add(out)))   # (1,1,128)
        return out


class model_aaindex(keras.Model):
    def __init__(self):
        super(model_aaindex,self).__init__()
        self.br_pep = CNN_peptide_aaindex()
        self.br_mhc = CNN_MHC_aaindex()
        self.flatten = layers.Flatten()
        self.fc1 = layers.Dense(128,activation='relu')
        self.fc2 = layers.Dense(1,activation='sigmoid')

    def call(self,input):
        x1,x2 = input[0],input[1]  # x1: (10,12,1)    x2: (46,12,1)
        out1 = self.flatten(self.br_pep(x1))
        out2 = self.flatten(self.br_mhc(x2))
        out = layers.concatenate([out1,out2])
        out = self.fc1(out)
        out = self.fc2(out)
        return out

    def model(self):
        x1 = keras.Input(shape=(10,12,1))
        x2 = keras.Input(shape=(46,12,1))
        return keras.Model(inputs=[x1,x2],outputs=self.call([x1,x2]))

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


def rescue_unknown_hla(hla, dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    #print(hla)
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + str(big_category[optimal][0])


def blosum50_new(peptide):
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    dic = MatrixInfo.blosum50
    matrix = np.zeros([21, 21])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            try:
                matrix[i, j] = dic[(amino[i], amino[j])]
            except KeyError:
                try:
                    matrix[i, j] = dic[(amino[j], amino[i])]
                except:
                    matrix[i, j] = -1

    encoded = np.empty([len(peptide), 21])  # (seq_len,21)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded


def paratope_dic(hla):
    df = hla
    dic = {}
    for i in range(df.shape[0]):
        hla = df['hla'].iloc[i]
        paratope = df['paratope'].iloc[i]
        dic[hla] = paratope
    return dic


def peptide_data(peptide):   # return numpy array [10,21,1]
    length = len(peptide)
    if length == 10:
        encode = blosum50_new(peptide)
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]
        encode = blosum(peptide)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode


def hla_data(hla, dic_inventory, hla_type):    # return numpy array [46,21,1]
    dic = paratope_dic(hla)
    try:
        seq = dic[hla_type]
    except KeyError:
        hla_type = rescue_unknown_hla(hla_type, dic_inventory)
        seq = dic[hla_type]
    encode = blosum(seq)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode


def construct(ori, hla, dic_inventory):
    series = []
    for i in range(ori.shape[0]):
        peptide = ori['peptide'].iloc[i]
        hla_type = ori['HLA'].iloc[i]
        immuno = np.array(ori['immunogenecity'].iloc[i]).reshape(1,-1)   # [1,1]

        encode_pep = peptide_data(peptide)    # [10,21]

        encode_hla = hla_data(hla, dic_inventory, hla_type)   # [46,21]
        series.append((encode_pep, encode_hla, immuno))
    return series

def draw_ROC(y_true,y_pred):

    fpr,tpr,_ = roc_curve(y_true,y_pred,pos_label=1)
    area_mine = auc(fpr,tpr)
    fig = plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
            lw=lw, label='ROC curve (area = %0.2f)' % area_mine)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

def draw_PR(y_true,y_pred):
    precision,recall,cutoff = precision_recall_curve(y_true,y_pred,pos_label=1)
    area_PR = auc(recall,precision)
    baseline = np.sum(np.array(y_true) == 1) / len(y_true)

    plt.figure()
    lw = 2
    plt.plot(recall,precision, color='darkorange',
            lw=lw, label='PR curve (area = %0.2f)' % area_PR)
    plt.plot([0, 1], [baseline, baseline], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve example')
    plt.legend(loc="lower right")
    plt.show()
    return cutoff


def pull_peptide(dataset):
    result = np.empty([len(dataset),10,21,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla(dataset):
    result = np.empty([len(dataset),46,21,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result


def pull_label(dataset):
    result = np.empty([len(dataset),1])
    for i in range(len(dataset)):
        result[i,:] = dataset[i][2]
    return result


def aaindex(peptide,after_pca):

    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)   # [12,21]
    encoded = np.empty([len(peptide), 12])  # (seq_len,12)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded


def peptide_data_aaindex(peptide,after_pca):   # return numpy array [10,21,1]
    length = len(peptide)
    if length == 10:
        encode = aaindex(peptide,after_pca)
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]
        encode = aaindex(peptide,after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode


def hla_data_aaindex(hla, dic_inventory, hla_type,after_pca):    # return numpy array [46,21,1]
    dic = paratope_dic(hla)
    try:
        seq = dic[hla_type]
    except KeyError:
        hla_type = rescue_unknown_hla(hla_type, dic_inventory)
        seq = dic[hla_type]
    encode = aaindex(seq,after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode


def construct_aaindex(ori, hla, dic_inventory,after_pca):
    series = []
    for i in range(ori.shape[0]):
        peptide = ori['peptide'].iloc[i]
        hla_type = ori['HLA'].iloc[i]
        immuno = np.array(ori['immunogenecity'].iloc[i]).reshape(1,-1)   # [1,1]

        encode_pep = peptide_data_aaindex(peptide,after_pca)    # [10,12]

        encode_hla = hla_data_aaindex(hla, dic_inventory, hla_type,after_pca)   # [46,12]
        series.append((encode_pep, encode_hla, immuno))
    return series

def pull_peptide_aaindex(dataset):
    result = np.empty([len(dataset),10,12,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla_aaindex(dataset):
    result = np.empty([len(dataset),46,12,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result


def pull_label_aaindex(dataset):
    result = np.empty([len(dataset),1])
    for i in range(len(dataset)):
        result[i,:] = dataset[i][2]
    return result

if __name__ == '__main__':
    # instantiate two models
    ResLikeCNN = model()
    ResLikeCNN_index = model_aaindex()

    # load the weight
    ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')
    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')

    # training performance
    ## ResLikeCNN, blosum encoding
    ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    dataset = construct(ori, hla, dic_inventory)   # [ (10,21,1),(46,21,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    label = pull_label(dataset)
    result1 = ResLikeCNN.predict(x=[input1, input2])

    ## ResLikeCNN, aaindex encoding
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    dataset = construct_aaindex(ori,hla,dic_inventory,after_pca)
    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    label = pull_label_aaindex(dataset)
    result2 = ResLikeCNN_index.predict(x=[input1,input2])

    ## Ensemble
    result3 = np.mean(np.concatenate([result1,result2],axis=1),axis=1)

    ## store the data
    result_train = pd.DataFrame({'label':label[:,0],'blosum':result1[:,0],'aaindex':result2[:,0],'ensemble':result3})
    result_train.to_csv('paper/figure2/result_train.txt',sep='\t')

    # same for validation
    ori_val = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_val = construct(ori_val, hla, dic_inventory)
    input1_val = pull_peptide(dataset_val)
    input2_val = pull_hla(dataset_val)
    label_val = pull_label(dataset_val)
    result1 = ResLikeCNN.predict(x=[input1_val, input2_val])

    dataset_val = construct_aaindex(ori_val, hla, dic_inventory,after_pca)
    input1_val = pull_peptide_aaindex(dataset_val)
    input2_val = pull_hla_aaindex(dataset_val)
    label_val = pull_label_aaindex(dataset_val)
    result2 = ResLikeCNN_index.predict(x=[input1_val,input2_val])

    result3 = np.mean(np.concatenate([result1,result2],axis=1),axis=1)
    result_val = pd.DataFrame({'label':label_val[:,0],'blosum':result1[:,0],'aaindex':result2[:,0],'ensemble':result3})
    result_val.to_csv('paper/figure2/result_val.txt',sep='\t')

    # let's draw figure 2a
    training = pd.read_csv('paper/figure2/result_train.txt',sep='\t')
    arrayP = [training['blosum'],training['aaindex'],training['ensemble']]
    arrayT = [training['label'],training['label'],training['label']]
    draw_combined_ROC(arrayP,arrayT)
    draw_combined_PR(arrayP,arrayT)

    # same for validation
    val = pd.read_csv('paper/figure2/result_val.txt',sep='\t')
    arrayP = [val['blosum'],val['aaindex'],val['ensemble']]
    arrayT = [val['label'],val['label'],val['label']]
    draw_combined_ROC(arrayP,arrayT)
    draw_combined_PR(arrayP,arrayT)

    # external testing dataset, as I always been doing
    # refer to combine_plot.py

    # draw a barplot to further illustrate the performance boost
    fig,ax = plt.subplots()
    y1 = [0.50,0.67,0.75,0.78,0.82,0.83]
    x1 = np.arange(6)
    labels1 = ['iedb', 'logistic regression', 'deephlapan', 'PIER(Blosum Encoding)', 'PIER(AAindex Encoding)','PIER(Ensemble)']
    ax.bar(x1,y1,color=['red', 'green', 'blue', 'cyan', 'magenta', 'orange'])
    ax.set_xticks(x1)
    ax.set_xticklabels(labels1,rotation=60)

    # dbfilter910.txt
    ori_score = pd.read_csv('data/db_filter910.txt',sep='\t')
    dataset_score = construct_aaindex(ori_score,hla,dic_inventory,after_pca)
    input1_score = pull_peptide_aaindex(dataset_score)
    input2_score = pull_hla_aaindex(dataset_score)
    label_score = pull_label_aaindex(dataset_score)
    scoring2 = ResLikeCNN_index.predict(x=[input1_score,input2_score])

    scoring3 = np.mean(np.concatenate([scoring1,scoring2],axis=1),axis=1)


    # HLA-A*0201, 9mer
    test1 = pd.read_csv('/Users/ligk2e/Desktop/tmp/tmp1_9.txt',sep='\t')
    peptide_score = test1.iloc[:,0].tolist()
    hla_score = ['HLA-A*6801'] * len(peptide_score)
    immuno_score = ['0'] * len(peptide_score)
    ori_score = pd.DataFrame({'peptide':peptide_score,'HLA':hla_score,'immunogenecity':immuno_score})
    dataset_score = construct(ori_score,hla,dic_inventory,)
    input1_score = pull_peptide(dataset_score)
    input2_score = pull_hla(dataset_score)
    label_score = pull_label(dataset_score)
    exp1_1 = ResLikeCNN.predict(x=[input1_score,input2_score])

    dataset_score = construct_aaindex(ori_score,hla,dic_inventory,after_pca)
    input1_score = pull_peptide_aaindex(dataset_score)
    input2_score = pull_hla_aaindex(dataset_score)
    label_score = pull_label_aaindex(dataset_score)
    exp1_2 = ResLikeCNN_index.predict(x=[input1_score,input2_score])
    exp1_3 = np.mean(np.concatenate([exp1_1,exp1_2],axis=1),axis=1)
    pd.Series([1 if i > 0.5 else 0 for i in exp1_3]).to_csv('/Users/ligk2e/Desktop/tmp/exp7.txt', index=None)


    # after having the comparison of netMHCpan and PIER in 5 plus HLA types, let's draw barplot
    '''
    HLA-A0201   total:52, netMHCpan correctly indicated: 21, PIER correctly indicated: 24
    HLA-A0101   total:20, netMHCpan correctly indicated:0, PIER correctly indicated:13
    HLA-A1101   total:42, netMHCpan correctly indicated:5,PIER correctly indeicated:26
    HLA-A3101   total:23, netMHCpan correctly indicated:0, PIER correctly indicated: 21
    HLA-A6801   total:20  netMHCpan correctly indicated:3, PIER correctly indicated:17
    HLA-B0801   total:29 netMHCpan correctly indicated:2, PIER correctly indicated:2
    HLA-B5801   total:19, netMHCpan correctly indicated:6, PIER correctly indicated:5
    '''
    import itertools
    fig,axes = plt.subplots(2)
    # subplot1
    width=0.25
    y1_total = [20,52,42]
    y1_net = [0,21,5]
    y1_PIER = [13,24,26]
    x1_total = np.arange(3)
    x1_net = [i+width for i in x1_total]
    x1_PIER = [i+2*width for i in x1_total]
    axes[0].bar(x1_total,y1_total,color='green',width=width,edgecolor='white',label='total antigen counts')
    axes[0].bar(x1_net,y1_net,color='blue',width=width,edgecolor='white',label='netMHCpan correctly identified')
    axes[0].bar(x1_PIER,y1_PIER,color='orange',width=width,edgecolor='white',label='PIER correctly identified')
    axes[0].set_xticks(x1_net)
    axes[0].set_xticklabels(['HLA-A*0101','HLA-A*0201','HLA-A*1101'])
    axes[0].set_ylabel('Antigen Counts')
    axes[0].set_ylim((0,60))
    x1 = list(itertools.chain(x1_total,x1_net,x1_PIER))
    y1 = list(itertools.chain(y1_total,y1_net,y1_PIER))
    for i in range(len(x1)):
        axes[0].text(x1[i]-width/4,y1[i]+2,y1[i])

    # subplot2
    width=0.25
    y2_total = [23,20,29]
    y2_net = [0,3,2]
    y2_PIER = [21,17,2]
    x2_total = np.arange(3)
    x2_net = [i+width for i in x2_total]
    x2_PIER = [i+2*width for i in x2_total]
    axes[1].bar(x2_total,y2_total,color='green',width=width,edgecolor='white',label='total antigen counts')
    axes[1].bar(x2_net,y2_net,color='blue',width=width,edgecolor='white',label='netMHCpan correctly identified')
    axes[1].bar(x2_PIER,y2_PIER,color='orange',width=width,edgecolor='white',label='PIER correctly identified')
    axes[1].set_xticks(x2_net)
    axes[1].set_xticklabels(['HLA-A*3101','HLA-A*6801','HLA-B*0801'])
    axes[1].set_ylabel('Antigen Counts')
    axes[1].set_ylim((0,37))
    x2 = list(itertools.chain(x2_total,x2_net,x2_PIER))
    y2 = list(itertools.chain(y2_total,y2_net,y2_PIER))
    for i in range(len(x2)):
        axes[1].text(x2[i]-width/4,y2[i]+2,y2[i])

    # add legend
    axes[0].legend(loc='upper left',bbox_to_anchor=(1,1))
    fig.tight_layout()
    plt.show()
    plt.savefig('paper/figure3/reveal.pdf')

    ## if we do transfer learning, the performance doesn't improve
    def draw_combined_ROC(arrayP, arrayT):
        fig = plt.figure()
        colormap = ['red', 'green', 'blue', 'orange']
        legend = ['FineTune_last_one_layer','Finetune_last_two_layer','Finetune_last_three_layer','PIER(Ensemble)']
        # colormap = ['cyan','magenta','orange']
        # legend = ['PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)']
        for i in range(len(arrayP)):
            fpr, tpr, _ = roc_curve(arrayT[i], arrayP[i], pos_label=1)
            area = auc(fpr, tpr)
            lw = 2
            plt.plot(fpr, tpr, color=colormap[i],
                     lw=lw, label='{0} (area = {1:0.2f})'.format(legend[i], area))
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc="lower right")
        plt.show()
        plt.savefig('paper/figure3/transfer_learning_ROC.pdf')


    def draw_combined_PR(arrayP, arrayT):
        fig = plt.figure()
        colormap = ['red', 'green', 'blue', 'orange']
        legend = ['FineTune_last_one_layer','Finetune_last_two_layer','Finetune_last_three_layer','PIER(Ensemble)']
        # colormap = ['cyan','magenta','orange']
        # legend = ['PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)']
        for i in range(len(arrayP)):
            precision, recall, _ = precision_recall_curve(arrayT[i], arrayP[i], pos_label=1)
            area = auc(recall, precision)
            lw = 2
            plt.plot(recall, precision, color=colormap[i],
                     lw=lw, label='{0} (area = {1:0.2f})'.format(legend[i], area))
        plt.plot([0, 1], [0.12, 0.12], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(loc="best")
        plt.show()
        plt.savefig('paper/figure3/transfer_learning_pr.pdf')

    # make sure ResLikeCNN has been initiated
    ResLikeCNN = model()
    ResLikeCNN.load_weights('paper/figure3/transfer_learning/fine_tune/freeze0123/')
    ResLikeCNN.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    result_1 = ResLikeCNN.predict(x=[input1_test,input2_test])   # Finetune_last_one_layer

    ResLikeCNN = model()
    ResLikeCNN.load_weights('paper/figure3/transfer_learning/fine_tune/freeze012/')
    ResLikeCNN.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    result_2 = ResLikeCNN.predict(x=[input1_test,input2_test])   # Finetune_last_two_layer


    ResLikeCNN = model()
    ResLikeCNN.load_weights('paper/figure3/transfer_learning/fine_tune/freeze01/')
    ResLikeCNN.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    result_3 = ResLikeCNN.predict(x=[input1_test,input2_test])   # Finetune_last_three_layer

    # Ensemble model load, by the way save your ensemble result
    ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    ResLikeCNN_index = model_aaindex()
    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')
    ResLikeCNN_index.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct_aaindex(ori_test, hla, dic_inventory,after_pca)
    input1_test = pull_peptide_aaindex(dataset_test)
    input2_test = pull_hla_aaindex(dataset_test)
    label_test = pull_label_aaindex(dataset_test)
    result_aaindex = ResLikeCNN_index.predict(x=[input1_test,input2_test])

    ResLikeCNN = model()
    ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')
    ResLikeCNN_index.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    result_blosum = ResLikeCNN.predict(x=[input1_test,input2_test])

    result_ensemble = np.mean(np.concatenate([result_aaindex,result_blosum],axis=1),axis=1)
    result_df = pd.DataFrame({'label':label_test[:,0],'aaindex':result_aaindex[:,0],'blosum':result_blosum[:,0],'ensemble':result_ensemble})
    result_df.to_csv('paper/model/ensemble_result.txt',sep='\t',index=None)

    transfer_df = pd.DataFrame({'FineTune_last_one':result_1[:,0],'FineTune_last_two':result_2[:,0],'FineTune_last_three':result_3[:,0]})
    transfer_df.to_csv('paper/figure3/transfer_df.txt',sep='\t',index=None)

    arrayP = [result_1,result_2,result_3,result_ensemble]
    arrayT = [label_test,label_test,label_test,label_test]
    draw_combined_ROC(arrayP,arrayT)
    draw_combined_PR(arrayP,arrayT)










