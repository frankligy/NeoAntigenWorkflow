import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import matplotlib.pyplot as plt
import numpy as np
from Bio.SubsMat import MatrixInfo
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler,MinMaxScaler,RobustScaler
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import collections

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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


class CNN_peptide_aaindex(layers.Layer):
    def __init__(self):
        super(CNN_peptide_aaindex,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(3,17),strides=(1,1))
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
        self.conv = layers.Conv2D(filters=16,kernel_size=(15,17),strides=(1,1)) # (32,1,16)
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
        x1 = keras.Input(shape=(10,9,1))
        x2 = keras.Input(shape=(46,9,1))
        return keras.Model(inputs=[x1,x2],outputs=self.call([x1,x2]))

def aaindex(peptide,after_pca):

    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)   # [12,21]
    encoded = np.empty([len(peptide), 17])  # (seq_len,12)
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
    result = np.empty([len(dataset),10,17,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla_aaindex(dataset):
    result = np.empty([len(dataset),46,17,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result


def pull_label_aaindex(dataset):
    result = np.empty([len(dataset),1])
    for i in range(len(dataset)):
        result[i,:] = dataset[i][2]
    return result



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


def dict_to_matrix():
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    dic = MatrixInfo.blosum60
    matrix = np.zeros([21,21])

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            try:
                matrix[i, j] = dic[(amino[i], amino[j])]
            except KeyError:
                try:
                    matrix[i, j] = dic[(amino[j], amino[i])]
                except:
                    matrix[i, j] = -1
    return matrix


def blosum(peptide):
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = dict_to_matrix()
    encoded = np.empty([len(peptide), 21])  # (seq_len,21)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded

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


def add_X(array):
    me = np.mean(array)
    array = np.append(array, me)
    return array


def read_index(path):
    with open(path, 'r') as f:
        data = f.readlines()
        array = []
        for line in data:
            line = line.lstrip(' ').rstrip('\n')
            line = re.sub(' +', ' ', line)

            items = line.split(' ')
            items = [float(i) for i in items]
            array.extend(items)
        array = np.array(array)
        array = add_X(array)
        Index = collections.namedtuple('Index',
                                       ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                                        'T', 'W', 'Y', 'V', 'X'])
        I = Index._make(array)
    return I, array  # namedtuple


def read_all_indices():
    table = np.empty([21, 566])
    for i in range(566):
        if len(str(i)) == 1:
            ii = '00' + str(i)
        elif len(str(i)) == 2:
            ii = '0' + str(i)
        else:
            ii = str(i)

        NA_list_str = ['472', '473', '474', '475', '476', '477', '478', '479', '480', '481', '520', '523', '524']
        NA_list_int = [int(i) for i in NA_list_str]
        if ii in NA_list_str: continue

        path = '/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/AAindex1/index{0}.txt'.format(ii)

        _, array = read_index(path)

        table[:, i] = array
    table = np.delete(table, NA_list_int, 1)
    return table


def scaling(table):  # scale the features
    table_scaled = RobustScaler().fit_transform(table)
    return table_scaled


def wrapper_read_scaling():
    table = read_all_indices()
    table_scaled = scaling(table)
    return table_scaled


def pca_get_components(result):
    pca= PCA()
    pca.fit(result)
    result = pca.explained_variance_ratio_
    sum_ = 0
    for index,var in enumerate(result):
        sum_ += var
        if sum_ > 0.95:
            return index    # 25 components



def pca_apply_reduction(result):   # if 95%, 12 PCs, if 99%, 17 PCs, if 90%,9 PCs
    pca = PCA(n_components=12)  # or strictly speaking ,should be 26, since python is 0-index
    new = pca.fit_transform(result)
    return new

if __name__ == '__main__':
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    dataset = construct_aaindex(ori, hla, dic_inventory,after_pca)   # [ (10,12,1),(46,12,1),(1,1)   ]
    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    label = pull_label_aaindex(dataset)

    ori_val = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_val = construct_aaindex(ori_val, hla, dic_inventory,after_pca)
    input1_val = pull_peptide_aaindex(dataset_val)
    input2_val = pull_hla_aaindex(dataset_val)
    label_val = pull_label_aaindex(dataset_val)

    ResLikeCNN_index = model_aaindex()

    ResLikeCNN_index.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    history = ResLikeCNN_index.fit(
        x=[input1,input2],   # feed a list into
        y=label,
        validation_data = ([input1_val,input2_val],label_val),
        batch_size=512,
        epochs=9,
        class_weight = {0:0.2,1:0.8}   # I have 20% positive and 80% negative in my training data
    )

    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',
                           sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct_aaindex(ori_test, hla, dic_inventory,after_pca)
    input1_test = pull_peptide_aaindex(dataset_test)
    input2_test = pull_hla_aaindex(dataset_test)
    label_test = pull_label_aaindex(dataset_test)
    # seperateCNNmodel.evaluate(x=[input1_test,input2_test],y=label_test,batch_size=512)
    result17 = ResLikeCNN_index.predict(x=[input1_test, input2_test])
    hard = [1 if i > 0.5 else 0 for i in result17]
    confusion_matrix(label_test, hard)
    f1_score(label_test, hard)
    accuracy_score(label_test, hard)
    draw_ROC(label_test, result17)
    draw_PR(label_test, result17)

    ResLikeCNN_index.save_weights('paper/figure5/different_encoding/aaindex17/')
    result_aaindex = pd.DataFrame({'aaindex9':result9[:,0],'aaindex17':result17[:,0]})
    result_aaindex.to_csv('paper/figure5/result_aaindex.txt',sep='\t',index=None)

    ResLikeCNN = model()

    dataset = construct(ori, hla, dic_inventory)   # [ (10,12,1),(46,12,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    label = pull_label(dataset)

    dataset_val = construct(ori_val, hla, dic_inventory)
    input1_val = pull_peptide(dataset_val)
    input2_val = pull_hla(dataset_val)
    label_val = pull_label(dataset_val)

    ResLikeCNN.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    history = ResLikeCNN.fit(
        x=[input1, input2],  # feed a list into
        y=label,
        validation_data=([input1_val, input2_val], label_val),
        batch_size=512,
        epochs=9,
        class_weight={0: 0.2, 1: 0.8}  # I have 20% positive and 80% negative in my training data
    )

    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)

    result_blosum30 = ResLikeCNN.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result_blosum30]
    confusion_matrix(label_test, hard)
    f1_score(label_test, hard)
    accuracy_score(label_test, hard)
    draw_ROC(label_test, result_blosum30)
    draw_PR(label_test, result_blosum30)

    ResLikeCNN.save_weights('paper/figure5/different_encoding/blosum30/')

    # test blosum60
    ResLikeCNN = model()
    result_blosum60 = ResLikeCNN.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result_blosum60]
    confusion_matrix(label_test, hard)
    f1_score(label_test, hard)
    accuracy_score(label_test, hard)
    draw_ROC(label_test, result_blosum60)
    draw_PR(label_test, result_blosum60)

    ResLikeCNN.save_weights('paper/figure5/different_encoding/blosum60/')

    result_blosum = pd.DataFrame({'blosum30':result_blosum30[:,0],'blosum60':result_blosum60[:,0]})
    result_blosum.to_csv('paper/figure5/result_blosum.txt',sep='\t',index=None)

    '''
    Let's start to draw figures
    '''

    result_blosum = pd.read_csv('paper/figure5/result_blosum.txt',sep='\t')
    result_aaindex = pd.read_csv('paper/figure5/result_aaindex.txt',sep='\t')
    result_ori = pd.read_csv('paper/model/ensemble_result.txt',sep='\t')


    def draw_combined_ROC(arrayP, arrayT):
        fig = plt.figure()
        colormap = ['red', 'green','magenta', 'orange']
        legend = ['AAindex(PC=9,variance=90%)','AAindex(PC=17,variance=99%)','AAindex(PIER,PC=12,variance=95%)',
                  'PIER(Ensemble)']
        for i in range(len(arrayP)):
            fpr, tpr, _ = roc_curve(arrayT[i], arrayP[i], pos_label=1)
            area = auc(fpr, tpr)
            lw = 2
            plt.plot(fpr, tpr, color=colormap[i],
                     lw=lw, label='{0}'.format(legend[i]))
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc="lower right")
        plt.show()
        plt.savefig('paper/figure5/aaindex_encoding_ROC.pdf')


    def draw_combined_PR(arrayP, arrayT):
        fig = plt.figure()
        colormap = ['red', 'green','magenta', 'orange']
        legend = ['AAindex(PC=9,variance=90%)','AAindex(PC=17,variance=99%)','AAindex(PIER,PC=12,variance=95%)',
                  'PIER(Ensemble)']
        for i in range(len(arrayP)):
            precision, recall, _ = precision_recall_curve(arrayT[i], arrayP[i], pos_label=1)
            area = auc(recall, precision)
            lw = 2
            plt.plot(recall, precision, color=colormap[i],
                     lw=lw, label='{0}'.format(legend[i]))
        plt.plot([0, 1], [0.12, 0.12], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(loc="best")
        plt.show()
        plt.savefig('paper/figure5/AAindex_encoding_PR.pdf')
    # panelA aaindex comparison
    arrayP = [result_aaindex['aaindex9'],result_aaindex['aaindex17'],result_ori['aaindex'],result_ori['ensemble']]
    arrayT = [result_ori['label'],result_ori['label'],result_ori['label'],result_ori['label']]
    draw_combined_PR(arrayP,arrayT)
    draw_combined_ROC(arrayP,arrayT)

    ##barplot
    '''
    aaindex9,aaindex17,aaindex12,emsemble
        ROC: 0.77, 0.79,0.82,0.83
        PR: 0.27,0.25,0.37,0.41
    '''

    fig,ax = plt.subplots()
    y = [0.77,0.79,0.82,0.83]
    x = np.arange(4)
    width=0.5
    ax.bar(x[0],y[0],color='red',width=width)
    ax.bar(x[1],y[1],color='green',width=width)
    ax.bar(x[2],y[2],color='magenta',width=width)
    ax.bar(x[3],y[3],color='orange',width=width)
    ax.set_ylim(0,1.0)
    ax.set_ylabel('AUROC')
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(['aaindex(PC=9)','aaindex(PC=17)','aaindex(PIER,PC=12)','PIER(Ensemble)'],rotation=20,fontsize=7)
    for i in range(4):
        ax.text(x[i]-width/4,y[i]+0.03,y[i])
    plt.show()
    plt.savefig('paper/figure5/aaindex_barplot_AUROC.pdf')

    fig,ax = plt.subplots()
    y = [0.27,0.25,0.37,0.41]
    x = np.arange(4)
    width=0.5
    ax.bar(x[0],y[0],color='red',width=width)
    ax.bar(x[1],y[1],color='green',width=width)
    ax.bar(x[2],y[2],color='magenta',width=width)
    ax.bar(x[3],y[3],color='orange',width=width)
    ax.set_ylim(0,0.51)
    ax.set_ylabel('AUPR')
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(['aaindex(PC=9)','aaindex(PC=17)','aaindex(PIER,PC=12)','PIER(Ensemble)'],rotation=20,fontsize=7)
    for i in range(4):
        ax.text(x[i]-width/4,y[i]+0.02,y[i])
    plt.show()
    plt.savefig('paper/figure5/aaindex_barplot_AUPR.pdf')


    # same for blosum
    def draw_combined_ROC(arrayP, arrayT):
        fig = plt.figure()
        colormap = ['red', 'green','cyan', 'orange']
        legend = ['BLOSUM30','BLOSUM60','BLOSUM50(PIER)',
                  'PIER(Ensemble)']
        for i in range(len(arrayP)):
            fpr, tpr, _ = roc_curve(arrayT[i], arrayP[i], pos_label=1)
            area = auc(fpr, tpr)
            lw = 2
            plt.plot(fpr, tpr, color=colormap[i],
                     lw=lw, label='{0}'.format(legend[i]))
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc="lower right")
        plt.show()
        plt.savefig('paper/figure5/blosum_encoding_ROC.pdf')


    def draw_combined_PR(arrayP, arrayT):
        fig = plt.figure()
        colormap = ['red', 'green','cyan', 'orange']
        legend = ['BLOSUM30','BLOSUM60','BLOSUM50(PIER)',
                  'PIER(Ensemble)']
        for i in range(len(arrayP)):
            precision, recall, _ = precision_recall_curve(arrayT[i], arrayP[i], pos_label=1)
            area = auc(recall, precision)
            lw = 2
            plt.plot(recall, precision, color=colormap[i],
                     lw=lw, label='{0}'.format(legend[i]))
        plt.plot([0, 1], [0.12, 0.12], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(loc="best")
        plt.show()
        plt.savefig('paper/figure5/blosum_encoding_PR.pdf')

    arrayP = [result_blosum['blosum30'],result_blosum['blosum60'],result_ori['blosum'],result_ori['ensemble']]
    arrayT = [result_ori['label'],result_ori['label'],result_ori['label'],result_ori['label']]
    draw_combined_PR(arrayP,arrayT)
    draw_combined_ROC(arrayP,arrayT)

    '''
    blosum30,blosum60,blosum50,ensemble
        ROC: 0.79,0.75,0.78,0.83
        PR: 0.29,0.27,0.33,0.41
    '''

    fig,ax = plt.subplots()
    y = [0.79,0.75,0.78,0.83]
    x = np.arange(4)
    width=0.5
    ax.bar(x[0],y[0],color='red',width=width)
    ax.bar(x[1],y[1],color='green',width=width)
    ax.bar(x[2],y[2],color='magenta',width=width)
    ax.bar(x[3],y[3],color='orange',width=width)
    ax.set_ylim(0,1.0)
    ax.set_ylabel('AUROC')
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(['BLOSUM30','BLOSUM60','BLOSUM50(PIER)','PIER(Ensemble)'],rotation=20,fontsize=7)
    for i in range(4):
        ax.text(x[i]-width/4,y[i]+0.03,y[i])
    plt.show()
    plt.savefig('paper/figure5/blosum_barplot_AUROC.pdf')

    fig,ax = plt.subplots()
    y = [0.29,0.27,0.33,0.41]
    x = np.arange(4)
    width=0.5
    ax.bar(x[0],y[0],color='red',width=width)
    ax.bar(x[1],y[1],color='green',width=width)
    ax.bar(x[2],y[2],color='magenta',width=width)
    ax.bar(x[3],y[3],color='orange',width=width)
    ax.set_ylim(0,0.51)
    ax.set_ylabel('AUPR')
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(['BLOSUM30','BLOSUM60','BLOSUM50(PIER)','PIER(Ensemble)'],rotation=20,fontsize=7)
    for i in range(4):
        ax.text(x[i]-width/4,y[i]+0.02,y[i])
    plt.show()
    plt.savefig('paper/figure5/blosum_barplot_AUPR.pdf')

    ### Finally, let's proof ensemble method makes sense
    import scipy.stats as sc
    p = sc.pearsonr(result_ori['blosum'].values,result_ori['aaindex'].values)
    fig,ax = plt.subplots()
    ax.scatter(result_ori['blosum'],result_ori['aaindex'],color='green',s=10)
    ax.set_ylabel('AAindex(PC=12) Encoding prediction score')
    ax.set_xlabel('BLOSUM50 Encoding prediction score')
    ax.legend(['Pearson Correlation: {0:.2f}'.format(p[0])],loc='best')
    plt.savefig('paper/figure5/correlation.pdf')

    m,b = np.polyfit(result_ori['blosum'].values,result_ori['aaindex'].values,deg=1)
    ax.plot(result_ori['blosum'],m*result_ori['blosum']+b,color='red')



















