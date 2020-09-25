'''
Accept peptide and MHC input,
compute immunogenicity


'''

import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler,MinMaxScaler,RobustScaler
import collections
import re



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

def computing_s(peptide,mhc):
    # print(peptide)
    # print(mhc)
    # print(type(peptide))
    # print(type(mhc))
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    ResLikeCNN_index = model_aaindex()
    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')

    peptide_score = [peptide]
    hla_score = [mhc]
    immuno_score = ['0']
    ori_score = pd.DataFrame({'peptide':peptide_score,'HLA':hla_score,'immunogenecity':immuno_score})
    dataset_score = construct_aaindex(ori_score,hla,dic_inventory,after_pca)
    input1_score = pull_peptide_aaindex(dataset_score)
    input2_score = pull_hla_aaindex(dataset_score)
    label_score = pull_label_aaindex(dataset_score)
    scoring = ResLikeCNN_index.predict(x=[input1_score,input2_score])
    return float(scoring)

def computing_m(peptide,mhc):    # multiple MHC query
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    ResLikeCNN_index = model_aaindex()
    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')

    hla_score = ['HLA-A*0203','HLA-A*0205','HLA-A*0206','HLA-A*0207','HLA-A*3001','HLA-A*6801','HLA-A*9235','HLA-A*9253','HLA-B*0602','HLA-B*1557','HLA-B*1801','HLA-B*3505','HLA-B*4001','HLA-B*4002','HLA-B*4044','HLA-B*4102','HLA-B*4104','HLA-B*4202','HLA-B*4601','HLA-B*5706','HLA-B*5712','HLA-B*8103','HLA-B*9234','HLA-C*0102','HLA-C*0304','HLA-C*0401','HLA-C*0517','HLA-C*0602','HLA-C*0756','HLA-C*1604','HLA-A*0101','HLA-A*0201','HLA-A*0224','HLA-A*0301','HLA-A*0362','HLA-A*1101','HLA-A*2301','HLA-A*2402','HLA-A*3003','HLA-A*9234','HLA-B*0702','HLA-B*0801','HLA-B*1402','HLA-B*1501','HLA-B*2703','HLA-B*2705','HLA-B*2709','HLA-B*2713','HLA-B*3501','HLA-B*3508','HLA-B*3901','HLA-B*4201','HLA-B*4402','HLA-B*4403','HLA-B*4405','HLA-B*5101','HLA-B*5301','HLA-B*5701','HLA-B*5703','HLA-B*5801','HLA-B*8102','HLA-C*1510']
    peptide_score = [peptide] * len(hla_score)
    immuno_score = ['0'] * len(hla_score)
    ori_score = pd.DataFrame({'peptide':peptide_score,'HLA':hla_score,'immunogenecity':immuno_score})
    dataset_score = construct_aaindex(ori_score,hla,dic_inventory,after_pca)
    input1_score = pull_peptide_aaindex(dataset_score)
    input2_score = pull_hla_aaindex(dataset_score)
    label_score = pull_label_aaindex(dataset_score)
    scoring = ResLikeCNN_index.predict(x=[input1_score,input2_score])
    ori_score['immunogenecity'] = scoring
    ori_score.sort_values(by=['immunogenecity'],ascending=False,inplace=True)
    top5 = ori_score.iloc[0:5]
    
    p = top5['peptide'].tolist()
    m = top5['HLA'].tolist()
    i = [item for item in top5['immunogenecity']]
    return p,m,i


def hla_convert(hla):
    hla = hla.replace('*','')
    return hla

def svg_path(hla):
    path = ["/static/{0}_positive_9.png".format(hla),"/static/{0}_negative_9.png".format(hla),"/static/{0}_positive_10.png".format(hla),"/static/{0}_negative_10.png".format(hla)]
    