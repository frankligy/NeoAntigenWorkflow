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

def add_GAN_sample(dataset,pseudo_p_total,pseudo_MHC_total):
    length = pseudo_p_total.shape[0]
    pseudo_label_total = np.ones([length,1])   # 12800,1
    pseudo_p_total = [np.transpose(item,(1,2,0)).astype(np.float64) for item in pseudo_p_total]  # 12800,1,10,12
    pseudo_MHC_total = [np.transpose(item,(1,2,0)).astype(np.float64) for item in pseudo_MHC_total]  # 12800,1,46,12
    for i in range(length):
        tup = (pseudo_p_total[i],pseudo_MHC_total[i],pseudo_label_total[i])
        dataset.append(tup)
    return dataset




if __name__ == '__main__':
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    dataset = construct_aaindex(ori, hla, dic_inventory,after_pca)   # [ (10,12,1),(46,12,1),(1,1)   ]
    # load the GAN generated positive samples
    import pickle
    with open('/Users/ligk2e/Desktop/tmp/pseudo_p_total.p','rb') as f1:
        pseudo_p_total = pickle.load(f1)
    with open('/Users/ligk2e/Desktop/tmp/pseudo_MHC_total.p','rb') as f2:
        pseudo_MHC_total = pickle.load(f2)
    dataset = add_GAN_sample(dataset,pseudo_p_total,pseudo_MHC_total)

    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    label = pull_label_aaindex(dataset)



    ori_val = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_val = construct_aaindex(ori_val, hla, dic_inventory,after_pca)
    input1_val = pull_peptide_aaindex(dataset_val)
    input2_val = pull_hla_aaindex(dataset_val)
    label_val = pull_label_aaindex(dataset_val)

    ResLikeCNN_index = model_aaindex()

    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN')

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
        epochs=7,
        class_weight = {0:0.2,1:0.8}   # I have 20% positive and 80% negative in my training data
    )

    # now let's test in external dataset
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',
                           sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct_aaindex(ori_test, hla, dic_inventory,after_pca)
    input1_test = pull_peptide_aaindex(dataset_test)
    input2_test = pull_hla_aaindex(dataset_test)
    label_test = pull_label_aaindex(dataset_test)
    # seperateCNNmodel.evaluate(x=[input1_test,input2_test],y=label_test,batch_size=512)
    result = ResLikeCNN_index.predict(x=[input1_test, input2_test])
    hard = [1 if i > 0.5 else 0 for i in result]
    confusion_matrix(label_test, hard)
    f1_score(label_test, hard)
    accuracy_score(label_test, hard)
    draw_ROC(label_test, result)
    draw_PR(label_test, result)


    # external neoantigen dataset
    ext_test = pd.read_csv('data/mannual_cancer_testing_fiter910.txt',sep='\t')
    dataset_ext = construct_aaindex(ext_test, hla, dic_inventory,after_pca)
    input1_ext = pull_peptide_aaindex(dataset_ext)
    input2_ext = pull_hla_aaindex(dataset_ext)
    label_ext = pull_label_aaindex(dataset_ext)
    result = ResLikeCNN_index.predict(x=[input1_ext, input2_ext])
    hard = [1 if i > 0.5 else 0 for i in result]
    confusion_matrix(label_ext, hard)
    f1_score(label_ext, hard)
    accuracy_score(label_ext, hard)
    draw_ROC(label_ext, result)
    draw_PR(label_ext, result)

    ResLikeCNN_index.save_weights('aaindex12_encoding_ReslikeCNN_reproduce/')


    ResLikeCNN_index.save_weights('GAN_augmented/')
    ResLikeCNN_index.save_weights('aaindex12_encoding_ReslikeCNN')
