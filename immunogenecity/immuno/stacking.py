'''
A stacking Ensemble model integrating random forest and ResNet,
Encoding strategies:
    1. AAindex1
    2. BLOSUM50
'''

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd
import os
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
from sklearn.ensemble import RandomForestClassifier
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler,MinMaxScaler,RobustScaler
import collections
from Bio.SubsMat import MatrixInfo
import matplotlib.pyplot as plt

# block1: pre-processing AAindex encoding
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

# block2: BLOSUM50 and hla preprossing

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


def dict_to_matrix():
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    dic = MatrixInfo.blosum50
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
        encode = blosum(peptide)
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

# block3: aaindex ResNet
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

def pull_peptide_aaindex(dataset):
    result = np.empty([len(dataset), 10, 12, 1])
    for i in range(len(dataset)):
        result[i, :, :, :] = dataset[i][0]
    return result

def pull_hla_aaindex(dataset):
    result = np.empty([len(dataset), 46, 12, 1])
    for i in range(len(dataset)):
        result[i, :, :, :] = dataset[i][1]
    return result

def pull_label_aaindex(dataset):
    result = np.empty([len(dataset), 1])
    for i in range(len(dataset)):
        result[i, :] = dataset[i][2]
    return result

# block4:
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
    precision,recall,_ = precision_recall_curve(y_true,y_pred,pos_label=1)
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

if __name__ == '__main__':
    # training data and validation data as a whole
    immuno_training = pd.read_csv('data/shuffle_training_test.txt',sep='\t')
    immuno_val = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')
    immuno = pd.concat([immuno_training,immuno_val])

    # hla
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    # for aaindex
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]

    # 10-fold cross-training
    idx = np.arange(immuno.shape[0])  # all index
    np.random.shuffle(idx)   # shuffle
    folds_idx = np.array_split(idx,10)
    folds_df = []
    for fold in folds_idx:
        sub = immuno.iloc[fold].set_index(pd.Index(np.arange(len(fold))))
        folds_df.append(sub)

    # First train random forest
    ## aaindex
    predictions = []
    for i in range(len(folds_df)):
        # ith fold will be prediction, other 9 will be training
        training = [data for m,data in enumerate(folds_df) if m != i]
        ori = pd.concat(training)
        dataset = construct_aaindex(ori,hla,dic_inventory,after_pca)
        X = np.empty((len(dataset), 12 * 56))  # 28581
        Y = ori['immunogenecity'].values

        for j, (x, y,_) in enumerate(dataset):
            x = x.reshape(-1)  # 10*12*1 ---> 120
            y = y.reshape(-1)  # 46*12*1 ---> 552
            X[j,:] = np.concatenate([x,y])   # 672
        clf = RandomForestClassifier()
        clf.fit(X, Y)

        ori_test = folds_df[i]
        testing_dataset = construct_aaindex(ori_test, hla, dic_inventory, after_pca)
        X_test = np.empty((len(testing_dataset), 12 * 56))
        Y_test = ori_test['immunogenecity'].values

        for i, (x, y, _) in enumerate(testing_dataset):
            x = x.reshape(-1)
            y = y.reshape(-1)
            X_test[i, :] = np.concatenate([x, y])
        y_pred = clf.predict_proba(X_test)
        ori_test['rf_aaindex'] = y_pred[:,1]
        predictions.append(ori_test)

    total_predictions = pd.concat(predictions)
    total_predictions.to_csv('total_prediction_rf_aaindex.txt',sep='\t',index=None)

    # BLOSUM
    predictions = []
    for i in range(len(folds_df)):
        # ith fold will be prediction, other 9 will be training
        training = [data for m,data in enumerate(folds_df) if m != i]
        ori = pd.concat(training)
        dataset = construct(ori,hla,dic_inventory)
        X = np.empty((len(dataset), 21 * 56))
        Y = ori['immunogenecity'].values

        for j, (x, y,_) in enumerate(dataset):
            x = x.reshape(-1)
            y = y.reshape(-1)
            X[j,:] = np.concatenate([x,y])
        clf = RandomForestClassifier()
        clf.fit(X, Y)

        ori_test = folds_df[i]
        testing_dataset = construct(ori_test, hla, dic_inventory)
        X_test = np.empty((len(testing_dataset), 21 * 56))
        Y_test = ori_test['immunogenecity'].values

        for i, (x, y, _) in enumerate(testing_dataset):
            x = x.reshape(-1)
            y = y.reshape(-1)
            X_test[i, :] = np.concatenate([x, y])
        y_pred = clf.predict_proba(X_test)
        ori_test['rf_blosum'] = y_pred[:,1]
        predictions.append(ori_test)

    total_predictions = pd.concat(predictions)
    total_predictions.to_csv('total_prediction_rf_blosum.txt',sep='\t',index=None)

    # aaindex_ResNet (load the model and the result)
    # blosum ResNet (load the model and the result)




######################################################################################################
    immuno_training = pd.read_csv('data/shuffle_training_test.txt',sep='\t')
    immuno_val = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')

    # hla
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    # for aaindex
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]

    # RF_aaindex
    dataset = construct_aaindex(immuno_training, hla, dic_inventory, after_pca)
    X = np.empty((len(dataset), 12 * 56))  # 28581
    Y = immuno_training['immunogenecity'].values

    for j, (x, y, _) in enumerate(dataset):
        x = x.reshape(-1)  # 10*12*1 ---> 120
        y = y.reshape(-1)  # 46*12*1 ---> 552
        X[j, :] = np.concatenate([x, y])  # 672
    clf = RandomForestClassifier()
    clf.fit(X, Y)

    from joblib import dump,load
    dump(clf,'RF_aaindex.joblib')

    ori_test = immuno_val
    testing_dataset = construct_aaindex(ori_test, hla, dic_inventory, after_pca)
    X_test = np.empty((len(testing_dataset), 12 * 56))
    Y_test = ori_test['immunogenecity'].values

    for i, (x, y, _) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i, :] = np.concatenate([x, y])
    y_pred_rf_aaindex = clf.predict_proba(X_test)

    # RF_blosum
    dataset = construct(immuno_training, hla, dic_inventory)
    X = np.empty((len(dataset), 21 * 56))  # 28581
    Y = immuno_training['immunogenecity'].values

    for j, (x, y, _) in enumerate(dataset):
        x = x.reshape(-1)  # 10*12*1 ---> 120
        y = y.reshape(-1)  # 46*12*1 ---> 552
        X[j, :] = np.concatenate([x, y])  # 672
    clf = RandomForestClassifier()
    clf.fit(X, Y)

    from joblib import dump,load
    dump(clf,'RF_blosum.joblib')

    ori_test = immuno_val
    testing_dataset = construct(ori_test, hla, dic_inventory)
    X_test = np.empty((len(testing_dataset), 21 * 56))
    Y_test = ori_test['immunogenecity'].values

    for i, (x, y, _) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i, :] = np.concatenate([x, y])
    y_pred_rf_blosum = clf.predict_proba(X_test)

    df_rf_validation_set = pd.DataFrame({'label':Y_test,'aaindex':y_pred_rf_aaindex[:,1],'blosum':y_pred_rf_blosum[:,1]})
    df_rf_validation_set.to_csv('df_rf_validation_set.txt',sep='\t',index=None)

    # Resnet aaindex
    ResLikeCNN_index = model_aaindex()
    ResLikeCNN_index.load_weights('aaindex12_encoding_ResLikeCNN_reproduce/')
    dataset = construct_aaindex(immuno_val, hla, dic_inventory,after_pca)   # [ (10,12,1),(46,12,1),(1,1)   ]
    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    label = pull_label_aaindex(dataset)
    y_pred_resnet_aaindex = ResLikeCNN_index.predict(x=[input1,input2])

    # ResNet blosum
    ResLikeCNN = model()
    ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')
    dataset = construct(immuno_val, hla, dic_inventory)   # [ (10,12,1),(46,12,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    label = pull_label(dataset)
    y_pred_resnet_blosum = ResLikeCNN.predict(x=[input1,input2])

    df_resnet_validation_set = pd.DataFrame({'label':Y_test,'aaindex':y_pred_resnet_aaindex[:,0],'blosum':y_pred_resnet_blosum[:,0]})
    df_resnet_validation_set.to_csv('df_resnet_validation_set.txt',sep='\t',index=None)

    ### Finally meta-classfier
    df_rf_validation_set = pd.read_csv('df_rf_validation_set.txt',sep='\t')
    df_resnet_validation_set = pd.read_csv('df_resnet_validation_set.txt',sep='\t')

    data = np.empty([df_rf_validation_set.shape[0],5])
    for i in range(data.shape[0]):
        data[i,0] = df_rf_validation_set.iloc[i]['aaindex']
        data[i,1] = df_rf_validation_set.iloc[i]['blosum']
        data[i,2] = df_resnet_validation_set.iloc[i]['aaindex']
        data[i,3] = df_resnet_validation_set.iloc[i]['blosum']
        data[i,4] = df_resnet_validation_set.iloc[i]['label']

    # perform a logistic regression
    X = data[:,0:3]
    y = data[:,4]
    from sklearn.linear_model import LogisticRegression
    clf = LogisticRegression()
    clf.fit(X,y)

    from joblib import dump,load
    dump(clf,'logisticRegression.joblib')

    result = clf.predict_proba(X)
    draw_ROC(y,result[:,1])   # 0.93    # 0.88
    draw_PR(y,result[:,1])    # 0.76    # 0.55

    # testing
    rf_serious = pd.read_csv('rf_serious.txt',sep='\t')
    ensemble = pd.read_csv('ensemble_result.txt',sep='\t')

    all4_data = np.empty([ensemble.shape[0],5])
    for i in range(all4_data.shape[0]):
        data[i,0] = rf_serious['rf_aaindex'].iloc[i]
        data[i,1] = rf_serious['rf_blosum'].iloc[i]
        data[i,2] = ensemble['aaindex'].iloc[i]
        data[i,3] = ensemble['blosum'].iloc[i]
        data[i,4] = ensemble['label'].iloc[i]

    result_test = clf.predict_proba(data[:,0:3])
    draw_ROC(data[:,4],result_test[:,1])
    draw_PR(data[:,4],result_test[:,1])

    # additional dataset

    def final_model(ori,hla,dic_inventory,after_pca=None):
        # data
        dataset1 = construct_aaindex(ori,hla,dic_inventory,after_pca)
        dataset2 = construct(ori,hla,dic_inventory)
        input1_b = pull_peptide(dataset2)
        input2_b = pull_hla(dataset2)
        input1_a = pull_peptide_aaindex(dataset1)
        input2_a = pull_hla_aaindex(dataset1)

        X1 = np.empty((len(dataset1), 12 * 56))  # 28581
        Y1 = ori['immunogenecity'].values
        for j, (x, y, _) in enumerate(dataset1):
            x = x.reshape(-1)  # 10*12*1 ---> 120
            y = y.reshape(-1)  # 46*12*1 ---> 552
            X1[j, :] = np.concatenate([x, y])  # 672

        X2 = np.empty((len(dataset2), 21 * 56))  # 28581
        Y2 = ori['immunogenecity'].values
        for j, (x, y, _) in enumerate(dataset2):
            x = x.reshape(-1)  # 10*12*1 ---> 120
            y = y.reshape(-1)  # 46*12*1 ---> 552
            X2[j, :] = np.concatenate([x, y])  # 672

        # load the model
        from joblib import load
        clf_rf_aaindex = load('RF_aaindex.joblib')
        clf_rf_blosum = load('RF_blosum.joblib')

        ResLikeCNN_index = model_aaindex()
        ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')
        ResLikeCNN = model()
        ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')

        # predict
        col1 = clf_rf_aaindex.predict_proba(X1)[:,1]
        col2 = clf_rf_blosum.predict_proba(X2)[:,1]
        col3 = ResLikeCNN_index.predict(x=[input1_a,input2_a])[:,0]
        col4 = ResLikeCNN.predict(x=[input1_b,input2_b])[:,0]
        print(col1.shape,col2.shape,col3.shape,col4.shape)
        data = np.stack([col1,col2,col3,col4],axis=1)

        # final
        clf_log = load('logisticRegression.joblib')
        final = clf_log.predict_proba(data)[:,1]
        return final

    addition = pd.read_csv('data/mannual_cancer_testing_fiter910.txt',sep='\t')
    a = final_model(addition,hla,dic_inventory,after_pca)
    draw_ROC(addition['immunogenecity'],a)
    draw_PR(addition['immunogenecity'],a)
    hard = [1 if i > 0.5 else 0 for i in a]
    confusion_matrix(addition['immunogenecity'],hard)
    accuracy_score(addition['immunogenecity'],hard)  # 0.65
    f1_score(addition['immunogenecity'],hard)   # 0.51
























































    class CustomCallback(keras.callbacks.Callback):
        def on_epoch_end(self, epoch, logs=None):
            if logs.get('accuracy') > 0.97:
                self.model.stop_training = True
    predictions = []
    for i in range(len(folds_df)):
        # ith fold will be prediction, other 9 will be training
        training = [data for m,data in enumerate(folds_df) if m != i]
        ori = pd.concat(training)
        dataset = construct_aaindex(ori,hla,dic_inventory,after_pca)
        input1 = pull_peptide_aaindex(dataset)
        input2 = pull_hla_aaindex(dataset)
        label = pull_label_aaindex(dataset)
        count0 = sum([True if item==0 else False for item in label[:,0]])
        count1 = sum([True if item==1 else False for item in label[:,0]])
        ratio = count1/(count0+count1)
        print(ratio,1-ratio)
        ResLikeCNN_index = model_aaindex()
        ResLikeCNN_index.compile(
            loss='binary_crossentropy',
            optimizer=keras.optimizers.Adam(lr=0.0001),
            metrics=['accuracy']
        )
        callback = keras.callbacks.EarlyStopping(monitor='loss', patience=3)
        history = ResLikeCNN_index.fit(
            x=[input1, input2],  # feed a list into
            y=label,
            batch_size=1024,
            epochs=200,
            class_weight={0: ratio, 1: 1-ratio},  # I have 20% positive and 80% negative in my training data
            callbacks=[CustomCallback(),callback]
        )

        ori_test = folds_df[i]
        testing_dataset = construct_aaindex(ori_test, hla, dic_inventory, after_pca)
        input1_test = pull_peptide_aaindex(testing_dataset)
        input2_test = pull_hla_aaindex(testing_dataset)
        label_test = pull_label_aaindex(testing_dataset)
        result = ResLikeCNN_index.predict(x=[input1_test, input2_test])
        ori_test['resnet_aaindex'] = result[:,0]
        predictions.append(ori_test)
    total_predictions = pd.concat(predictions)
    total_predictions.to_csv('total_prediction_resnet_aaindex.txt',sep='\t',index=None)




