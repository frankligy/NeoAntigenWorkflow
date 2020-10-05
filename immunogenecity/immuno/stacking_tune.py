'''
Tune the parameters of random forest and ANN
'''

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pandas as pd
import os
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler,MinMaxScaler,RobustScaler
from sklearn.model_selection import cross_val_score
import collections
from Bio.SubsMat import MatrixInfo
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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

    # 10-fold cross-validation
    idx = np.arange(immuno.shape[0])  # all index
    np.random.shuffle(idx)   # shuffle
    folds_idx = np.array_split(idx,10)
    folds_df = []
    for fold in folds_idx:
        sub = immuno.iloc[fold].set_index(pd.Index(np.arange(len(fold))))
        folds_df.append(sub)

    # First tune random forest
    ## aaindex
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    dataset = construct_aaindex(immuno, hla, dic_inventory, after_pca)
    X = np.empty((len(dataset), 12 * 56))  # 28581
    Y = immuno['immunogenecity'].values
    for i, (x, y, _) in enumerate(dataset):
        x = x.reshape(-1)  # 10*12*1 ---> 120
        y = y.reshape(-1)  # 46*12*1 ---> 552
        X[i, :] = np.concatenate([x, y])  # 672
    tree = [10,50,100,200]
    scores = []
    for i in tree:
        clf = RandomForestClassifier(n_estimators=i)
        score = cross_val_score(clf,X,Y,cv=10,scoring='accuracy')
        scores.append(score.mean())
    '''
    accuracy:
    [0.8658192090395481, 0.8731010671688638, 0.8746076585059634, 0.8757376020087886]
    ideas:
    given the long time it will cost, probably default value 100 is a decent number of trees
    '''
    clf = RandomForestClassifier()
    clf.fit(X,Y)
    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = construct_aaindex(ori_test,hla,dic_inventory,after_pca)
    X_test = np.empty((len(testing_dataset),12*56))
    Y_test = ori_test['immunogenecity'].values

    for i,(x,y,_) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i,:] = np.concatenate([x,y])
    y_pred_ann_aaindex = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred_ann_aaindex[:,1])
    draw_PR(Y_test,y_pred_ann_aaindex[:,1])

    ## BLOSUM50
    dataset = construct(immuno,hla,dic_inventory)
    X = np.empty((len(dataset), 21 * 56))
    Y = immuno['immunogenecity'].values
    for i, (x, y, _) in enumerate(dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X[i, :] = np.concatenate([x, y])  # 1176
    tree = [10,50,100,200]
    scores = []
    for i in tree:
        clf = RandomForestClassifier(n_estimators=i)
        score = cross_val_score(clf,X,Y,cv=10,scoring='accuracy')
        scores.append(score.mean())
    '''
    accuracy:
    [0.8657564344005022, 0.874074074074074, 0.8760514752040175, 0.876585059635907]
    ideas:
    the same, choose default value 100
    '''
    clf = RandomForestClassifier()
    clf.fit(X,Y)
    # let's have a test
    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = construct(ori_test,hla,dic_inventory)
    X_test = np.empty((len(testing_dataset),21*56))
    Y_test = ori_test['immunogenecity'].values

    for i,(x,y,_) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i,:] = np.concatenate([x,y])
    y_pred_ann_blosum = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred_ann_blosum[:,1])
    draw_PR(Y_test,y_pred_ann_blosum[:,1])


    # ANN (multi-layer perceptron in scikit-learn)
    ## aaindex
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    dataset = construct_aaindex(immuno, hla, dic_inventory, after_pca)
    X = np.empty((len(dataset), 12 * 56))  # 28581
    Y = immuno['immunogenecity'].values
    for i, (x, y, _) in enumerate(dataset):
        x = x.reshape(-1)  # 10*12*1 ---> 120
        y = y.reshape(-1)  # 46*12*1 ---> 552
        X[i, :] = np.concatenate([x, y])  # 672
    hidden = [16,32,64,128]
    scores = []
    for i in hidden:
        clf = MLPClassifier(hidden_layer_sizes=(i,),max_iter=500,early_stopping=True)
        score = cross_val_score(clf,X,Y,cv=10,scoring='accuracy')
        scores.append(score.mean())
    '''
    accuracy:
    [0.8519460138104206, 0.8534526051475204, 0.8536095417451349, 0.8576271186440678]
    let's try 128 as hidden_layer_size
    '''
    clf = MLPClassifier(hidden_layer_sizes=(128,),max_iter=500,early_stopping=True)
    clf.fit(X,Y)  # 25
    clf.loss_curve_
    '''
    [0.4668304867568887, 0.3460529378300128, 0.31972987917553347, 0.2952449079642726, 0.2831002537316494, 
    0.27636176488426945, 0.26015989043934223, 0.24933265643487557, 0.24093220946298638, 0.22859117598708353, 
    0.22079888485300375, 0.21276676696297417, 0.20455013463303684, 0.19579674945347855, 0.1915102361384166, 
    0.17913638494855647, 0.18090662568980403, 0.16516021327116226, 0.160718476654002, 0.1587605899495024, 
    0.15126680306897233, 0.14966821340923941, 0.1430078797582131, 0.13576514987668298, 0.1319773838652472]
    '''
    # let's have a test
    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = construct_aaindex(ori_test,hla,dic_inventory,after_pca)
    X_test = np.empty((len(testing_dataset),12*56))
    Y_test = ori_test['immunogenecity'].values

    for i,(x,y,_) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i,:] = np.concatenate([x,y])
    y_pred_ann_aaindex = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred_ann_aaindex[:,1])
    draw_PR(Y_test,y_pred_ann_aaindex[:,1])

    ## BLOSUM
    dataset = construct(immuno,hla,dic_inventory)
    X = np.empty((len(dataset), 21 * 56))
    Y = immuno['immunogenecity'].values
    for i, (x, y, _) in enumerate(dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X[i, :] = np.concatenate([x, y])  # 1176
    hidden = [16,32,64,128]
    scores = []
    for i in hidden:
        clf = MLPClassifier(hidden_layer_sizes=(i,),max_iter=500,early_stopping=True)
        score = cross_val_score(clf,X,Y,cv=10,scoring='accuracy')
        scores.append(score.mean())
    '''
    [0.8537037037037036, 0.8549591964846203, 0.8540489642184557, 0.8580037664783429]
    choose 128 as well
    '''
    clf = MLPClassifier(hidden_layer_sizes=(128,),max_iter=500,early_stopping=True)
    clf.fit(X,Y)  # 60
    clf.loss_curve_

    # let's have a test
    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = construct(ori_test,hla,dic_inventory)
    X_test = np.empty((len(testing_dataset),21*56))
    Y_test = ori_test['immunogenecity'].values

    for i,(x,y,_) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i,:] = np.concatenate([x,y])
    y_pred_ann_blosum = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred_ann_blosum[:,1])
    draw_PR(Y_test,y_pred_ann_blosum[:,1])






