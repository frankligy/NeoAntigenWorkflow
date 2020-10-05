#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:45:20 2020

@author: ligk2e
"""


import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from Bio.SubsMat import MatrixInfo
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC,LinearSVC
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler,MinMaxScaler,RobustScaler
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import pickle
import shelve
import matplotlib.pyplot as plt

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
        #self.new = self.padding()
        self.new = self.padding_oneside()
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

    def padding_oneside(self):
        len_values = [tup[0].shape[0] for tup in self.middle]
        #max_length = max(len_values)  
        max_length = 56      
        # padding
        bucket = []
        for item in self.middle:

            length = item[0].shape[0]
            gap = max_length - length
                    

            padding_right = torch.empty([gap,21]).fill_(-1.0)
            final = torch.cat([item[0],padding_right],dim=0)
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
    def blosum50_new(peptide):
        amino = 'ARNDCQEGHILKMFPSTWYV-'
        dic = MatrixInfo.blosum50
        matrix = np.zeros([21,21])
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                try:
                    matrix[i,j] = dic[(amino[i],amino[j])] 
                except KeyError:
                    try:
                        matrix[i,j] = dic[(amino[j],amino[i])]
                    except:
                        matrix[i,j] = -1
                    
        encoded = torch.empty([len(peptide),21])       # (seq_len,21)       
        for i in range(len(peptide)):

            encoded[i,:] = torch.from_numpy(matrix[:,amino.index(peptide[i])])
                
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
            X = dataset.blosum50_new(cat).float()   # 2-d tensor
            #X = dataset.onehot_classic(cat).float()
            y = torch.tensor(immuno).long()  # 0-d tensor
            lis.append((X,y))
        return lis

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

def paratope_dic(hla):
    df = hla
    dic = {}
    for i in range(df.shape[0]):
        hla = df['hla'].iloc[i]
        paratope = df['paratope'].iloc[i]
        dic[hla] = paratope
    return dic

if __name__ == '__main__':
    # load the data

    ori = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    training_dataset = dataset(ori,hla,dic_inventory)


    # convert to scikit-learn format, 21*56 features
    X = np.empty((len(training_dataset),21*56))   # 28579 * 1176
    Y = np.empty(len(training_dataset))           # 28579
    for i,(x,y) in enumerate(training_dataset):
        x = x.reshape(-1)   # flatten row by row
        X[i,:] = x
        Y[i] = y 

    # testing 
    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = dataset(ori_test,hla,dic_inventory)
    X_test = np.empty((len(testing_dataset),21*56))
    Y_test = np.empty(len(testing_dataset))
    for i,(x,y) in enumerate(testing_dataset):
        x = x.reshape(-1)
        X_test[i,:] = x
        Y_test[i] = y


    # logistic regression
    clf = LogisticRegression(max_iter=1000).fit(X,Y)
    predictions = clf.predict(X_test)
    y_pred = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred[:,1])   
    draw_PR(Y_test,y_pred[:,1])  
    confusion_matrix(Y_test,predictions)
    f1_score(Y_test,predictions)
    accuracy_score(Y_test,predictions)

    # serialize the y_pred and Y_test
    s = shelve.open('logistic')
    s['y_pred'] = y_pred
    s['Y_test'] = Y_test
    s.close()

    # SVM, linear
    clf = LinearSVC(random_state=0).fit(X,Y)
    predictions = clf.predict(X_test)
    y_pred = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred[:,1])   
    draw_PR(Y_test,y_pred[:,1])  
    confusion_matrix(Y_test,predictions)
    f1_score(Y_test,predictions)
    accuracy_score(Y_test,predictions)

    # Random Forest
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier()
    clf.fit(X,Y)
    y_pred_blosum = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred_blosum[:,1])
    draw_PR(Y_test,y_pred_blosum[:,1])
    predictions = clf.predict(X_test)
    confusion_matrix(Y_test,predictions)
    f1_score(Y_test,predictions)
    accuracy_score(Y_test,predictions)

    ## It performs so well, let's look at aaindex performance
    import re
    import collections
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    dataset = construct_aaindex(ori,hla,dic_inventory,after_pca)
    X = np.empty((len(dataset), 12 * 56))  # 28581
    Y = ori['immunogenecity'].values
    for i, (x, y,_) in enumerate(dataset):
        x = x.reshape(-1)  # 10*12*1 ---> 120
        y = y.reshape(-1)  # 46*12*1 ---> 552
        X[i,:] = np.concatenate([x,y])   # 672

    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = construct_aaindex(ori_test,hla,dic_inventory,after_pca)
    X_test = np.empty((len(testing_dataset),12*56))
    Y_test = ori_test['immunogenecity'].values
    for i,(x,y,_) in enumerate(testing_dataset):
        x = x.reshape(-1)
        y = y.reshape(-1)
        X_test[i,:] = np.concatenate([x,y])

    clf = RandomForestClassifier()
    clf.fit(X,Y)
    y_pred_aaindex = clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred[:,1])
    draw_PR(Y_test,y_pred[:,1])
    predictions = clf.predict(X_test)
    confusion_matrix(Y_test,predictions)
    f1_score(Y_test,predictions)
    accuracy_score(Y_test,predictions)

    random_forest_no_depth_restriction = pd.DataFrame({'label':Y_test,'rf_blosum':y_pred_blosum[:,1],'rf_aaindex':y_pred_aaindex[:,1]})
    random_forest_no_depth_restriction.to_csv('rf_no_depth.txt',sep='\t',index=None)

    # record when you don't use max_depth
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier(max_depth=50)
    clf.fit(X,Y)
    y_pred_50= clf.predict_proba(X_test)
    draw_ROC(Y_test,y_pred_50[:,1])
    draw_PR(Y_test,y_pred_50[:,1])
    predictions = clf.predict(X_test)
    confusion_matrix(Y_test,predictions)
    f1_score(Y_test,predictions)
    accuracy_score(Y_test,predictions)

    # add them to random_forest folder
    random_forest_no_depth_restriction['aaindex_max2'] = y_pred_2[:,1]
    random_forest_no_depth_restriction['aaindex_max3'] = y_pred_3[:,1]
    random_forest_no_depth_restriction['aaindex_max4'] = y_pred_4[:,1]
    random_forest_no_depth_restriction['aaindex_max10'] = y_pred_10[:,1]
    random_forest_no_depth_restriction['aaindex_max50'] = y_pred_50[:,1]

    random_forest_no_depth_restriction.to_csv('rf_no_depth.txt',sep='\t',index=None)



    # ANN
    class ANN(nn.Module):
        def __init__(self):
            super(ANN,self).__init__()
            self.fc1 = nn.Linear(1176,128)
            self.relu = nn.ReLU()
            self.fc2 = nn.Linear(128,2)
            self.sigmoid = nn.Sigmoid()

        def forward(self,x):
            out = self.sigmoid(self.fc2(self.relu(self.fc1(x))))
            return out


    class c1(Dataset):
        def __init__(self,X,Y):
            self.x = X
            self.y = Y

        def __len__(self):
            return self.x.shape[0]

        def __getitem__(self,idx):
            x = torch.from_numpy(self.x[idx])
            y = torch.tensor(self.y[idx])
            return (x,y)


    ori = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    training_dataset = dataset(ori,hla,dic_inventory)


    # convert to scikit-learn format, 21*56 features
    X = np.empty((len(training_dataset),21*56))   # 28579 * 1176
    Y = np.empty(len(training_dataset))           # 28579
    for i,(x,y) in enumerate(training_dataset):
        x = x.reshape(-1)   # flatten row by row
        X[i,:] = x
        Y[i] = y

    # testing
    ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    testing_dataset = dataset(ori_test,hla,dic_inventory)
    X_test = np.empty((len(testing_dataset),21*56))
    Y_test = np.empty(len(testing_dataset))
    for i,(x,y) in enumerate(testing_dataset):
        x = x.reshape(-1)
        X_test[i,:] = x
        Y_test[i] = y


    c1_instance = c1(X,Y)
    c1_loader = DataLoader(c1_instance,batch_size=128,shuffle=True)
    model = ANN()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    loss_f=nn.CrossEntropyLoss()
    for epoch in range(30):
        for i in c1_loader:
            x = i[0].float()
            y = i[1].long().reshape(-1)
            optimizer.zero_grad()

            y_pred = model(x)
            print(y_pred.size(),y.size())
            loss = loss_f(y_pred,y)
            loss.backward()
            optimizer.step()


    c2_instance = c1(X_test,Y_test)
    c2_loader = DataLoader(c2_instance,batch_size=654)
    model.eval()
    with torch.no_grad():
        for i in c2_loader:
            x = i[0].float()
            y = i[1].long().reshape(-1)
            y_pred = model(x)

    # serialize the y_pred and Y_test
    s = shelve.open('ANN')
    s['y_pred'] = y_pred.numpy()
    s['Y_test'] = Y_test
    s.close()

    draw_ROC(Y_test,y_pred[:,1])
    draw_PR(Y_test,y_pred[:,1])











 
    
    
    

