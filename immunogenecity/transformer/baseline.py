#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:45:20 2020

@author: ligk2e
"""

from utils import *
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC,LinearSVC

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
ori_test =pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/ineo_testing_filter910.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
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


# SVM, linear
clf = LinearSVC(random_state=0).fit(X,Y)
predictions = clf.predict(X_test)
y_pred = clf.predict_proba(X_test)
draw_ROC(Y_test,y_pred[:,1])   
draw_PR(Y_test,y_pred[:,1])  
confusion_matrix(Y_test,predictions)
f1_score(Y_test,predictions)
accuracy_score(Y_test,predictions)




 
    
    
    

