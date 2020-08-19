#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:27:25 2020

@author: ligk2e
"""



import shelve
import torch
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score
import matplotlib.pyplot as plt
import pandas as pd



s = shelve.open('/Users/ligk2e/Desktop/immunogenecity/testing13')
y_pred = s['y_pred']
predictions = s['predictions']
y_mine = s['y']

confusion_matrix(y_mine,predictions)
f1_score(y_mine,predictions)



diff = y_pred[:,1] - y_pred[:,0]





# draw ROC curve
fpr_mine,tpr_mine,_ = roc_curve(y_mine,diff,pos_label=1)
area_mine = auc(fpr_mine,tpr_mine)
fig = plt.figure()
lw = 2
plt.plot(fpr_mine, tpr_mine, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % area_mine)
plt.plot(fpr, tpr, color='black',
         lw=lw, label='ROC curve (area = %0.2f)' % area)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()
fig.savefig('/Users/ligk2e/Desktop/immunogenecity/comp.svg')

# draw PR curve
precision,recall,_ = precision_recall_curve(y_mine,diff,pos_label=1)
area_PR = auc(recall,precision)
baseline = np.sum(np.array(y_mine) == 1) / len(y_mine)

plt.figure()
lw = 2
plt.plot(recall,precision, color='darkorange',
         lw=lw, label='PR curve (area = %0.2f)' % area_PR)
plt.plot([0, 1], [0.24, 0.24], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('PR curve example')
plt.legend(loc="lower right")
plt.show()


#### benchmark with IEDB
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

result = pd.read_csv('/Users/ligk2e/Desktop/immunogenecity/result.csv')
testing_data = pd.read_excel('/Users/ligk2e/Desktop/immunogenecity/testing_more.xlsx')
testing_data = clean_data_frame(testing_data)
testing_data.to_csv('/Users/ligk2e/Desktop/immunogenecity/testing_more_clean.txt',sep='\t',index=None)

result.sort_values('peptide',inplace=True)
testing_data.sort_values('peptide',inplace=True)

combine = pd.DataFrame({'peptide':result['peptide'],'predictions':result['score'],'y':testing_data['immunogenecity']})

# now we have the data

# draw ROC curve
fpr,tpr,_ = roc_curve(combine['y'].values,combine['predictions'].values,pos_label=1)
area = auc(fpr,tpr)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % area)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()

# draw PR curve
precision,recall,_ = precision_recall_curve(combine['y'].values,combine['predictions'].values,pos_label=1)
area_PR = auc(recall,precision)
baseline = np.sum(np.array(y) == 1) / len(y)

plt.figure()
lw = 2
plt.plot(recall,precision, color='darkorange',
         lw=lw, label='PR curve (area = %0.2f)' % area_PR)
plt.plot([0, 1], [0.1976, 0.1976], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('PR curve example')
plt.legend(loc="lower right")
plt.show()

### benchmark with netCTLpan1.1
df = pd.read_excel('/Users/ligk2e/Desktop/immunogenecity/testing_more.xlsx')
y = df['immunogenecity'].values
y_pred = df['netCTLpan'].values.astype('float64')

# draw ROC curve
fpr,tpr,_ = roc_curve(y,y_pred,pos_label=1)
area = auc(fpr,tpr)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % area)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()

# draw PR curve
precision,recall,_ = precision_recall_curve(y,y_pred,pos_label=1)
area_PR = auc(recall,precision)
baseline = np.sum(np.array(y) == 1) / len(y)

plt.figure()
lw = 2
plt.plot(recall,precision, color='darkorange',
         lw=lw, label='PR curve (area = %0.2f)' % area_PR)
plt.plot([0, 1], [0.1976, 0.1976], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('PR curve example')
plt.legend(loc="lower right")
plt.show()


# benchmark with deephlapan
df = pd.read_csv('/Users/ligk2e/Desktop/immunogenecity/value.txt',sep='\t',names=['y','y_pred'])
y = df['y'].values
y_pred = df['y_pred'].values
# draw ROC curve
fpr,tpr,_ = roc_curve(y,y_pred,pos_label=1)
area = auc(fpr,tpr)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % area)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()

# draw PR curve
precision,recall,_ = precision_recall_curve(y,y_pred,pos_label=1)
area_PR = auc(recall,precision)
baseline = np.sum(np.array(y) == 1) / len(y)

plt.figure()
lw = 2
plt.plot(recall,precision, color='darkorange',
         lw=lw, label='PR curve (area = %0.2f)' % area_PR)
plt.plot([0, 1], [0.1976, 0.1976], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('PR curve example')
plt.legend(loc="lower right")
plt.show()


#new dataset with deephlapan
df1 = pd.read_csv('/Users/ligk2e/Desktop/immunogenecity/ineo_testing_new.txt',sep='\t')
df2 = pd.read_csv('/Users/ligk2e/Desktop/immunogenecity/ineo_testing_new_final_predicted_result.csv')
y=df1['immunogenecity'].values
y_pred=df2['immunogenic score'].values

hard_y = [1 if i > 0.5 else 0 for i in y_pred]

counter = 0
for i in range(len(y)):
    if y[i] == hard_y[i]: counter += 1

confusion_matrix(y,hard_y)
f1_score(y,hard_y)
# accuracy 58.92%

# draw ROC curve
fpr,tpr,_ = roc_curve(y,y_pred,pos_label=1)
area = auc(fpr,tpr)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='black',
         lw=lw, label='ROC curve (area = %0.2f)' % area)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()

# draw PR curve
precision,recall,_ = precision_recall_curve(y,y_pred,pos_label=1)
area_PR = auc(recall,precision)
baseline = np.sum(np.array(y) == 1) / len(y)

plt.figure()
lw = 2
plt.plot(recall,precision, color='darkorange',
         lw=lw, label='PR curve (area = %0.2f)' % area_PR)
plt.plot([0, 1], [0.24, 0.24], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('PR curve example')
plt.legend(loc="lower right")
plt.show()










    




