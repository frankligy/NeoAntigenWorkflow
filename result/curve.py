#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:27:25 2020

@author: ligk2e
"""



import shelve
import torch
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix
import matplotlib.pyplot as plt
import pandas as pd



s = shelve.open('/Users/ligk2e/Desktop/immunogenecity/testing')
y_pred = s['y_pred']
predictions = s['predictions']
y = s['y']

diff = []
for i in range(len(y_pred)):
    diff.append(ast.literal_eval(np.format_float_scientific(y_pred[i,1])) - ast.literal_eval(np.format_float_scientific(y_pred[i,0])))



diff = np.subtract(y_pred[:,1],y_pred[:,0])


bundle = [i for i in zip(y_pred[:,1],predictions,y)]
sorted_bundle = np.array(sorted(bundle,key=lambda x: x[0]))

metrics = rank(sorted_bundle[:,[1,2]],step=10)


# rank

def rank(xy,step=10):
    metrics = []
    for i in range(len(xy)//step):
        query = xy[0:step*(i+1)]       
        confusion_mat = confusion_matrix(query[:,1],query[:,0],labels=[0,1])
        fpr = confusion_mat[0,1]/np.sum(confusion_mat[:,1])
        tpr = confusion_mat[0,0]/np.sum(confusion_mat[:,0])
        
        recall = tpr
        precision = confusion_mat[0,1]/np.sum(confusion_mat[0:])
        metrics.append((fpr,tpr,recall,precision))
    return metrics
        







# draw ROC curve
fpr,tpr,_ = roc_curve(y,y_pred[:,1],pos_label=1)
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
precision,recall,_ = precision_recall_curve(y,y_pred[:,1],pos_label=1)
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

### benchmark with









    




