#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 18:43:41 2020

@author: ligk2e
"""


import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt





################## Functions ###############################
def check(key,counts_h):  # check the all counts for a given event, counts_h contains all counts for matched control
    for i in range(counts_h.shape[0]):
        event = counts_h.iloc[i][0]
        if key in event:
            cond = counts_h.iloc[i].isnull().values.any()
            detail = counts_h.iloc[i].tolist()[1:]
            break
    return cond,detail
    
# dfori is EventAnnotation file, df is PSI matrix file after metadata analysis    
def getPercent(df,dfgroup,dfori,write=False): # get percentage information for each event (Among healthy and tumor samples, how many of them are non-zero)
    dic = {}
    for i in range(dfgroup.shape[0]):
        id_, label = dfgroup[0].tolist()[i],dfgroup[2].tolist()[i]
        if label == 'BC':     # label needs to be changed
            #id_ = 'BC'+':'+id_
            try: dic['BC'].append(id_)
            except KeyError:
                dic['BC'] = []
                dic['BC'].append(id_)
        elif label == 'Control':
            #id_ = 'Control'+':'+id_
            try: dic['Control'].append(id_)
            except KeyError:
                dic['Control'] = []
                dic['Control'].append(id_)
    
    num_t = len(dic['BC'])
    num_h = len(dic['Control'])
 
    percentArray_h = [] 
    percentArray_t = [] 
    dicPercent = {}        
    for j in range(dfori.shape[0]):
        event = dfori['UID'].tolist()[j]    # next 'UID' was converted to 'UID.1'
        
        allvalue_t = list(map(lambda x:round(x,2),dfori.iloc[j][dic['BC']].tolist()))
        nonzero_t = len(np.nonzero(allvalue_t)[0])
        nan_t = sum(np.isnan(allvalue_t))
        zero_t = num_t - non_zero_t - nan_t
        percent_t = nonzero_t/num_t  # percentage of non-zeros, you cound change it to zeros or nans
        percentArray_t.append(percent_t)
            
        #allvalue = dfori[dic['R1-V7']].iloc[j].tolist()
        allvalue_h = list(map(lambda x:round(x,2),dfori.iloc[j][dic['Control']].tolist()))
        nonzero_h = len(np.nonzero(allvalue_h)[0])
        nan_h = sum(np.isnan(allvalue_h))
        zero_h = num_h - non_zero_h - nan_h
        percent_h = nonzero_h/num_h  # percentage of non-zeros, you cound change it to zeros or nans
        percentArray_h.append(percent_h)
        dicPercent[event]=(percent_t,allvalue_t,percent_h,allvalue_h)
        
    col0,col1,col2,col3 = [],[],[],[]
    for k in range(df.shape[0]):
        splice = df['UID'].tolist()[k]
        per_t,all_t,per_h,all_h = dicPercent[splice][0],dicPercent[splice][1],dicPercent[splice][2],dicPercent[splice][3]
        col0.append(per_t)
        col1.append(all_t)
        col2.append(per_h)
        col3.append(all_h)
    df['tumor_percent'] = col0
    df['tumor_distribution'] = col1
    df['healthy_percent'] = col2
    df['healthy_distribution'] = col3
    if write==True: df.to_csv('PSImatrix_percentage.txt',sep='\t',index=None)
    return df


df = pd.read_csv('/Users/ligk2e/Downloads/pseudoPSImatrix.txt',sep='\t')
dfgroup = pd.read_csv('/Users/ligk2e/Downloads/groups.txt',sep='\t',header=None)
df_ori = pd.read_csv('/Users/ligk2e/Downloads/BC_vs_Control.txt',sep='\t')
df_9 = getPercent(df,dfgroup,df_ori,True)  
    
    
    
    
    
    
    
    
    
    
    






################## Code Chunks #############################
'''
chunk1: give counts file, get seperate counts for healthy sample and tumor samples that are of interest.
code below is for breast cancer TCGA dataset
'''    
groups = pd.read_csv('/Users/ligk2e/Downloads/groups.txt',sep='\t',header=None)   
counts = pd.read_csv('/Users/ligk2e/Desktop/project_breast/R1-V6/counts.TCGA-BRCA.txt',sep='\t')
fields = counts.columns.tolist()  # all colnames

healthy,tumor = [],[]
for i in range(groups.shape[0]):
    label = groups.iloc[i][2]
    id_ = groups.iloc[i][0]
    if label=='Control': healthy.append(id_)   # need change label name
    elif label=='BC': tumor.append(id_)        # need change label name
    
result = []       # which fileds will be grouped to healthy
for j in fields:
    truth = 0    # assume don't pick up any one
    for k in healthy:
        if k in j:
            truth=1
            break
    result.append(truth)

result[0] = 1       # first filed is ID, so will always include it
index = np.nonzero(result)[0].tolist()    # a list contains all indices where the column will retain
                                # np.nonzero return a tuple, only pick the first one [0]


number_col = np.arange(1224)    # change colname to num-based, need change total number  
counts.columns=number_col

counts_h = counts[index]   

counts_h.to_csv('/Users/ligk2e/Desktop/project_breast/R1-V6/counts.TCGA-BRCA.healthy.txt',sep='\t',index=None)


'''
chunk2: piece together all MS-supported splicing event, write as a new dataframe

'''
import pandas as pd
import ast
pep = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/peptides.txt',sep='\t')['Sequence'].tolist()
idx8,idx9,idx10,idx11=[],[],[],[]
for item in pep:
    length = len(item)
    if length==8:
        mer8_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_8_mark.txt',sep='\t')
        mer8 = mer8_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer8):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-B57:01']
                if [item] in candidate: idx8.append(idx)

    if length==9:
        mer9_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_9_mark.txt',sep='\t')
        mer9 = mer9_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer9):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-B57:01']
                if [item] in candidate: idx9.append(idx)
       
    if length==10:
        mer10_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_10_mark.txt',sep='\t')
        mer10 = mer10_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer10):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-B57:01']
                if [item] in candidate: idx10.append(idx)
   
    if length==11:
        mer11_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_11_mark.txt',sep='\t')
        mer11 = mer11_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer11):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-B57:01']
                if [item] in candidate: idx11.append(idx)

        
new8_df = mer8_df.iloc[idx8]
new9_df = mer9_df.iloc[idx9] 
new10_df = mer10_df.iloc[idx10] 
new11_df = mer11_df.iloc[idx11] 
new = pd.concat([new8_df,new9_df,new10_df,new11_df])

new.to_csv('/Users/ligk2e/Desktop/new.txt',sep='\t',index=None)

'''
chunk3: get percentage information for each event (Among healthy and tumor samples, how many of them are non-zero)
'''
    
    
    
    