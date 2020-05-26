#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 16:31:10 2020

@author: ligk2e
"""

import pandas as pd
import os
import sys
import numpy as np

df = pd.read_csv('/Users/ligk2e/Downloads/BC_vs_Control.txt',sep='\t')
groups = pd.read_csv('/Users/ligk2e/Downloads/groups.txt',sep='\t',header=None)

dic = {}
for i in range(groups.shape[0]):
    idn= groups.iloc[i][2]
    TCGA = groups.iloc[i][0]
    try: dic[idn].append(idn + ':' + TCGA)
    except KeyError: 
        dic[idn] = []
        dic[idn].append(idn + ':' + TCGA)

col1,col2,col3,col4 = [],[],[],[]        
for j in range(df.shape[0]):
    BC = np.nanmean(df.iloc[j][dic['BC']].tolist())
    CT = np.nanmean(df.iloc[j][dic['Control']].tolist())
    diff = BC-CT
    cond = True if diff > 0.01 else False
    col1.append(BC)
    col2.append(CT)
    col3.append(diff)
    col4.append(cond)
df['BC_mean'] = col1
df['CT_mean'] = col2
df['dPSI'] = col3
df['cond'] = col4

df_new = df[df['cond']]

df_out = df_new[['UID','dPSI','BC_mean','CT_mean']]
df_out.to_csv('/Users/ligk2e/Downloads/pseudoPSImatrix.txt',sep='\t',index=None)


with open('/Users/ligk2e/Downloads/test.out.txt','r') as f3:
    next(f3)
    punchline = f3.readline().rstrip('\n').split(' ')
    TMn = int(punchline[-1])
    


def diffNovelFromNotInvolved(df):
    col = []

    for i in range(df.shape[0]):
        cond = True
        exam_match_whole_tran = df['exam_match_whole_tran'].iloc[i]   # will be a list already, 1*1
        for item in ast.literal_eval(exam_match_whole_tran):
            if item: 
                cond = False
                break
        col.append(cond)
    print(len(col),df.shape[0],col)
    print(len(pd.Series(col)))
    try:df_novel = df[pd.Series(col)]
    except:
        print(len(col),df.shape[0],col)
        raise Exception('jkjk')
    return df_novel  

df = pd.read_csv('/Users/ligk2e/Desktop/project_LUAD/fjdfjd.txt',sep='\t')
df_1 = diffNovelFromNotInvolved(df)


a = df['exam_match_whole_tran'].iloc[0]



# df_out, groups,df
def getPercent(df,dfgroup,dfori,write=False):
    dic = {}
    for i in range(dfgroup.shape[0]):
        id_, label = dfgroup[0].tolist()[i],dfgroup[2].tolist()[i]
        if label == 'BC':
            id_ = 'BC'+':'+id_
            try: dic['BC'].append(id_)
            except KeyError:
                dic['BC'] = []
                dic['BC'].append(id_)
        elif label == 'Control':
            id_ = 'Control'+':'+id_
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
        truth_t = [True if item == 0 else False for item in allvalue_t]
        allzero_t = sum(truth_t)   # how many zeros in this cluster
        percent_t = allzero_t/num_t  # percentage of zeros
        percentArray_t.append(percent_t)
            
        #allvalue = dfori[dic['R1-V7']].iloc[j].tolist()
        allvalue_h = list(map(lambda x:round(x,2),dfori.iloc[j][dic['Control']].tolist()))
        truth_h = [True if item == 0 else False for item in allvalue_h]
        allzero_h = sum(truth_h)   # how many zeros in this cluster
        percent_h = allzero_h/num_h   # percentage of zeros
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
    df['tumor_zero_percent'] = col0
    df['tumor_distribution'] = col1
    df['healthy_zero_percent'] = col2
    df['healthy_distribution'] = col3
    return df

df_9 = getPercent(df_out,groups,df)
df_9.to_csv('/Users/ligk2e/Downloads/pseudoPSImatrix_percentage.txt',sep='\t',index=None)
 

