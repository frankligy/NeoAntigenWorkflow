#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 20:39:42 2020

@author: ligk2e
"""
import pandas as pd
import os
os.chdir('/Users/ligk2e/Desktop/project_LUAD')

dfori = pd.read_csv('LUAD_Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-filtered-75p.txt',sep='\t')
dfgroup = pd.read_csv('groups.txt',sep='\t',header=None,names=['TCGA-ID','group','label'])


dic = {}
for i in range(dfgroup.shape[0]):
    id_, label = dfgroup['TCGA-ID'].tolist()[i],dfgroup['label'].tolist()[i]
    if label == 'R1-V7':
        try: dic['R1-V7'].append(id_)
        except KeyError:
            dic['R1-V7'] = []
            dic['R1-V7'].append(id_)
    elif label == 'Healthy':
        try: dic['Healthy'].append(id_)
        except KeyError:
            dic['Healthy'] = []
            dic['Healthy'].append(id_)

num_t = len(dic['R1-V7'])
num_h = len(dic['Healthy'])
 
percentArray_h = [] 
percentArray_t = [] 
dicPercent = {}        
for j in range(dfori.shape[0]):
    event = dfori['UID'].tolist()[j]    # next 'UID' was converted to 'UID.1'
    
    allvalue_t = dfori.iloc[j][dic['R1-V7']].tolist()
    truth_t = [True if item == 0 else False for item in allvalue_t]
    allzero_t = sum(truth_t)   # how many zeros in this cluster
    percent_t = allzero_t/num_t  # percentage of zeros
    percentArray_t.append(percent_t)
        
    #allvalue = dfori[dic['R1-V7']].iloc[j].tolist()
    allvalue_h = dfori.iloc[j][dic['Healthy']].tolist()
    truth_h = [True if item == 0 else False for item in allvalue_h]
    allzero_h = sum(truth_h)   # how many zeros in this cluster
    percent_h = allzero_h/num_h   # percentage of zeros
    percentArray_h.append(percent_h)
    dicPercent[event]=(percent_t,percent_h)

    
a = [True if i==1 else False for i in percentArray_h]



dfquery = pd.read_csv('PSI.R1-V7_vs_Healthy.txt',sep='\t')
col0,col1 = [],[]
for k in range(dfquery.shape[0]):
    splice = dfquery['UID'].tolist()[k]
    per_t,per_h = dicPercent[splice][0],dicPercent[splice][1]
    col0.append(per_t)
    col1.append(per_h)
dfquery['tumor_zero_percent'] = col0
dfquery['healthy_zero_percent'] = col1
dfquery.to_csv('see.txt',sep='\t',index=None)


    

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random

fig = plt.figure()
plt.plot(np.arange(100),random.choices(percentArray_t,k=100),color='g',alpha=0.3)
plt.plot(np.arange(100),random.choices(percentArray_h,k=100),color='r',alpha = 0.3)


plt.plot(np.arange(100),percentArray_t[0:100],color='g',alpha=0.3)
plt.plot(np.arange(100),percentArray_h[0:100],color='r',alpha = 0.3)


sraTable1 = pd.read_csv('SraRunTable-GTEX1.txt',sep=',')
sraTable2 = pd.read_csv('SraRunTable-GTEX2.txt',sep=',')
sraTable3 = pd.read_csv('SraRunTable-GTEX3.txt',sep=',')
sraData = pd.read_csv('GTEx_EventAnnotation.txt',sep='\t')

dicSRA = {}
for i in range(sraTable1.shape[0]):
    accID, tissue= sraTable1['Run'].tolist()[i],sraTable1['body_site'].tolist()[i]
    try: dicSRA[tissue].append(accID)
    except KeyError:
        dicSRA[tissue] = []
        dicSRA[tissue].append(accID)

for i in range(sraTable2.shape[0]):
    accID, tissue= sraTable2['Run'].tolist()[i],sraTable2['body_site'].tolist()[i]
    try: dicSRA[tissue].append(accID)
    except KeyError:
        dicSRA[tissue] = []
        dicSRA[tissue].append(accID)
    
for i in range(sraTable3.shape[0]):
    accID, tissue= sraTable3['Run'].tolist()[i],sraTable3['body_site'].tolist()[i]
    try: dicSRA[tissue].append(accID)
    except KeyError:
        dicSRA[tissue] = []
        dicSRA[tissue].append(accID)   
        
import shelve
with shelve.open('bytefile') as db:   # DbfilenameShelf object
    db['dicSRA'] = dicSRA

import pickle
with open('dicSRA.p','wb') as file2:
    pickle.dump(dicSRA,file2)


with shelve.open('bytefile') as db:
    dicSRA = db['dicSRA']
    
with open('dicSRA.p','rb') as file1:
    dicSRA = pickle.load(file1)


dicTissueExp = {}
for i in range(sraData.shape[0]):
    print('this is the {0}th run'.format(i))
    event = sraData['UID'].tolist()[i]
    for tissue,accID in dicSRA.items():
        try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
        except KeyError:            
            dicTissueExp[event] = {}
            dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
        
        
with shelve.open('bytefile') as db:
    dicSRA = db['dicSRA']
    dicTissueExp = db['dicTissueExp']

with open('dicTissueExp.p','rb') as file2:
    dicTissueExp = pickle.load(file2)    

def inspectGTEx(df)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




    
            
            