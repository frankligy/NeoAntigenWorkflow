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
truth = pd.Series(a)
guide = dfori[truth]
guide.to_csv('feature_EventAnnotation.txt',sep='\t',index=None)


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


def inspectGTEx(dicTissueExp,event):
    tissueExp = dicTissueExp[event]
    for tissue, expression in tissueExp.items():
        
        
                            
    
    
with bz2.BZ2File('dicSRA.pbz2','wb') as f: 
    cPickle.dump(dicSRA, f)  

with bz2.BZ2File('dicSRA.pbz2','rb') as f1:
    data = cPickle.load(f1)    

import multiprocessing    
manager = multiprocessing.Manager()
Global = manager.Namespace()
Global.x = 10
Global.y = 'hello'
Global._z = 12.3    # this is an attribute of the proxy
print(Global)
Namespace(x=10, y='hello')    
    
    
    
import pandas as pd

sraData = pd.read_csv('GTEx_EventAnnotation.txt',sep='\t')

import pickle
with open('dicSRA.p','rb') as file1:
    dicSRA=pickle.load(file1)

dicTissueExp = {}
semaphore = 0

def search(ns):
    semaphore = ns.semaphore
    dicTissueExp = ns.dicTissueExp
    sraData = ns.sraData
    dicSRA = ns.dicSRA
    while semaphore < sraData.shape[0]:
        print('this is the {0}th run'.format(semaphore))
        event = sraData['UID'].tolist()[semaphore]
        for tissue,accID in dicSRA.items():
            try: dicTissueExp[event][tissue] = sraData.iloc[semaphore][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
            except KeyError:          
                dicTissueExp[event] = {}
                dicTissueExp[event][tissue] = sraData.iloc[semaphore][[accID+'_1.bed' for accID in dicSRA[tissue]]].values 
        semaphore += 1

def process_task(ns,lock): 
    for _ in range(sraData.shape[0]//4+1): 
        lock.acquire() 
        search(ns)
        lock.release() 

import multiprocessing
with multiprocessing.Manager() as manager:
    ns = manager.Namespace()
    ns.sraData = sraData
    ns.dicSRA = dicSRA
    ns.dicTissueExp = dicTissueExp
    ns.semaphore = semaphore



lock = multiprocessing.Lock()
p1 = multiprocessing.Process(target=process_task,args=(ns,lock,))
p1.start()
p1.join()
p2 = multiprocessing.Process(target=process_task,args=(ns,lock,))
p2.start()
p2.join()
p3 = multiprocessing.Process(target=process_task,args=(ns,lock,))
p3.start()
p3.join()
p4 = multiprocessing.Process(target=process_task,args=(ns,lock,))
p4.start()
p4.join()


p1.join()
p2.join()
p3.join()
p4.join()



import bz2
import _pickle as cpickle
with bz2.BZ2File('dicTissueExp1.pbz2','wb') as f1:
    cpickle.dump(dicTissueExp,f1)       
    
    
with bz2.BZ2File('dicTissueExp.pbz2','rb') as f2:
    dicTissueExp = cpickle.load(f2)    

import pickle    
from time import process_time
start = process_time()   
with open('dicTissueExp.p','rb') as f3:
    dicTissueExp = pickle.load(f3)
end = process_time()
print('consume {}s'.format(end-start))
    
a = [0,0.1,0.3,0.5,0.0,0.0,0.12]
fig = plt.figure()
plt.plot(np.arange(len(a)),a)
plt.ylim(-0.2,1.0)
plt.plot(np.arange(len(a)),np.zeros(len(a)),linestyle='dashed')

plt.bar(np.arange(len(a)),a,width=0.5,label='tissue')
plt.legend()


dicTissueExp['TSPAN6:ENSG00000000003:I3.1-E4.1|ENSG00000000003:E3.4-E4.1']['Cells - Cultured fibroblasts']
inspectGTEx('KYAT3:ENSG00000137944:E4.1-I4.1_88965413|ENSG00000137944:E4.1-E5.1')
inspectGTEx('NSUN5P2:ENSG00000106133:I2.1-E3.1|ENSG00000106133:E2.1-E3.1')   # inconsistent between GTEx and TCGA

tissue = list(dicSRA.keys())

def inspectGTEx(event,tissue='all',plot=True):
    flag = 0
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    import numpy as np
    global dicTissueExp
    if tissue=='all':
        tissueExp = dicTissueExp[event]
        for tis,exp in tissueExp.items():
            exp = exp.astype('float64')
            exp=exp[np.logical_not(np.isnan(exp))]
            if exp.size == 0: print('{0} data incomplete'.format(tis))
            elif np.any(exp):   # have non-zero element
                if plot==True:
                    fig = plt.figure()
                    plt.bar(np.arange(len(exp)),exp,width=0.2,label=tis)
                    plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                    plt.legend()
                    plt.savefig('./figures/{1}.pdf'.format(event,tis),bbox_inches='tight')
                    plt.close(fig)
                else: continue
            else: 
                flag += 1
                print('No expression in {}'.format(tis))
            
    else:
        expression = dicTissueExp[event][tissue]
        exp = expression.astype('float64')
        exp=exp[np.logical_not(np.isnan(exp))]
        if exp.size == 0: print('{0} data incomplete'.format(tissue))
        elif np.any(exp):   # have non-zero element
            plt.bar(np.arange(len(exp)),exp,width=0.2,label='tissue')
            plt.legend()
            plt.savefig('./{}.pdf'.format(tissue),bbox_inches='tight')
            plt.show()
            print(expression)
    return flag  
        
"""        
['Artery - Tibial',
 'Adipose - Subcutaneous',
 'Breast - Mammary Tissue',
 'Heart - Atrial Appendage',
 'Uterus',
 'Thyroid',
 'Adrenal Gland',
 'Lung',
 'Spleen',
 'Brain - Cortex',
 'Adipose - Visceral (Omentum)',
 'Colon - Sigmoid',
 'Skin - Not Sun Exposed (Suprapubic)',
 'Small Intestine - Terminal Ileum',
 'Whole Blood',
 'Esophagus - Gastroesophageal Junction',
 'Muscle - Skeletal',
 'Cells - Cultured fibroblasts',
 'Liver',
 'Minor Salivary Gland',
 'Skin - Sun Exposed (Lower leg)',
 'Stomach',
 'Ovary',
 'Heart - Left Ventricle',
 'Esophagus - Muscularis',
 'Nerve - Tibial',
 'Colon - Transverse',
 'Brain - Hypothalamus',
 'Brain - Nucleus accumbens (basal ganglia)',
 'Artery - Aorta',
 'Pituitary',
 'Testis',
 'Artery - Coronary',
 'Pancreas',
 'Brain - Cerebellum',
 'Kidney - Cortex',
 'Vagina',
 'Esophagus - Mucosa',
 'Prostate',
 'Brain - Caudate (basal ganglia)',
 'Brain - Putamen (basal ganglia)',
 'Cells - EBV-transformed lymphocytes',
 'Brain - Cerebellar Hemisphere',
 'Brain - Frontal Cortex (BA9)',
 'Brain - Spinal cord (cervical c-1)',
 'Brain - Amygdala',
 'Brain - Substantia nigra',
 'Brain - Hippocampus',
 'Brain - Anterior cingulate cortex (BA24)']        
"""        
        
        
        
        
        
        
        
        
        
        
            
            