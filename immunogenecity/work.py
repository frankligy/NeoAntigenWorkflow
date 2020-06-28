#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 12:09:54 2020

@author: ligk2e
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/ligk2e/Desktop/immunogenecity')
from utility import *
import scipy.stats as sc
from mlxtend.evaluate import permutation_test
import collections

positive = pd.read_excel('/Users/ligk2e/Desktop/immunogenecity/HC neoantigens.xlsx',skiprows=[0],index=None,nrows=295)
negative = pd.read_excel('/Users/ligk2e/Desktop/immunogenecity/negative_collection/negative_filtered.xlsx',index=None)

pep_pos = positive['Mut peptide'].tolist()
pep_neg = negative['peptides'].tolist()

# process the pep_pos and pep_neg
pep_neg = [item.replace(' ','').rstrip(' ') for item in pep_neg]
pep_pos = [item.replace(' ','').rstrip(' ') for item in pep_pos]

pep_overlap = list(set(pep_pos).intersection(set(pep_neg)))   # get the same overlap, 44 as Venny


# update negative

new_df = pd.DataFrame({'peptides':pep_neg})

negative.update(new_df)

# eliminate those false positive
cond = []
for i in range(positive.shape[0]):
    pep = positive.iloc[i]['Mut peptide']
    if pep in pep_overlap:
        cond.append(False)
    else:
        cond.append(True)
positive['cond'] = cond
positive = positive[positive['cond']]
positive = positive.drop(columns=['cond'])


# let's extract the final pep_pos and pep_neg
pep_pos = positive['Mut peptide'].tolist()   # 251
pep_neg = negative['peptides'].tolist()     # 404

# delete the redundant one 
pep_pos = list(set(pep_pos))
pep_neg = list(set(pep_neg))

pep_neg = [item.rstrip(' ') for item in pep_neg]
pep_pos = [item.rstrip(' ') for item in pep_pos]


# let's look at the length distribution
pep_pos_len = [len(item) for item in pep_pos]
pep_neg_len = [len(item) for item in pep_neg]

plt.hist(pep_pos_len,bins=[7,8,9,10,11,12])
plt.hist(pep_neg_len,bins=[7,8,9,10,11,12])


# let's look into for 9mer, the properties
# hydrophobicity, polarity, bulkiness

# first, hydrophobicity
hydro_9mer_pos = []
for i in range(len(pep_pos)):
    if pep_pos_len[i] == 9:
        tmp = get_hydrophobicity(pep_pos[i])
        hydro_9mer_pos.append(tmp)
        
hydro_9mer_neg = []
for i in range(len(pep_neg)):
    if pep_neg_len[i] == 9:
        tmp = get_hydrophobicity(pep_neg[i])
        hydro_9mer_neg.append(tmp)
        
sc.ranksums(hydro_9mer_pos,hydro_9mer_neg)

hydro_10mer_pos = []
for i in range(len(pep_pos)):
    if pep_pos_len[i] == 10:
        tmp = get_hydrophobicity(pep_pos[i])
        hydro_10mer_pos.append(tmp)
        
hydro_10mer_neg = []
for i in range(len(pep_neg)):
    if pep_neg_len[i] == 10:
        tmp = get_hydrophobicity(pep_neg[i])
        hydro_10mer_neg.append(tmp)
        
sc.ranksums(hydro_10mer_pos,hydro_10mer_neg)
sc.ranksums(hydro_10mer_pos+hydro_9mer_pos,hydro_10mer_neg+hydro_9mer_neg)
        
    
# 9mer and 10mer aren't different regarding hydrophobicity
# let's see all the case regardless of length
hydro_pos = []
for i in range(len(pep_pos)):

    tmp = get_hydrophobicity(pep_pos[i])
    hydro_pos.append(tmp)
        
hydro_neg = []
for i in range(len(pep_neg)):
    tmp = get_hydrophobicity(pep_neg[i])
    hydro_neg.append(tmp)
    
sc.ranksums(hydro_pos,hydro_neg)

# let's use correlation to explore the effect of hydrophobicity
dic_pos_freq = {}
for aa in 'ACDEFGHIKLMNPQRSTVWY':
    accum = 0
    for item in pep_pos:
        accum += sum([True if i==aa else False for i in item])
    tmp = accum/sum(pep_pos_len)
    dic_pos_freq[aa] = tmp
    
dic_neg_freq = {}
for aa in 'ACDEFGHIKLMNPQRSTVWY':
    accum = 0
    for item in pep_neg:
        accum += sum([True if i==aa else False for i in item])
    tmp = accum/sum(pep_neg_len)
    dic_neg_freq[aa] = tmp
    
dic_ratio_freq = {}
for aa in 'ACDEFGHIKLMNPQRSTVWY':
    dic_ratio_freq[aa] = dic_pos_freq[aa]/dic_neg_freq[aa]

# let's test if dic_ratio_freq is correlated with hydrophobicity
hydro_index = []
for aa in 'ACDEFGHIKLMNPQRSTVWY':   
    hydro_index.append(properties[aa][0])
    
bulk_index = []
for aa in 'ACDEFGHIKLMNPQRSTVWY':   
    bulk_index.append(properties[aa][1])
    
polar_index = []
for aa in 'ACDEFGHIKLMNPQRSTVWY':   
    polar_index.append(properties[aa][2])

    
sc.spearmanr(list(dic_ratio_freq.values()),hydro_index)
sc.spearmanr(list(dic_ratio_freq.values()),bulk_index)
sc.spearmanr(list(dic_ratio_freq.values()),polar_index)

# using one-sided permutation test to determine if there are any significant differences between pos and neg as to each mm
def my_permutation_test(aa):
    binary_pos = []
    for item in pep_pos:
        for i in item:
            if i == aa: binary_pos.append(1)
            else: binary_pos.append(0)
    binary_neg = []
    for item in pep_neg:
        for i in item:
            if i==aa: binary_neg.append(1)
            else: binary_neg.append(0)
    p_value = permutation_test(binary_pos,binary_neg,method='approximate',num_rounds=10000,func=lambda binary_pos,binary_neg: ((sum(binary_pos)/len(binary_pos))/(sum(binary_neg)/len(binary_neg))),seed=0)                         
    print('Queried aa:',aa)
    print('Observed differences: %.5f' % ((sum(binary_pos)/len(binary_pos))/(sum(binary_neg)/len(binary_neg))))
    print('P value: %.2f' % p_value)  

for aa in 'ACDEFGHIKLMNPQRSTVWY':  
    my_permutation_test(aa)   

'''
significant:
Queried aa: G
Observed differences: 1.34687
P value: 0.01
Queried aa: Q
Observed differences: 1.29419
P value: 0.04
Queried aa: R
Observed differences: 1.24567
P value: 0.04

almost significant:
Queried aa: W
Observed differences: 1.43122
P value: 0.06
Queried aa: V
Observed differences: 1.15401
P value: 0.06
Queried aa: H
Observed differences: 1.24601
P value: 0.08
'''

# position importance by using Kullback-Leibler
# for 9mer
mer9 = []
for i in range(len(pep_pos)):
    if pep_pos_len[i] == 9:
        mer9.append(list(pep_pos[i]))
mat_mer9 = np.array(mer9)
mat_mer9_freq = np.zeros([20,9])
aalist = 'ACDEFGHIKLMNPQRSTVWY'
for i in range(len(aalist)): 
    for j in range(0,9):
        mat_mer9_freq[i,j] = collections.Counter(mat_mer9[:,j])[aalist[i]]/mat_mer9.shape[0]
        
    
mer9_neg = []
for i in range(len(pep_neg)):
    if pep_neg_len[i] == 9:
        mer9_neg.append(list(pep_neg[i]))
mat_mer9_neg = np.array(mer9_neg)
mat_mer9_freq_neg = np.zeros([20,9])
aalist = 'ACDEFGHIKLMNPQRSTVWY'
for i in range(len(aalist)): 
    for j in range(0,9):
        mat_mer9_freq_neg[i,j] = collections.Counter(mat_mer9_neg[:,j])[aalist[i]]/mat_mer9_neg.shape[0]       
        
    
# after having two frequency matrix, let's calculate KL divergence
importance = np.zeros([9,])
for j in range(mat_mer9_freq.shape[1]):
    #importance[j] = sc.entropy(mat_mer9_freq[:,j],mat_mer9_freq_neg[:,j]) 
    importance[j] = KL_divergence(mat_mer9_freq[:,j],mat_mer9_freq_neg[:,j])   
'''
array([0.13706572, 0.28142314, 0.18352944, 0.19483134, 0.20482319,
       0.1238152 , 0.11597121, 0.13078297, 0.11547343])

'''               

# for 10mer
pep_neg = [item.rstrip(' ') for item in pep_neg]
pep_pos = [item.rstrip(' ') for item in pep_pos]
pep_pos_len = [len(item) for item in pep_pos]
pep_neg_len = [len(item) for item in pep_neg]

mer10 = []
for i in range(len(pep_pos)):
    if pep_pos_len[i] == 10:
        print(pep_pos[i])
        mer10.append(list(pep_pos[i]))
mat_mer10 = np.array(mer10)
mat_mer10_freq = np.zeros([20,10])
aalist = 'ACDEFGHIKLMNPQRSTVWY'
for i in range(len(aalist)): 
    for j in range(0,10):
        mat_mer10_freq[i,j] = collections.Counter(mat_mer10[:,j])[aalist[i]]/mat_mer10.shape[0]
        
    
mer10_neg = []
for i in range(len(pep_neg)):
    if pep_neg_len[i] == 10:
        mer10_neg.append(list(pep_neg[i]))
mat_mer10_neg = np.array(mer10_neg)
mat_mer10_freq_neg = np.zeros([20,10])
aalist = 'ACDEFGHIKLMNPQRSTVWY'
for i in range(len(aalist)): 
    for j in range(0,10):
        mat_mer10_freq_neg[i,j] = collections.Counter(mat_mer10_neg[:,j])[aalist[i]]/mat_mer10_neg.shape[0]       
        
    
# after having two frequency matrix, let's calculate KL divergence
importance10 = np.zeros([10,])
for j in range(mat_mer10_freq.shape[1]):
    #importance[j] = sc.entropy(mat_mer9_freq[:,j],mat_mer9_freq_neg[:,j]) 
    importance10[j] = KL_divergence(mat_mer10_freq[:,j],mat_mer10_freq_neg[:,j])      
    
'''
array([0.2158908 , 0.33372127, 0.20192451, 0.41047172, 0.2295425 ,
       0.20682632, 0.2817167 , 0.14172125, 0.14058824, 0.04602353])
'''    

        































        