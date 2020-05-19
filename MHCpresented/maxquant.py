#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 10:09:24 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop/project_test')
import pandas as pd
import ast

neo_antigen11 = pd.read_csv('NeoJunction_11_mark.txt',sep='\t',usecols=['MHCIresult'])
neo_antigen_list11 = neo_antigen11['MHCIresult'].tolist()

def getb(list_):
    sb,wb = [],[]
    for item in list_:
        if item == 'No candidates': continue
        else: 
            info = ast.literal_eval(item)
            pep = info['HLA-A29:02']
            pep_sb,pep_wb = pep[0],pep[1]
            sb.extend(pep_sb)
            wb.extend(pep_wb)
    b = sb + wb
    return b

b11 = getb(neo_antigen_list11)
b11_st = set(b11)

# write to fasta
def toFasta(list_):
    with open('./resultMHC/whole.fa','w') as file1:
        for index,item in enumerate(list_):
            file1.write('>mer{0}\n'.format(index+1))
            file1.write('{0}\n'.format(item))
            
toFasta(b)
        
        
junction = pd.read_csv('./resultMHC/NeoJunction_11_new.txt',sep='\t',usecols=['exam_seq'])   
junction_intact = [item.replace(',','') for item in junction['exam_seq'].tolist()]

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


bucket = []
for item in junction_intact:
    item1,item2,item3 = item[0:],item[1:],item[2:]
    for i in [item1,item2,item3]:
        try: AA = str(Seq(i,generic_dna).translate(to_stop=False))
        except: print('preexisting bug!!!! subexon doesn\'t exist')
        else:
            AA_list = AA.split('*')
            for aa in AA_list:
                if len(aa) >= 8:bucket.append(aa)

toFasta(bucket)

counter = 0
for i in b:
    for j in bucket:
        if i in j: 
            counter += 1
            break

# we prove that all the 11 mer are in whole.fasta file





peptide = pd.read_csv('/Users/ligk2e/Downloads/peptide.txt',sep='\t',header=None,names=['peptide'])
experiment = peptide['peptide'].tolist()

experiment_st = set(experiment)

b3 = b_st.intersection(experiment_st)
b4 = list(b3)

count = 0
for item in b:
    for pep in experiment:
        if item in pep: count += 1


mer11 = list(filter(lambda x: len(x)==11,experiment))

for item in b:
    for rr in mer11:
        if item == rr: print(item)


peptide_paper = pd.read_csv('/Users/ligk2e/Downloads/peptide_paper.txt',sep='\t',header=None,names=['peptide'])
experiment_paper = peptide_paper['peptide'].tolist()

experiment_paper_st = set(experiment_paper)

b5 = list(experiment_st.intersection(experiment_paper_st))

b6 = list(b_st.intersection(experiment_paper_st))












