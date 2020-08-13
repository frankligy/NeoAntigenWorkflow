#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:25:58 2020

@author: ligk2e
"""
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SubsMat import MatrixInfo
import json
import numpy as np
from collections import Counter
import pandas as pd

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split,Subset


def json2fsa(hla):
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope/{}.json'.format(hla),'r') as f:
        data = json.load(f)
        
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/{}.fsa'.format(hla),'w') as f:        
        for item in data:
            key = list(item.keys())[0]
            value = list(item.values())[0]
            f.write('>{}\n'.format(key))
            f.write('{}\n'.format(value))

def multiple_json2fsa():
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope/inventory_new.txt','r') as f1:
        for line in f1:
            line = line.rstrip('\n')
            json2fsa(line)
'''          
run clustal-omega:
    download binary: http://www.clustal.org/omega/
    chmod 777 ./clustal-omega-1.2.3-macosx
    ./clustal-omega-1.2.3-macosx -i "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/HLA-A*0101.fsa" -o "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/HLA-A*0101.aligned.fasta" --auto -v            
   
    
run multiple sequentially:
    cat /Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope/inventory_new.txt | while read line; do
    ./clustal-omega-1.2.3-macosx -i "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/${line}.fsa" -o "/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/${line}.aligned.fasta" --auto -v; done         
    
only hla that has more than 1 paratope will be processed in clustal-omega   
'''

def single_paratope(hla):
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/{}.fsa'.format(hla),'r') as f:
        seq = f.readlines()[1].rstrip('\n')
    return hla,seq


def matrix2concensus(mat):
    final = ''
    for j in range(mat.shape[1]):
        most = Counter(mat[:,j]).most_common(1)[0][0]  # if most_common(2): [('A', 3), ('C', 1)]
        if most == '-':
            most = Counter(mat[:,j]).most_common(2)[1][0]
        final += most
    return final
    
    

def msa_paratope(hla):
    alignment = AlignIO.read(open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/{}.aligned.fasta'.format(hla)),"fasta")
    msa = []
    for record in alignment:
        msa.append(list(record.seq))   # another part is record.id
    mat = np.array(msa)
    final = matrix2concensus(mat)
    return hla,final
    
    
        
    

def hla_paratope():

            
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa/inventory_single.txt','r') as single:
        singles = single.readlines()  # each one will contain '\n'
        singles = [item.rstrip('\n') for item in singles]
    
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_paratope_fsa_aligned/inventory_msa.txt','r') as multiple:
        multiples = multiple.readlines()
        multiples = [item.rstrip('\n') for item in multiples]
        
    
    with open('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable.txt','w') as f:
        for item in singles:
            hla,seq = single_paratope(item)   
            f.write('{0}\t{1}\n'.format(hla,seq))
        for item in multiples:
            hla,seq = msa_paratope(item)
            f.write('{0}\t{1}\n'.format(hla,seq))
            
########################################################################################################            

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

def convert_hla(hla):
    cond = True
    hla = hla.replace(':','')
    if len(hla) < 9: cond = False   # HLA-A3
    elif len(hla) == 9:   # HLA-A3002
        f = hla[0:5]  # HLA-A
        e = hla[5:]   # 3002
        hla = f+'*'+e
    return hla,cond

def convert_hla_series(df):
    new = []
    col = []
    for i in df['HLA']:
        hla,cond = convert_hla(i)
        col.append(cond)
        if cond == True: new.append(hla)
    df = df.loc[pd.Series(col)]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    df['HLA'] = new
    return df
        
def test_no_space(series):
    for i in series:
        if ' ' in i:
            print('damn')

   
'''
a = pd.read_excel('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data.xlsx')
a1 = clean_data_frame(a)
test_no_space(a1.iloc[:,0])
test_no_space(a1.iloc[:,1])
a1.iloc[:,2].dtype

a2 = convert_hla_series(a1)
a2.to_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data.txt',sep='\t',index=None)

then use:
    
{ cat data.txt | head -n 1; cat data.txt | tail -n +2 | sort -u -k1,2; } > data_new.txt, only 32669 training data left

ori = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data_new.txt',sep='\t')
hla = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable.txt',sep='\t',header=None,names=['hla','paratope'])
inventory = hla['hla']
dic_inventory = dict_inventory(inventory)


'''
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
        self.new = self.padding()
        
    def __len__(self):
        return len(self.new)
    
    def __getitem__(self,idx):
        return self.new[idx]
    
    
    def padding(self):
        len_values = [tup[0].shape[0] for tup in self.middle]
        max_length = max(len_values)
        
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
    
    def paratope_dic(self):
        df = self.hla
        self.dic = {}
        for i in range(df.shape[0]):
            hla = df['hla'].iloc[i]
            paratope = df['paratope'].iloc[i]
            self.dic[hla] = paratope
    
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
            print(i)
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
            X = dataset.blosum50(cat).float()   # 2-d tensor
            y = torch.tensor(immuno).long()  # 0-d tensor
            lis.append((X,y))
        return lis
            
            
            
            

    



















        
    
    