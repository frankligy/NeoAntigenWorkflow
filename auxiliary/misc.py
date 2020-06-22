#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 20:00:16 2020

@author: ligk2e
"""

import os
import pandas as pd
import numpy as np

counts = pd.read_csv('/Users/ligk2e/Desktop/project_breast/R1-V6/counts.TCGA-BRCA.txt',sep='\t')

events = pd.read_csv('/Users/ligk2e/Desktop/TCGA-breast/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-filtered-names-75p.txt',sep='\t')

tumorPSI = events['TCGA-E2-A10A-01'].values

tumorPSI[np.isnan(tumorPSI)] = 0.0

tumorPSI = tumorPSI.tolist()


UID = events['UID'].tolist()

part = []
for i in range(len(UID)):
    uid = UID[i]
    x = uid.split('|')
    try: x[0].split(':')[3]
    except IndexError: event = x[0].split(':')[2]
    else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])
    query = uid.split('|')[1].split(':')[0] + ':' + event
    part.append(query)
    
counts_id = counts['AltAnalyze_ID'].tolist()

counts_id_useful = [item.split('=')[0] for item in counts_id]


new = [name[:15] for name in counts.columns]
counts.columns = new


query_id_counts = counts['TCGA-E2-A10A-01'].tolist()



    
dic = {}
for i in range(len(counts_id_useful)):
    ensg = counts_id_useful[i].split(':')[0]
    event = counts_id_useful[i].split(':')[1:]   
    event = ''.join(event)
    info = query_id_counts[i]
    try:
        dic[ensg].append((event,info))
    except KeyError:
        dic[ensg] = []
        dic[ensg].append((event,info))



counts_extract = []
for i in range(len(part)):
    flag = False
    part_ensg = part[i].split(':')[0]
    searchSpace = dic[part_ensg]
    for j in searchSpace:
        match = part_ensg+':'+j[0]
        if part[i] == match:
            counts_extract.append(j[1])
            flag = True
            break
    if flag == False:
        counts_extract.append(0)


final = pd.DataFrame({'UID':UID,'PSI':tumorPSI,'count':counts_extract})

cond = []
for i in range(final.shape[0]):
    psi = final['PSI'].iloc[i]
    if psi == 0.0:
        cond.append(False)
    else:
        cond.append(True)
final['cond'] = cond
final = final[final['cond']]
final = final.drop(columns=['cond'])


counts_h = pd.read_csv('/Users/ligk2e/Desktop/project_breast/R1-V6/counts.TCGA-BRCA.healthy.txt',sep='\t')
def check(key,counts_h):  # check the all counts for a given event, counts_h contains all counts for matched control
    for i in range(counts_h.shape[0]):
        event = counts_h.iloc[i]['0']
        if key in event:
            detail = counts_h.iloc[i].tolist()[1:]
            average = np.mean(detail)
            if average <= 10: cond = True
            else: cond = False
            break
    return cond

col = []
for i in range(final.shape[0]):
    print(i,check(event,counts_h))
    event = ':'.join(final['UID'].iloc[i].split('|')[0].split(':')[1:])
    col.append(check(event,counts_h))
final['col'] = cond
final = final[final['col']]
final = final.drop(columns=['col'])



final.to_csv('/Users/ligk2e/Desktop/TCGA-breast/TCGA-E2-A10A-01.ori.txt',sep='\t',index=None)







def grabEnsemblTranscriptTable1(EnsTID):
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/{0}?expand=1".format(EnsTID)     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})     
    try: decoded = r.json()
    except: 
        print('JSON unknoen error')
    else:
        try: translationStartIndex = decoded['Translation']['start']
        except KeyError: print('{0} is not translatable'.format(EnsTID)) # might be invalid transcriptID or they don't have tranStartIndex(don't encode protein)
        else: return translationStartIndex

import requests
grabEnsemblTranscriptTable('ENST00000420190')



gtfEnsembl91 = pd.read_csv('/Users/ligk2e/Downloads/Homo_sapiens.GRCh38.100.chr.gtf',sep='\t',skiprows=5,header=None)
newName = ['seqname','source','feature','start','end','score','strand','frame','attribute']
gtfEnsembl91.columns = newName

cond = []
for i in range(gtfEnsembl91.shape[0]):
    print(i)
    anno = gtfEnsembl91.iloc[i]['feature']
    if anno == 'start_codon' or anno == 'stop_codon':
        cond.append(True)
    else:
        cond.append(False)
        
gtfEnsembl91['cond'] = cond
gtfEnsembl91_retained = gtfEnsembl91[cond]  
gtfEnsembl91 = gtfEnsembl91_retained.drop(columns=['cond'])

tranID = []
for i in range(gtfEnsembl91.shape[0]):
    print(i)
    attr = gtfEnsembl91.iloc[i]['attribute']
    attr_list = attr.split(';')
    for item in attr_list:
        if item.startswith(' transcript_id'):   # the preceding space!!
            enstID = item.split('"')[1].rstrip('"')
            tranID.append(enstID)
gtfEnsembl91['tranID'] = tranID

gtfEnsembl91.to_csv('/Users/ligk2e/Desktop/TCGA-breast/gtfEnsembl91.txt',sep='\t',index =None)
    


def dictGTFconstruct():
    global dicGTF
    dicGTF = {}
    for i in range(gtfEnsembl91.shape[0]):
        feature = gtfEnsembl91.iloc[i]['feature']
        start = gtfEnsembl91.iloc[i]['start']
        tranID = gtfEnsembl91.iloc[i]['tranID']
        if feature == 'start_codon':
            dicGTF[tranID] = start

dictGTFconstruct()
    
def grabEnsemblTranscriptTable(EnsTID):
    try:
        translationStartIndex = dicGTF[EnsTID]
    except KeyError:
        print('{0} does not have valid EnsTID'.format(EnsTID))
    else: return translationStartIndex
        
        
    
    





























