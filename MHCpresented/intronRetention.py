#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 20:39:59 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop/project')
import sys
import pandas as pd
import re
import requests
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from urllib.error import HTTPError

def intron(event,EnsGID,junction,dict_exonCoords,dictExonList):
    merBucket = []
    if event.startswith('E'):   # only consider this situation, since intron retentions are predominantly associated with a translatable preceding exon, ending up with a early stop codon
        # namely: E2.4-I2.1, E22.1-I22.1
        former = event.split('-')[0]  # E2.4, E22.1
        attrs = dict_exonCoords[EnsGID][former] # chr, strand, start, end
        strand = attrs[1]        
        allTransDict = dictExonList[EnsGID] 
        for tran,exonlist in allTransDict.items():
            if former in exonlist: 
                try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                else:
                    if type(tranStartIndex) == int:
                        if strand == '+':
                            intronStartIndex = int(attrs[3]) + 1
                            remainder = (intronStartIndex - tranStartIndex) % 3   # how many nts remaining before the first nt in intron
                            # if 0, means the first nt in intron will be the first nt in codon triplet,
                            # if 1, means the first nt in intron will be the second nt in codon triplet,
                            # if 2, means the first nt in intron will be the third nt in codon triplet.
                            if remainder == 0: junctionSlice = junction[junction.find(',') - 24:].replace(',','')
                            elif remainder == 1: junctionSlice = junction[junction.find(',') - 25:].replace(',','')
                            elif remainder == 2: junctionSlice = junction[junction.find(',') - 26:].replace(',','')
                        elif strand == '-':
                            intronStartIndex = int(attrs[3]) - 1
                            remainder = (tranStartIndex - intronStartIndex) % 3 
                            if remainder == 0: junctionSlice = junction[junction.find(',') + 24:].replace(',','')
                            elif remainder == 1: junctionSlice = junction[junction.find(',') + 25:].replace(',','')
                            elif remainder == 2: junctionSlice = junction[junction.find(',') + 26:].replace(',','')
                
                        peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                        merArray = extractNmer(peptide,9) 
                        merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
        
    elif event.startswith('I'): merBucket = [[]]
    return merBucket

def extractNmer(peptide,N):
    starIndex = peptide.find('*') 
   # print(starIndex)
    merArray = []
    if starIndex == -1 and len(peptide) >= N:
        for i in range(0,len(peptide)-N+1,1):
            mer = peptide[i:i+N]
            merArray.append(mer)
    if starIndex >= N:
        peptideTrun = peptide[:starIndex]
        #print(peptideTrun)
        for j in range(0,len(peptideTrun)-N+1,1):
            mer = peptideTrun[j:j+N]
            merArray.append(mer)
    return merArray



# https://rest.ensembl.org/documentation/info/lookup
def grabEnsemblTranscriptTable(EnsTID):
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/{0}?expand=1".format(EnsTID)     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})     
    decoded = r.json()
    try: translationStartIndex = decoded['Translation']['start']
    except KeyError: print('{0} is not translatable'.format(EnsTID))
    else: return translationStartIndex
# for except branch, we don't specify return command, so it will return NoneType    
    
    
    
    
    

    
    
    
    
    