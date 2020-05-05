#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:32:27 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop/project_breast/R1-V6/')
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from decimal import Decimal as D
import pickle
import regex
import re
from time import process_time
import collections
from urllib.error import HTTPError
import requests
#from backgroundGTEx import backgroundGTEx,convertExonList,merPassGTExCheck,spread,getNmerNormal
#from intronRetention import intron,grabEnsemblTranscriptTable


class Meta():  #inspect an object: dir(), vars(), instanceName.__dict__, mannually set __repr__ or __str__
    def __init__(self, df):
        self.df = df
    
    def __repr__(self):
        return {'df':self.df}
    
    def retrieveJunctionSite(self):
        exam_seq,back_seq = [],[]
        for i in range(self.df.shape[0]):
            temp = uid(self.df,i)
            EnsID = list(temp.keys())[0].split(':')[1]
            exam_site = list(temp.values())[0][0]
            back_site = list(temp.values())[0][1]
            exam_site_1 = subexon_tran(exam_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
            exam_site_2 = subexon_tran(exam_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
            exam_seq_join = ','.join([exam_site_1,exam_site_2])
            exam_seq.append(exam_seq_join)
            back_site_1 = subexon_tran(back_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
            back_site_2 = subexon_tran(back_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
            back_seq_join = ','.join([back_site_1,back_site_2])
            back_seq.append(back_seq_join)
            
        self.df['exam_seq'] = exam_seq
        self.df['back_seq'] = back_seq
        
    def matchWithExonlist(self,df_exonlist,dict_exonCoords):
        col1,col2 = [],[]        
        for i in range(self.df.shape[0]):
            temp=uid(self.df,i)
            EnsID=list(temp.keys())[0].split(':')[1]
            Exons_examined = exon_extract(temp,0,EnsID)
            Exons_back = exon_extract(temp,1,EnsID)
            col1.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_examined))
            col2.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_back))
        self.df['exam_whole_transcripts'] = col1
        self.df['back_whole_transcripts'] = col2

    def getORF(self):
        col0,col1 = [],[]
        for index,value in enumerate(['exam_whole_transcripts','back_whole_transcripts']): 
            col = eval('col'+str(index))
            for transcripts in list(self.df[value]):
                if not transcripts: col.append(None)  # None's type is Nonetype
                tempArray = []
                for transcript in transcripts:
                    if not transcript: tempArray.append('')
                    else:
                        maxTran = transcript2peptide(transcript)
                        tempArray.append(maxTran)
                col.append(tempArray)
        self.df['exam_ORF_tran'] = col0
        self.df['back_ORF_tran'] = col1
        
    def getORFaa(self):
        col0,col1 = [],[]
        for index,value in enumerate(['exam_ORF_tran','back_ORF_tran']):
            col = eval('col'+str(index))
            for transcripts in list(self.df[value]):
                if not transcripts: col.append(None)
                tempArray = []
                for transcript in transcripts:
                    if not transcript: tempArray.append('')
                    else:
                        maxAA = str(Seq(transcript,generic_dna).translate(to_stop=False))
                        tempArray.append(maxAA)
                col.append(tempArray)
        self.df['exam_ORF_aa'] = col0
        self.df['bacK_ORF_aa'] = col1
                    

                  
                
            
            

class NeoJ(Meta):   
    def __init__(self,df,N):
        super().__init__(df)   # super() means Meta
        self.mer = N    # length of mer that we are aiming to get

    def phaseTranslateJunction(self):
        col0,col1 = [],[]
        for i in range(self.df.shape[0]):
            wholeTrans = list(self.df['exam_whole_transcripts'])[i]
            ORFtrans = list(self.df['exam_ORF_tran'])[i]
            junctionOri = list(self.df['exam_seq'])[i]  # have comma to delimit former and latter part
            junction = list(self.df['exam_seq'])[i].replace(',','')
            phaseArray,peptideArray = [],[]
            if not wholeTrans:             # the EnsGID doesn't exist in mRNA-exon file
                phaseArray.append('#')
                peptideArray.append('#')
                print('This gene {0} is absent in mRNA-ExonID file\n'.format(list(self.df['UID'])[i]))
            for j in range(len(wholeTrans)):
                if not wholeTrans[j]: 
                    phaseArray.append('') 
                    peptideArray.append('') # splicing evnet can not match to this transcript
                else:
                    startJun = wholeTrans[j].find(junction[:-1]) # trim the right-most overhang for simplicity
                    endJun = startJun + len(junction[:-1])
                    startORF = wholeTrans[j].find(ORFtrans[j])
                    endORF = startORF + len(ORFtrans[j])
                    if startJun > endORF or endJun < startORF:
                        phase = '*'  # means junction site will not be translated
                        peptide = '*'
                        phaseArray.append(phase)
                        peptideArray.append(phase)
                        print('The {}th transcript of {}, even though junction site\
                        could match with, but optimal ORF suggests that junction site will not be translated'.format(j,list(self.df['UID'])[i]))
                    else: 
                        former = junctionOri.find(',')  # length of former part, also the index of splicesite
                        phase = abs(startJun + former - startORF) % 3 
                        N = self.mer
                         
                        if phase == 0: 
                            front=0 if former - ((N-1)*3) < 0 else former - (N-1)*3
                            
                            junctionSlice = junction[front: former + ((N-1)*3)].replace(',','')
                            print(junctionSlice)
                        elif phase == 1: 
                            front=0 if former - ((N-1)*3+1) <0 else former - ((N-1)*3+1)
                            junctionSlice = junction[front: former + ((N-1)*3+2)].replace(',','')
                            
                        elif phase == 2: 
                            front=0 if former- ((N-1)*3+2) <0 else former- ((N-1)*3+2)
                            junctionSlice = junction[front: former + ((N-1)*3+1)].replace(',','')
                 
            
                        peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                        print(peptide)
                    
                        
                        
                        
                        
                        

                        phaseArray.append(phase)
                        peptideArray.append(peptide)
            
            col0.append(phaseArray) 
            col1.append(peptideArray)
        self.df['phase'] = col0 # ['',1,'*',2] OR ['#']
        self.df['peptide'] = col1  # ['',MSTTG,'*',TYCTT] OR ['#']
                
    def getNmer(self):
        col = []
        for i in range(self.df.shape[0]):
            peptides = list(self.df['peptide'])[i]
            merArray = []
            condition = any([False if pep=='' else True for pep in peptides]) # if False, means ['','',''], means can not match with any of existing transcript
            if not condition: merArray.append('MANNUAL')
            else:
                for peptide in peptides:
                    if peptide == '': merArray.append('')
                    elif peptide == '*': merArray.append('*')
                    else: 
                        tempArray = extractNmer(peptide,self.mer)
                        merArray.append(tempArray)
            truthTable1 = [False if mer == '*' or mer == [] else True for mer in merArray]
            if not any(truthTable1): merArray = ['Match but no translation or early stop']
            col.append(merArray)  # [ '','*',[], ['MTDJDKFDJF','FJDKFJDF'] ]
        self.df['{0}mer'.format(self.mer)] = col
        
    def getFinalMers(self):
        col = []
        for i in range(self.df.shape[0]):
            uid, merArray = list(self.df['UID'])[i], list(self.df['{0}mer'.format(self.mer)])[i]
            col.append(merPassGTExCheck(dictGTEx,uid,merArray))
        self.df['matchedFinalMer'] = col
                 
    def mannual(self,check=False):
        col = []
        for i in range(self.df.shape[0]):
            merArray = []
            uid,junction,Nmer = self.df['UID'].tolist()[i],self.df['exam_seq'].tolist()[i],self.df['{0}mer'.format(self.mer)].tolist()[i]
            if Nmer == ['MANNUAL']: 
                EnsGID = uid.split('|')[0].split(':')[1]
                x = uid.split('|')
                try: x[0].split(':')[3]   # these try-excepts aim to consider fusion gene situation
                except IndexError: event = x[0].split(':')[2]
                else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])     # E22-ENSG:E31  
                try: x[1].split(':')[2]
                except IndexError: backEvent = x[1].split(':')[1]
                else: backEvent = x[1].split(':')[2]                         
                if re.search(r'I.+_\d+',event) or 'U' in event:  # they belong to NewExon type: including UTR type and blank type(the annotation is blank)
                # for NewExon type, no need to infer translation phase, just use junction sequence, mer should span junction site
                    junctionIndex = junction.find(',')
                    Nminus1 = self.mer - 1 
                    junctionSlice = junction[junctionIndex-((Nminus1*3)-1):junctionIndex+(Nminus1*3)].replace(',','') # no need to worry out of index, slicing operation will automatically change it to 'end' if overflow
                    merArray = dna2aa2mer(junctionSlice,self.mer)   # merArray is a nested list
                    if check == True: merArray = merPassGTExCheck(dictGTEx,uid,merTumor)
                if re.search(r'E\d+\.\d+_\d+',event) or 'ENSG' in event:  # they belong to newSplicing site or fusion gene
                # for this type, just check if the former part (E1.2-E3.4, here E1.2 is former part) exist in known transcript, then infer the phase of translation
                    print(event)
                    merArray = newSplicingSite(event,EnsGID,junction,self.mer,dictExonList,dict_exonCoords)   # nested list
                    if check==True: merArray = merPassGTExCheck(dictGTEx,uid,merTumor)
                    if merArray == [[]]: merArray = ['newSplicing or fusion gene, already checked, query subexon not present in known transcript']
                if re.search(r'^I\d+\.\d+-',event) or re.search(r'-I\d+\.\d+$',event):
                    merTumor = intron(event,EnsGID,junction,dict_exonCoords,dictExonList,self.mer)  # nested list
                    if check==True: merArray = merPassGTExCheck(dictGTEx,uid,merTumor)
                    if merArray == [[]]: merArray = ['intron retention, already checked, either former subexon not present in known transcript or matched transcript is not Ensembl-compatible ID']
                if re.search(r'E\d+\.\d+-E\d+\.\d+$',event):   # novel splicing: CLASP1:ENSG00000074054:E25.1-E26.1, just don't match with any existing one
                    # here we extract the backEvent, which is E24.1-E27.1
                    print(event)
                    backEvent = backEvent.replace('-','|')  # now it is E24.1|E27.1
                    merArray = novelOrdinal(event,backEvent,EnsGID,junction,dictExonList,dict_exonCoords,self.mer)
                    if check==True: merArray = merPassGTExCheck(dictGTEx,uid,merTumor)
                    if merArray ==[[]]: merArray = ['novel ordinal splicing event, already checked, its background event can not match up with any existing transcript either']
                    
            col.append(merArray)
        self.df['mannual'] = col
        
    def collectNeoAntigen(self):
        list_ = []
        known,novel = self.df['matchedFinalMer'].tolist(),self.df['mannual'].tolist()
        whole = known+novel
        filterBucket = [[],['newSplicing or fusion gene, already checked, either not translated or be eliminated by background check'],['intron retention, already checked, either former subexon not present in known transcript or matched transcript is not Ensembl-compatible ID']]
        for i in whole:
            if not i in filterBucket: list_ += i
        self.mhcNeoAntigens = list(set(list_))
        
def novelOrdinal(event,backEvent,EnsGID,junction,dictExonList,dict_exonCoords,N): # E25.1-E26.1
    latter = event.split('-')[1]  # E26.1
    attrs = dict_exonCoords[EnsGID][latter] #  chr, strand, start, end
    strand,start = attrs[1], attrs[2]
    merBucket = []
    allTransDict = dictExonList[EnsGID]
    for tran,exonlist in allTransDict.items():
        if backEvent in exonlist:
            try:tranStartIndex = grabEnsemblTranscriptTable(tran)
            except HTTPError: continue
            else:
                if type(tranStartIndex) == int:
                     if strand == '+':
                         remainder = (int(start) - tranStartIndex) % 3   
                         if remainder == 0: 
                             front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                         elif remainder == 1: 
                             front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                         elif remainder == 2: 
                             front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                             junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                     elif strand == '-':

                         remainder = (tranStartIndex - int(start)) % 3 
                         if remainder == 0: 
                             front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                         elif remainder == 1: 
                             front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                             junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                         elif remainder == 2: 
                             front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                             junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                
                     peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                     merArray = extractNmer(peptide,N) 
                     merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
           
    return merBucket
                            
                    

def newSplicingSite(event,EnsGID,junction,N,dictExonList,dict_exonCoords):   # one stop solution for newSplicingSite type
    try: event.split('-')[0].split('_')[1]   # see if the former one has trailing part E1.2_8483494
    except IndexError: 
        former = event.split('-')[0]  # no trailing part, former part just E1.2, means newSplicingSite is in latter part
        query = event.split('-')[1].split('_')[0]
        trailing = event.split('-')[1].split('_')[1]
    else: 
        query = event.split('-')[0].split('_')[0]  # has trailing part, former part should get rid of trailing part
        trailing = event.split('-')[0].split('_')[1]
    finally:
        merBucket = []   
        attrs = dict_exonCoords[EnsGID][query] # chr, strand, start, end
        strand = attrs[1]        
        allTransDict = dictExonList[EnsGID] 
        for tran,exonlist in allTransDict.items():
            if query in exonlist: 
                try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                else:
                    if type(tranStartIndex) == int:
                        if strand == '+':
                            remainder = (int(trailing) - tranStartIndex) % 3   
                            if remainder == 0: 
                                front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                        elif strand == '-':

                            remainder = (tranStartIndex - int(trailing)) % 3 
                            if remainder == 0: 
                                front = 0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front = 0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[junction.find(',') - ((N-1)*3+2):junction.find(',')+((N-1)*3+1)].replace(',','')
                
                        peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                        merArray = extractNmer(peptide,N) 
                        merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
           
    return merBucket

def getPhasePeptide(EnsGID,wholeTrans,ORFtrans,partialJunction,fullJunction):   # used in newSplicingSite situation
    phaseArray,peptideArray = [],[]
    if not wholeTrans:             # the EnsGID doesn't exist in mRNA-exon file
        phaseArray.append('#')
        peptideArray.append('#')
        print('This gene {0} is absent in mRNA-ExonID file\n'.format(EnsGID))
    for j in range(len(wholeTrans)):
        if not wholeTrans[j]: 
             phaseArray.append('') 
             peptideArray.append('')# splicing evnet can not match to this transcript
        else:
            startJun = wholeTrans[j].find(partialJunction[:-1]) # trim the right-most overhang for simplicity
            endJun = startJun + len(partialJunction[:-1])
            startORF = wholeTrans[j].find(ORFtrans[j])
            endORF = startORF + len(ORFtrans[j])
            if startJun > endORF or endJun < startORF:
                phase = '*'  # means junction site will not be translated
                peptide = '*'
                phaseArray.append(phase)
                peptideArray.append(phase)
                print('The {}th transcript of {}, even though junction site\
                could match with, but optimal ORF suggests that junction site will not be translated'.format(j,EnsGID))
            else: 
                
                
        
                former = fullJunction.find(',')  # length of former part, also the index of splicesite
                phase = abs(startJun + former - startORF) % 3 
                N = self.mer
                 
                if phase == 0: junctionSlice = junction[former - ((N-1)*3): former + ((N-1)*3)].replace(',','')
                elif phase == 1: junctionSlice = junction[former - ((N-1)*3+1):former + ((N-1)*3+2)].replace(',','')
                elif phase == 2: junctionSlice = junction[former- ((N-1)*3+2): former + ((N-1)*3+1)].replace(',','')
         
    
                peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                
                phaseArray.append(phase)
                peptideArray.append(peptide)
    return phaseArray,peptideArray

                    
def dna2aa2mer(dna,N):
    manner1,manner2,manner3 = dna[0:],dna[1:],dna[2:]
    merBucket = []
    for manner in [manner1,manner2,manner3]:
        #print(manner)
        try:
            aa = str(Seq(manner,generic_dna).translate(to_stop=False))
        except: 
            print('There are *** in the junction site, previous bugs')
        else:    
            mer = extractNmer(aa,N)
            merBucket.append(mer)
    if merBucket == []: merBucket = [[]]
    return merBucket # a nested list of 9mer
              
                
             
def extractNmer(peptide,N):  # it already considers '' and '*'
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
                
        
        

def transcript2peptide(cdna_sequence):   # actually to ORF
    reading_manners = []
    reading_manners.append(cdna_sequence[0:])
    reading_manners.append(cdna_sequence[1:])
    reading_manners.append(cdna_sequence[2:])
    frag_comp_array = []
    for manner in reading_manners:       
        pos = []
        for m in re.finditer(r'(TAA|TGA|TAG)',manner):   # for multiple instances
            if m.start() % 3 == 0:
                pos.append(m.start())
        frag_array = pos_to_frags(pos,manner)
        for frag in frag_array:
            if 'ATG' not in frag or len(frag) == 0:
                continue
            else:
                for n in re.finditer(r'ATG',frag):
                    if (len(frag) - n.start()) % 3 == 0:
                        frag_comp = frag[n.start():]
                        frag_comp_array.append(frag_comp)
                        break   # might have multiple 'ATG' so it is necessary to break when find first 'ATG'
                    else:
                        continue
                    
    #######################  # We think if you only has longer length(0-7) but add_score is not higher than original one, you are FAlSE
    max_seq = ''
    max_length = 0
    max_item_score = 0
    for item in frag_comp_array:
        temp1 = len(item)
        add_score = score_GC(item) + score_coding_bias(item)
        if (temp1 - max_length) >= 8:
            max_length = temp1
            max_item_score = add_score
            max_seq = item
        elif (temp1 - max_length) >= 0 and (temp1 - max_length) < 8:
            if add_score >= max_item_score:
                max_length = temp1
                max_item_score = add_score
                max_seq = item
#           else:
#                print('equal length but less likely to be a true ORF or longer length but less likely to be a true ORF',add_score,max_item_score) 
    max_seq_tran = max_seq
    return max_seq_tran

def score_GC(sequence):
    GC_content = 0
    length_seq = len(sequence)
    for nt in sequence:
        if nt == 'G' or nt == 'C':
            GC_content += 1
    GC_percent = GC_content / length_seq
    return GC_percent
            
def score_coding_bias(sequence):
    # coding frequency table is from GenScript webpage
    usage_dict = {'TTT':16.9,'TTC':20.4,'TTA':7.2,'TTG':12.6,'TAT':12.0,'TAC':15.6,'TAA':0.7,'TAG':0.5,
                  'CTT':12.8,'CTC':19.4,'CTA':6.9,'CTG':40.3,'CAT':10.4,'CAC':14.9,'CAA':11.8,'CAG':34.6,
                  'ATT':15.7,'ATC':21.4,'ATA':7.1,'ATG':22.3,'AAT':16.7,'AAC':19.5,'AAA':24.0,'AAG':32.9,
                  'GTT':10.9,'GTC':14.6,'GTA':7.0,'GTG':28.9,'GAT':22.3,'GAC':26.0,'GAA':29.0,'GAG':40.8,
                  'TCT':14.6,'TCC':17.4,'TCA':11.7,'TCG':4.5,'TGT':9.9,'TGC':12.2,'TGA':1.3,'TGG':12.8,
                  'CCT':17.3,'CCC':20.0,'CCA':16.7,'CCG':7.0,'CGT':4.7,'CGC':10.9,'CGA':6.3,'CGG':11.9,
                  'ACT':12.8,'ACC':19.2,'ACA':14.8,'ACG':6.2,'AGT':11.9,'AGC':19.4,'AGA':11.5,'AGG':11.4,
                  'GCT':18.6,'GCC':28.5,'GCA':16.0,'GCG':7.6,'GGT':10.8,'GGC':22.8,'GGA':16.3,'GGG':16.4} 
    # do a normaliztion for each triplet, then for all the triplet's sum, divided by the number of triplet
    min_freq = 4.5
    max_freq = 40.8
    norm_usage_dict = {}
    for codon,freq in usage_dict.items():
        norm_usage_dict[codon] = float((D(freq) - D(min_freq)) / (D(max_freq) - D(min_freq)))        
    length_seq = len(sequence)
    num_triplet = length_seq/3
    i = 0   
    score = 0
    while i < length_seq - 2:
        triplet = sequence[i:i+3:1]
        score_tri = norm_usage_dict[triplet]
        score += score_tri
        i += 3
    score_bias = score/num_triplet # scale by the number of triplet in the sequence
    return score_bias


def pos_to_frags(pos,sequence):
    frag_array = []
    if pos:       
        frag_array.append(sequence[0:pos[0]])
        i = 0
        while i < len(pos)-1:
            frag_array.append(sequence[pos[i]+3:pos[i+1]])
            i += 1
        last_seq = sequence[pos[-1]+3:]
        if not any(codon in last_seq for codon in ['TAA','TAG','TGA']):
            frag_array.append(sequence[pos[-1]+3:])
    return frag_array

def exon_extract(temp,pos,EnsID):
    Exons = list(temp.values())[0][pos].split('-')[0] + '|' + list(temp.values())[0][pos].split('-')[1]
    return Exons

def core_match(df_exonlist,dict_exonCoords,EnsID,Exons):
   
    try:
        df_certain = df_exonlist[df_exonlist['EnsGID'] == EnsID]
    except: full_transcript_store = []  # EnsGID is absent in df_exonlist
    full_transcript_store = []
    for item in list(df_certain['Exons']):
        full_transcript=''
        if Exons in item:
            Exonlist = item.split('|')
            for j in range(len(Exonlist)):
                coords = dict_exonCoords[EnsID][Exonlist[j]]
                strand = coords[1]
                judge = check_exonlist_general(Exonlist,j,strand)
                if strand == '+' and judge:   
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                elif strand == '+' and not judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                elif strand == '-' and judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                elif strand == '-' and not judge:
                    frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                    # expression of minus strand, need to draw an illustrator to visulize that.
                full_transcript += frag
            full_transcript = full_transcript.replace('\n','')
            full_transcript_store.append(full_transcript)   
        else: 
            full_transcript_store.append('')
    return full_transcript_store  # ['','ATTTTT','TTTGGCC'], # [] if EnsGID is not present in exonlist

def check_exonlist_general(exonlist,index,strand):
    dict = {}
    for subexon in exonlist:
        exon_num = subexon.split('.')[0]
        subexon_num = subexon.split('.')[1]
        if exon_num in dict:
            dict[exon_num].append(subexon_num)
        else:
            dict[exon_num] = []
            dict[exon_num].append(subexon_num)  # E14 > 1,2,4,5
    # check
    query = exonlist[index]
    query_exon_num = query.split('.')[0]   #E14.1
    query_subexon_num = int(query.split('.')[1])   #2, it is a int
    if strand == '+':
        if str(query_subexon_num + 1) in dict[query_exon_num]: return False
        else: return True
    else:
        if str(query_subexon_num + 1) in dict[query_exon_num]: return False
        else: return True
        
def neoJunctions(df):
    dfNeoJunction = df[((df['dPSI'] >= 0.3) & (df['avg-Others'] <= 0.23))]
    return dfNeoJunction  
    
def uid(df, i):
    uid = list(df['UID'])[i]       
    gene = uid.split(':')[0]
    dict = {}
    gene = gene + ':' + uid.split('|')[1].split(':')[0] # slicing the ENSG in background event
    x = uid.split('|')
    try: x[0].split(':')[3]
    except IndexError: event = x[0].split(':')[2]
    else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])
    finally: dict[gene] = [event]
    
    try: x[1].split(':')[2]
    except IndexError: backEvent = x[1].split(':')[1]
    else: backEvent = str(x[1].split(':')[1])+':'+str(x[1].split(':')[2])
    finally: dict[gene].append(backEvent)

    #{'gene:ENSid':[E22-E33,E34-E56]}
    # if fusion gene: E22-ENSG:E31
    return dict

def subexon_tran(subexon,EnsID,dict_exonCoords,dict_fa,flag): # flag means if it is site_1 or site_2
    try:
        attrs = dict_exonCoords[EnsID][subexon]
        exon_seq = query_from_dict_fa(dict_fa,attrs[2],attrs[3],EnsID,attrs[1])  
    except KeyError:
        if ':' in subexon:   #fusion gene
            fusionGeneEnsID = subexon.split(':')[0] # this kind of subexon must be site2
            fusionGeneExon = subexon.split(':')[1]
            if  '_' in fusionGeneExon:   # ENSG:E2.1_473843893894
                suffix = fusionGeneExon.split('_')[1]
                subexon = fusionGeneExon.split('_')[0]
                attrs = dict_exonCoords[fusionGeneEnsID][subexon]
                exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],fusionGeneEnsID,attrs[1])
            else:    # ENSG:E2.1
                attrs = dict_exonCoords[fusionGeneEnsID][fusionGeneExon]
                exon_seq = query_from_dict_fa(dict_fa,attrs[2],attrs[3],fusionGeneEnsID,attrs[1])
        else:
            try:   #E2.1_67878789798
                suffix = subexon.split('_')[1]
            except IndexError:
                exon_seq = '***********************'   
                print('{0} does not include in {1} exonlists'.format(subexon,EnsID))
            else:
                subexon = subexon.split('_')[0]
                try:
                    attrs = dict_exonCoords[EnsID][subexon]
                except KeyError:
                    chrUTR,strandUTR = utrAttrs(EnsID,dict_exonCoords)
                    exon_seq = utrJunction(suffix,EnsID,strandUTR,chrUTR,flag)  
                    print('{0} observes an UTR event {1}'.format(EnsID,subexon))
                else:
                    if flag == 'site2':           
                        exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],EnsID,attrs[1])  # chr,strand, start,end
                    elif flag == 'site1':
                        exon_seq = query_from_dict_fa(dict_fa,attrs[2],suffix,EnsID,attrs[1])
    return exon_seq

def utrJunction(site,EnsGID,strand,chr_,flag):  # U0.1_438493849, here 438493849 means the site
    if flag == 'site1' and strand == '+':  # U0.1_438493849 - E2.2
        otherSite = int(site) - 100 + 1   # extract UTR with length = 100
        exon_seq = query_from_dict_fa(dict_fa,otherSite,site,EnsGID,strand)
    elif flag == 'site1' and strand == '-':
        otherSite = int(site) + 100 - 1 
        exon_seq = query_from_dict_fa(dict_fa,site,otherSite,EnsGID,strand)
        if exon_seq == '': 
            exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
    elif flag == 'site2' and strand == '+':  # E5.3 - U5.4_48374838
        otherSite = int(site) + 100 -1
        exon_seq = query_from_dict_fa(dict_fa,site,otherSite,EnsGID,strand)
    elif flag == 'site2' and strand == '-':
        otherSite = int(site) - 100 + 1
        exon_seq = query_from_dict_fa(dict_fa,otherSite,site,EnsGID,strand)
    return exon_seq

def utrAttrs(EnsID,dict_exonCoords):
    exonDict = dict_exonCoords[EnsID] 
    attrs = next(iter(exonDict.values()))
    chr_,strand = attrs[0],attrs[1]
    return chr_,strand



def retrieveSeqFromUCSCapi(chr_,start,end):
    import requests
    import xmltodict
    url = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment={0}:{1},{2}'.format(chr_,start,end)
    response = requests.get(url)
    status_code = response.status_code
    assert status_code == 200
    my_dict = xmltodict.parse(response.content)
    exon_seq = my_dict['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n','').upper()
    return exon_seq
    
    
    
    
    
    
def fasta_to_dict(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            temp_list = []
            EnsID = title.split('|')[0]
            chro = title.split('|')[1]
            start = title.split('|')[2]
            end = title.split('|')[3]
            temp_list=[chro,start,end,seq]
            dict_fa[EnsID] = temp_list
    return dict_fa
        
def query_from_dict_fa(dict_fa,abs_start,abs_end,EnsID,strand):
    if strand == '+':        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse,generic_dna).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000 # endpoint in python is non-inclusive
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1,generic_dna)
        exon_seq = str(s.reverse_complement())
    return exon_seq

def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        next(file)
        for line in file:
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5])
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}
    return dict_exonCoords

def backgroundGTEx(N,dict_fa,dictExonList,dict_exonCoords):
    '''
    Run following command in your terminal, get sampleType.txt, wholeBloodID,allsampleID.txt
    
    cut -f 1,7 GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | awk 'BEGIN {FS="\t";} {if (NR>1) {print $1 "\t" $2;}}' > sampleType.txt
    awk '{if ($2 ~/Whole/) {print $1}}' sampleType.txt > wholeBloodID.txt 
    head -n 3 GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct | tail -n 1 | tr '\t' '\n' | awk 'NR>2' > allsampleID.txt
    
    '''
    
    wholeBloodID = list(pd.read_csv('BreastID.txt',sep='\n',header=None,names=['id'])['id'])
    allsampleID = list(pd.read_csv('allsampleID.txt',sep='\n',header=None,names=['id'])['id'])
    overlapID = []
    for id_ in allsampleID:
        if id_ in wholeBloodID: overlapID.append(id_)
    
    if os.path.exists('tranTPMsim.p'):
        with open('tranTPMsim.p','rb') as file1:
            tranTPMsim = pickle.load(file1)
    else:            
        tranTPMexp = pd.read_csv('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sep='\t',skiprows=2,usecols=overlapID)
        tranTPMid = pd.read_csv('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sep='\t',skiprows=2,usecols=['transcript_id','gene_id'])
        tranTPM = pd.concat([tranTPMid,tranTPMexp],axis = 1, sort = False)
        ave = tranTPM.mean(axis = 1, skipna = True, numeric_only = True)
        tranTPM = pd.concat([tranTPM,ave],axis=1,sort=False)
        tranTPMsim = tranTPM.drop(overlapID,axis = 1)
        tranTPMsim = tranTPMsim.rename(columns={0:'mean'})  # change the last average column's name to 'mean'

        with open('tranTPMsim.p','wb') as file2:
            pickle.dump(tranTPMsim,file2)
    '''
    transcript_id  gene_id    mean
      EnsTID        EnsGID    0.08
      
    '''
    # convert to a dict
    if os.path.exists('dictGTEx.p'):
        with open('dictGTEx.p','rb') as file3:
            dictGTEx = pickle.load(file3)
    else:            
        dictGTEx = {}
        for i in range(tranTPMsim.shape[0]):
            EnsGID = tranTPMsim.iat[i,1][:15]
            EnsTID = tranTPMsim.iat[i,0][:15]
            mean = tranTPMsim.iat[i,2]
            if mean > 0.3: highExpression = True
            else: highExpression = False
            if highExpression: merArrayNormal = getNmerNormal(EnsGID,EnsTID,N,dict_fa,dictExonList,dict_exonCoords)  # merArrayNormal is all 9mer for each transcript that expresses in normal tissue
            try: 
                dictGTEx[EnsGID][0][EnsTID] = highExpression
                try: merArrayNormal   
                except NameError: pass  # if no merArrayNormal, means that this transcript is not highly expressed
                else:dictGTEx[EnsGID][1]['wholeMerArrayNormal'].append(merArrayNormal)
            except KeyError: dictGTEx[EnsGID] = [{EnsTID:highExpression},{'wholeMerArrayNormal':[]}]
        with open('dictGTEx.p','wb') as file4:
            pickle.dump(dictGTEx,file4)

        # {'EnsGID:[{EnsTID:False,EnsID:True},{'wholeMerArrayNormal':[['DJFD','FKJDF'],['DJFD','FKJDF']]}]}
    return dictGTEx
 
def convertExonList(df):
    dictExonList = {}
    for i in range(df.shape[0]):
        EnsGID = df.iat[i,0]
        EnsTID = df.iat[i,1]
        exonList = df.iat[i,3]
        try: dictExonList[EnsGID][EnsTID] = exonList
        except KeyError: dictExonList[EnsGID] = {EnsTID:exonList}
        # {EnsGID:{EnsTID:exonlist,EnsTID:exonlist}}
    return dictExonList
    


def merPassGTExCheck(dictGTEx,uid,merArray):   # pass nested list [[],[],[]]
    EnsGID = uid.split('|')[0].split(':')[1]
    normalMer = dictGTEx[EnsGID][1]['wholeMerArrayNormal']
    arrangedNormal = spread(normalMer)
    arrangedSplicing = spread(merArray)
    uniqueInSplicing = list(arrangedSplicing - arrangedNormal)
    return uniqueInSplicing
    

def spread(list_):
    ret = []
    for i in list_:
        ret.extend(i) if isinstance(i,list) else ret.append(i)
    ret = set(ret)
    return ret



    

    
        
    

    
    
    
def flattenNestedList(list_,mode=0,clean=False): # mode = 0 means returning flatten list, mode = 1 means converting to set, mode = 2 means converting to list(set())
    # input a nested list, only extract element from non-empty sub list.
    flatten = []
    if not clean:
        for sublist in list_:
            flatten = [item for item in sublist]
        if mode == 0: return flatten
        elif mode == 1: return set(flatten)
        elif mode == 2: return list(set(flatten))
    elif clean:   # clean will remove '*',[],'MANNUAL'
        for sublist in list_:
            if isinstance(sublist,list) and sublist:  # if the sub is a non-empty list
                flatten = [item for item in sublist]
        if mode == 0: return flatten
        elif mode == 1: return set(flatten)
        elif mode == 2: return list(set(flatten))
            




           
def getNmerNormal(EnsID,EnsTID,N,dict_fa,dictExonList,dict_exonCoords):
    try:exonlist = dictExonList[EnsID][EnsTID] 
    except KeyError: print('{0} of {1} doesn\'t exist in mRNA-ExonID file'.format(EnsTID,EnsID))                             
    # obtain the full-length sequence for this exonlist to check its resultant 9mer
    else:
        Exonlist = exonlist.split('|')  # [E1.1,E1.2,E3.4....]
        full_transcript = ''
        for j in range(len(Exonlist)):
            coords = dict_exonCoords[EnsID][Exonlist[j]]
            strand = coords[1]
            judge = check_exonlist_general(Exonlist,j,strand)
            if strand == '+' and judge:   
                frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
            elif strand == '+' and not judge:
                frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
            elif strand == '-' and judge:
                frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
            elif strand == '-' and not judge:
                frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
            # expression of minus strand, need to draw an illustrator to visulize that.
            full_transcript += frag
        full_transcript = full_transcript.replace('\n','')
    
    
        # continue to get ORF
        ORF = transcript2peptide(full_transcript)
        peptide = str(Seq(ORF,generic_dna).translate(to_stop=False))
        merArray = extractNmer(peptide,N)  # all the 9mer the normal
        return merArray   # ['SRTTIFDFJD','FJDJFKDJFKDJ','FJKDFJKDF']
    
def intron(event,EnsGID,junction,dict_exonCoords,dictExonList,N):
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
                            if remainder == 0: 
                                front=0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front=0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+1)].replace(',','')
                        elif strand == '-':
                            intronStartIndex = int(attrs[3]) - 1
                            remainder = (tranStartIndex - intronStartIndex) % 3 
                            if remainder == 0: 
                                front=0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3)].replace(',','')
                            elif remainder == 1: 
                                front=0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+2)].replace(',','')
                            elif remainder == 2: 
                                front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                junctionSlice = junction[front:junction.find(',')+((N-1)*3+1)].replace(',','')
                
                        peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                        merArray = extractNmer(peptide,N) 
                        merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
        
    elif event.startswith('I'): merBucket = [[]]
    return merBucket


# https://rest.ensembl.org/documentation/info/lookup
def grabEnsemblTranscriptTable(EnsTID):
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/{0}?expand=1".format(EnsTID)     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})     
    decoded = r.json()
    try: translationStartIndex = decoded['Translation']['start']
    except KeyError: print('{0} is not translatable'.format(EnsTID)) # might be invalid transcriptID or they don't have tranStartIndex(don't encode protein)
    else: return translationStartIndex
# for except branch, we don't specify return command, so it will return NoneType   

def toFasta(list_,N):
    with open('./resultMHC/queryNetMHC_{0}.fa'.format(N),'w') as file1:
        for index,item in enumerate(list_):
            file1.write('>{0} mer\n'.format(index+1))
            file1.write('{0}\n'.format(item.strip('\'')))


if __name__ == "__main__":

    startTime = process_time()
    df = pd.read_csv('PSI.R1-V6.MergedResult.MergedFiles.2_vs_Others.txt',sep='\t') 
    df_exonlist = pd.read_csv('mRNA-ExonIDs.txt',sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    dict_exonCoords = exonCoords_to_dict('Hs_Ensembl_exon.txt','\t')
    dict_fa = fasta_to_dict('Hs_gene-seq-2000_flank.fa')
    metaBaml = Meta(df) #Instantiate Meta object
    dfNeoJunction = neoJunctions(metaBaml.df)
    NeoJBaml = NeoJ(dfNeoJunction,9) #Instantiate NeoJ object
    NeoJBaml.retrieveJunctionSite()
    NeoJBaml.matchWithExonlist(df_exonlist,dict_exonCoords)
    NeoJBaml.getORF()
    NeoJBaml.getORFaa()
    NeoJBaml.phaseTranslateJunction()
    NeoJBaml.getNmer()
    
    dictExonList = convertExonList(df_exonlist)  
    #dictGTEx = backgroundGTEx(NeoJBaml.mer,dict_fa,dictExonList,dict_exonCoords)
    
    #NeoJBaml.getFinalMers()
    NeoJBaml.mannual()
    #NeoJBaml.collectNeoAntigen()
    #toFasta(NeoJBaml.mhcNeoAntigens,NeoJBaml.mer)
    
    
    NeoJBaml.df.to_csv('./resultMHC/NeoJunction_{0}.txt'.format(NeoJBaml.mer),sep='\t',header=True,index = False)
    
    endTime = process_time()
    print('Time Usage: {} seconds'.format(endTime-startTime))
    #print(len(NeoJBaml.mhcNeoAntigens))
    




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




