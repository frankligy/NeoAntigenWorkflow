#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:32:27 2020

@author: ligk2e
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bz2
import _pickle as cpickle
import pickle
import re
from time import process_time
from urllib.error import HTTPError
import requests
import xmltodict
import multiprocessing as mp
import subprocess
from pathos.multiprocessing import ProcessingPool as Pool
import getopt
from decimal import Decimal as D




class Meta():  #inspect an object: dir(), vars(), instanceName.__dict__, mannually set __repr__ or __str__
    def __init__(self, df):
        
        self.df = df
    
    def __repr__(self):
        return {'df':self.df}
    
    def getPercent(self,dfgroup,dfori,name,write=False):
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
            
            allvalue_t = list(map(lambda x:round(x,2),dfori.iloc[j][dic['R1-V7']].tolist()))
            truth_t = [True if item == 0 else False for item in allvalue_t]
            allzero_t = sum(truth_t)   # how many zeros in this cluster
            percent_t = allzero_t/num_t  # percentage of zeros
            percentArray_t.append(percent_t)
                
            #allvalue = dfori[dic['R1-V7']].iloc[j].tolist()
            allvalue_h = list(map(lambda x:round(x,2),dfori.iloc[j][dic['Healthy']].tolist()))
            truth_h = [True if item == 0 else False for item in allvalue_h]
            allzero_h = sum(truth_h)   # how many zeros in this cluster
            percent_h = allzero_h/num_h   # percentage of zeros
            percentArray_h.append(percent_h)
            dicPercent[event]=(percent_t,allvalue_t,percent_h,allvalue_h)
            
        col0,col1,col2,col3 = [],[],[],[]
        for k in range(self.df.shape[0]):
            splice = self.df['UID'].tolist()[k]
            per_t,all_t,per_h,all_h = dicPercent[splice][0],dicPercent[splice][1],dicPercent[splice][2],dicPercent[splice][3]
            col0.append(per_t)
            col1.append(all_t)
            col2.append(per_h)
            col3.append(all_h)
        self.df['tumor_zero_percent'] = col0
        self.df['tumor_distribution'] = col1
        self.df['healthy_zero_percent'] = col2
        self.df['healthy_distribution'] = col3
        if write==True: self.df.to_csv('see{0}.txt'.format(name),sep='\t',index=None)
            
    def retrieveJunctionSite(self,dict_exonCoords,dict_fa):
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
                        print('The {}th transcript of {}, even though junction site could match with, but optimal ORF suggests that junction site will not be translated'.format(j,list(self.df['UID'])[i]))
                    else: 
                        former = junctionOri.find(',')  # length of former part, also the index of splicesite
                        phase = abs(startJun + former - startORF) % 3 
                        N = self.mer
                         
                        if phase == 0: 
                            front=0 if former - ((N-1)*3) < 0 else (former - (N-1)*3)
                            
                            junctionSlice = junction[front: former + ((N-1)*3)].replace(',','')
                            #print(junctionSlice)
                        elif phase == 1: 
                            front=0 if former - ((N-1)*3+1) <0 else (former - ((N-1)*3+1))
                            junctionSlice = junction[front: former + ((N-1)*3+2)].replace(',','')
                            
                        elif phase == 2: 
                            front=0 if former - ((N-1)*3+2) <0 else (former - ((N-1)*3+2))
                            junctionSlice = junction[front: former + ((N-1)*3+1)].replace(',','')
                 
            
                        peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))

                        phaseArray.append(phase)
                        peptideArray.append(peptide)
            
            col0.append(phaseArray) 
            col1.append(peptideArray)
        self.df['phase'] = col0 # ['',1,'*',2] OR ['#']
        self.df['peptide'] = col1  # ['',MSTTG,'*',TYCTT] OR ['#']
                
    def getNmer(self):
        col = []
        for i in range(self.df.shape[0]):
            print(i,'first round-Process{}'.format(os.getpid()))
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
            truthTable1 = [False if mer == '*' or mer == [] or mer== '' else True for mer in merArray]
            if not any(truthTable1): merArray = ['Match but no translation or early stop']
            col.append(merArray)  # [ '','*',[], ['MTDJDKFDJF','FJDKFJDF'] ]
        self.df['{0}mer'.format(self.mer)] = col
        
                 
    def mannual(self):
        col = []
        for i in range(self.df.shape[0]):
            print(i,'second mannual round-process{}'.format(os.getpid()))
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
                if 'ENSG' in event:  # E4.3--ENSG00000267881:E2.1
                    # nearly the same as newSpliicng, query is former subexon, trailing is replaced as breaking point(start coordinate of former subexon + length of former part - 1)
                    merArray = tranSplicing(event,EnsGID,junction,self.mer,dictExonList,dict_exonCoords)
                    if merArray == [[]]: merArray = ['transplicing, already checked, query subexon not present in known transcript']
                if re.search(r'I.+_\d+',event) or 'U' in event:  # they belong to NewExon type: including UTR type and blank type(the annotation is blank)
                # for NewExon type, no need to infer translation phase, just use junction sequence, mer should span junction site
                    junctionIndex = junction.find(',')
                    Nminus1 = self.mer - 1 
                    front = 0 if junctionIndex-(Nminus1*3+2) < 0 else junctionIndex-(Nminus1*3+2)
                    junctionSlice = junction[front:junctionIndex+(Nminus1*3+2)].replace(',','') # no need to worry out of index, slicing operation will automatically change it to 'end' if overflow
                    merArray = dna2aa2mer(junctionSlice,self.mer)   # merArray is a nested list
                if re.search(r'E\d+\.\d+_\d+',event):  # they belong to newSplicing site or fusion gene
                # for this type, just check if the former part (E1.2-E3.4, here E1.2 is former part) exist in known transcript, then infer the phase of translation
                    #print(event)
                    merArray = newSplicingSite(event,EnsGID,junction,self.mer,dictExonList,dict_exonCoords)   # nested list
                    if merArray == [[]]: merArray = ['newSplicing, already checked, query subexon not present in known transcript']

                if re.search(r'^I\d+\.\d+-',event) or re.search(r'-I\d+\.\d+$',event):
                    merArray = intron(event,EnsGID,junction,dict_exonCoords,dictExonList,self.mer)  # nested list
                    if merArray == [[]]: merArray = ['intron retention, already checked, either former subexon not present in known transcript or matched transcript is not Ensembl-compatible ID']
                if re.search(r'E\d+\.\d+-E\d+\.\d+$',event):   # novel splicing: CLASP1:ENSG00000074054:E25.1-E26.1, just don't match with any existing one
                    # here we extract the backEvent, which is E24.1-E27.1
                    #print(event)
                    backEvent = backEvent.replace('-','|')  # now it is E24.1|E27.1
                    merArray = novelOrdinal(event,backEvent,EnsGID,junction,dictExonList,dict_exonCoords,self.mer)
                    if merArray ==[[]]: merArray = ['novel ordinal splicing event, already checked, its background event can not match up with any existing transcript either']
                    
            col.append(merArray)
        self.df['mannual'] = col
        

    
    def netMHCresult(self,HLA,pathSoftWare,mode,sb=0.5,wb=2.0):
        col = []
        for i in range(self.df.shape[0]):
            #print(i)
            merList = self.df['{0}mer'.format(self.mer)].tolist()[i]
            if merList == ['MANNUAL']: merList = self.df['mannual'].tolist()[i]
            #print(merList)
            netMHCpan.pepFile(merList) # get query.pep file
            machine = netMHCpan('./resultMHC/temp/query{}.pep'.format(os.getpid()),HLA,pathSoftWare,self.mer,mode,sb,wb)
            dic = machine.seperator()

            col.append(dic)
        self.df['{0}result'.format(mode)]=col
            

class netMHCpan():
    # ../netMHCpan -p test.pep -BA -xls -a HLA-A01:01,HLA-A02:01 -xlsfile my_NetMHCpan_out.xls
    def __init__(self,intFile,HLA,pathSoftWare,length,mode,sb=0.5,wb=2.0):
        self.intFile=intFile
        self.HLA=HLA
        self.pathSoftWare=pathSoftWare
        self.sb = sb    # float, 0.5
        self.wb = wb    # float, 2.0
        self.length = length
        self.filter = wb  # float, 2.0
        self.mode=mode
    
    @staticmethod
    def pepFile(lis):
        result = []
        [result.extend(item) for item in lis if isinstance(item,list)]
        result = list(set(result))
        with open('./resultMHC/temp/query{}.pep'.format(os.getpid()),'w') as f1:
            [f1.write('{}\n'.format(mer)) for mer in result]
        # now we have query.pep file
    
    def seperator(self):
        if self.mode=='MHCI': 
            self.runSoftWareI()
            dic = self.postFileI()
            return dic
        elif self.mode=='MHCII':
            self.runSoftWareII()
            dic = self.postFileII()
            return dic 
            

    
    def runSoftWareI(self):
        with open('./resultMHC/temp/resultI{}.txt'.format(os.getpid()),'w') as f3:
            subprocess.run([self.pathSoftWare,'-p',self.intFile, '-BA','-a',self.HLA,'-rth', str(self.sb), '-rlt', str(self.wb), '-l',str(self.length),'-t',str(self.wb)],stdout=f3)        
       # will generate a result file  
    
    def runSoftWareII(self):
        with open('./resultMHC/temp/resultII{}.txt'.format(os.getpid()),'w') as f4:
            subprocess.run([self.pathSoftWare,'-f',self.intFile,'-inptype', '1', '-a',self.HLA,'-length',str(self.length)],stdout=f4)
    
    def postFileII(self):
        with open('./resultMHC/temp/resultII{}.txt'.format(os.getpid()),'r') as f5,open('./resultMHC/temp/resultII_parse{}.txt'.format(os.getpid()),'w') as f6:
            for line in f5:
                if line.startswith('#') or line.startswith('-') or line.strip('\n') == '':continue
                elif re.search(r'^\w+',line): continue
                elif re.search(r'Pos',line): continue
                else: f6.write(line)
        try:df=pd.read_csv('./resultMHC/temp/resultII_parse{}.txt'.format(os.getpid()),sep='\s+',header=None,index_col=0,names=[str(i+1) for i in range(11)])
        except pd.errors.EmptyDataError: dic = 'No candidates'
        else:
            hlaAllele = df['2'].tolist()  # HLA Allele
            mer = df['3'].tolist()   # kmer amino acid
            level = df['11'].tolist()  # <=SB or <=WB
            hla = self.HLA
    
            hlaList = hla.split(',')
            hlaNum = len(hlaList) 
            dic = {}   # store all the candidates with binding affinity
            for i in range(hlaNum):
                sb,wb=[],[]   # all strong binding neoantigen, all weak binding neoantigen
                hlaQuery = hlaList[i]
                occurence = [k for k in range(len(hlaAllele)) if hlaAllele[k] == hlaQuery]
                for j in occurence:
                    if level[j]=='<=SB': sb.append(mer[j])
                    elif level[j]=='<=WB':wb.append(mer[j])
                dic[hlaList[i]] = (sb,wb)
            self.neoantigen = dic
        return dic   
    
    
    def postFileI(self):
        # HLA = 'HLA-A01:01,HLA-A03:01,HLA-B07:02,HLA-B27:05,HLA-B58:01'
        with open('./resultMHC/temp/resultI{}.txt'.format(os.getpid()),'r') as f1, open('./resultMHC/temp/resultI_parse{}.txt'.format(os.getpid()),'w') as f2:
            for line in f1:
                if line.startswith('#') or line.startswith('-') or line.strip('\n') == '': continue
                elif re.search(r'^\w+',line): continue
                elif re.search(r'Pos',line): continue
                else: f2.write(line)
        try:df = pd.read_csv('./resultMHC/temp/resultI_parse{}.txt'.format(os.getpid()),sep='\s+', header=None,index_col=0)  
        except pd.errors.EmptyDataError: dic = 'No candidates'   
        else:
            hlaAllele = df[1].tolist()  # HLA Allele
            mer = df[2].tolist()   # kmer amino acid
            level = df[17].tolist()  # SB or WB
            hla = self.HLA
            
            hlaList = hla.split(',')
            hlaNum = len(hlaList) 
            dic = {}   # store all the candidates with binding affinity
            for i in range(hlaNum):
                sb,wb=[],[]   # all strong binding neoantigen, all weak binding neoantigen
                hlaQuery = netMHCpan.hlaType(hlaList[i])  # HLA-A01:01 to HLA-A*01:01
                occurence = [k for k in range(len(hlaAllele)) if hlaAllele[k] == hlaQuery]
                [sb.append(mer[j]) if level[j]=='SB' else wb.append(mer[j]) for j in occurence]
                dic[hlaList[i]] = (sb,wb)
            self.neoantigen = dic
        return dic


            
            
            
            
            
    @staticmethod
    def hlaType(hla):   # convert HLA-A01:01 to HLA-A*01:01
        index1 = re.search(r'^HLA-[A-Z]+',hla).span()[1]   # here, index will be 5, which is '0'
        former,latter = hla[0:index1],hla[index1:]
        hlaNew = former + '*' + latter
        return hlaNew
            

        



    

        
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
                            
def tranSplicing(event,EnsGID,junction,N,dictExonList,dict_exonCoords):  # E4.3-ENSG00000267881:E2.1
    query = event.split('-')[0]
    formerLength = len(junction.split(',')[0])
    merBucket = []
    try: attrs = dict_exonCoords[EnsGID][query]
    except KeyError: #E6.1_3854692-ENSG00000267881:E2.1
        print('complicated situation(both trailing and transplicing):{0},{1}'.format(EnsGID,event))
        trailing = query.split('_')[1]  #  3854692
        breaking = int(trailing)+1
    else:
        strand = attrs[1]
        start = attrs[2]
        breaking = int(start)+formerLength
    finally:
        allTranDict = dictExonList[EnsGID]
        for tran,exonlist in allTranDict.items():
            if query in exonlist: 
                try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                else:
                    if type(tranStartIndex) == int:
                        if strand == '+':
                            remainder = (int(breaking) - tranStartIndex) % 3   
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
    
                            remainder = (tranStartIndex - int(breaking)) % 3 
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
    #print(event)
    try: event.split('-')[0].split('_')[1]   # see if the former one has trailing part E1.2_8483494
    except IndexError: 
        former = event.split('-')[0]  # no trailing part, former part just E1.2, means newSplicingSite is in latter part
        query = former
        attrs = dict_exonCoords[EnsGID][query]
        start = attrs[2]
        trailing = int(start) + len(junction.split(',')[0])   # the coordinates of latter part position1
    else: 
        query = event.split('-')[0].split('_')[0]  # has trailing part, former part should get rid of trailing part
        trailing = int(event.split('-')[0].split('_')[1]) + 1    # the coordinates of latter part position1
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
        
def neoJunctions_oldMethod_noGTEx(df,colname):
    dfNeoJunction = df[((df['dPSI'] >= 0.3) & (df[colname] <= 0.23))]
    return dfNeoJunction

def neoJunction_oldMethod_checkGTEx(df,colname):
    condition = []
    for i in range(df.shape[0]):
        event = df.iloc[i]['UID']
        healthy = df.iloc[i][colname]
        dPSI = df.iloc[i]['dPSI']
        if dPSI >= 0.3 and healthy <= 0.23: 
            try:inspect = inspectGTEx(event,plot=False)
            except KeyError: cond = True    # absent in normal GTEx data set
            else: cond = True if inspect==49 else False  # =49 means have no expression in all tissue
        else: cond = False
        condition.append(cond)
    condition = pd.Series(condition)
    df_Neo = df[condition]
    return df_Neo

def neoJunction_newMethod_checkGTEx(df):
    condition = []
    for i in range(df.shape[0]):
        event = df.iloc[i]['UID']
        per_h = df.iloc[i]['healthy_zero_percent']
        if per_h > 0.95:
            try:inspect = inspectGTEx(event,plot=False)
            except KeyError: cond = True    # absent in normal GTEx data set
            else: cond = True if inspect==49 else False  # =49 means have no expression in all tissue
        else: cond = False
        condition.append(cond)
    condition = pd.Series(condition)
    df_Neo = df[condition]
    return df_Neo


def neoJunction_newMethod_noGTEx(df):
    condition = []
    for i in range(df.shape[0]):
        event = df.iloc[i]['UID']
        per_h = df.iloc[i]['healthy_zero_percent']
        if per_h > 0.95:
            cond = True
        else: cond = False
        condition.append(cond)
    condition = pd.Series(condition)
    df_Neo = df[condition]
    return df_Neo
        
def neoJunction_testMethod_noGTEx(df):
    return df        
        
def neoJunction_testMethod_checkGTEx(df):
    condition = []
    for i in range(df.shape[0]):
        event = df.iloc[i]['UID']

        try:inspect = inspectGTEx(event,plot=False)
        except KeyError: cond = True    # absent in normal GTEx data set
        else: cond = True if inspect==49 else False  # =49 means have no expression in all tissue
        condition.append(cond)
    condition = pd.Series(condition)
    df_Neo = df[condition]
    return df_Neo          
        
    
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
    


    

def spread(list_):
    ret = []
    for i in list_:
        ret.extend(i) if isinstance(i,list) else ret.append(i)
    ret = set(ret)
    return ret

    
def intron(event,EnsGID,junction,dict_exonCoords,dictExonList,N):
    merBucket = []
    if event.startswith('E'):   # only consider this situation, since intron retentions are predominantly associated with a translatable preceding exon, ending up with a early stop codon
        # namely: E2.4-I2.1, E22.1-I22.1
        former = event.split('-')[0]  # E2.4, E22.1
        latter = event.split('-')[1]
        #print(event,EnsGID,former)
        try: attrs_former = dict_exonCoords[EnsGID][former] # chr, strand, start, end
        except KeyError: print('preexisting bug, Event {0} in {1} doesn\'t have {2}!!!!!'.format(event,EnsGID,former))
        else:
            attrs_latter = dict_exonCoords[EnsGID][latter]
            strand = attrs_latter[1]        
            allTransDict = dictExonList[EnsGID] 
            for tran,exonlist in allTransDict.items():
                if former in exonlist: 
                    #print(tran,'*******************************************************')
                    try:tranStartIndex = grabEnsemblTranscriptTable(tran)
                    except HTTPError: continue   # for instance, BC4389439 it is not a Ensemble-valid transcript ID, just pass this transcript
                    else:
                        if type(tranStartIndex) == int:
                            if strand == '+':
                                intronStartIndex = int(attrs_latter[3])
                                remainder = (intronStartIndex - tranStartIndex) % 3   # how many nts remaining before the first nt in intron
                                # if 0, means the first nt in intron will be the first nt in codon triplet,
                                # if 1, means the first nt in intron will be the second nt in codon triplet,
                                # if 2, means the first nt in intron will be the third nt in codon triplet.
                                if remainder == 0: 
                                    front=0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                    junctionSlice = junction[front:].replace(',','')
                                elif remainder == 1: 
                                    front=0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                    junctionSlice = junction[front:].replace(',','')
                                elif remainder == 2: 
                                    front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                    junctionSlice = junction[front:].replace(',','')
                            elif strand == '-':
                                intronStartIndex = int(attrs_latter[3]) - 1
                                remainder = (tranStartIndex - intronStartIndex) % 3 
                                if remainder == 0: 
                                    front=0 if junction.find(',') - ((N-1)*3)<0 else junction.find(',') - ((N-1)*3)
                                    junctionSlice = junction[front:].replace(',','')
                                elif remainder == 1: 
                                    front=0 if junction.find(',') - ((N-1)*3+1)<0 else junction.find(',') - ((N-1)*3+1)
                                    junctionSlice = junction[front:].replace(',','')
                                elif remainder == 2: 
                                    front=0 if junction.find(',') - ((N-1)*3+2)<0 else junction.find(',') - ((N-1)*3+2)
                                    junctionSlice = junction[front:].replace(',','')
                    
                            peptide = str(Seq(junctionSlice,generic_dna).translate(to_stop=False))
                            merArray = extractNmer(peptide,N) 
                            merBucket.append(merArray)
        if merBucket == []: merBucket = [[]]  # means no match for the former subexon.
        
    elif event.startswith('I'): merBucket = [[]]
    return merBucket


# https://rest.ensembl.org/documentation/info/lookup
def grabEnsemblTranscriptTable(EnsTID):
    global counter
    print(EnsTID,'**********************************************************')
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/{0}?expand=1".format(EnsTID)     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})     
    try: decoded = r.json()
    except: 
        counter += 1
        print('JSON unknoen error {}'.format(counter))
    else:
        try: translationStartIndex = decoded['Translation']['start']
        except KeyError: print('{0} is not translatable'.format(EnsTID)) # might be invalid transcriptID or they don't have tranStartIndex(don't encode protein)
        else: return translationStartIndex
# for except branch, we don't specify return command, so it will return NoneType   

def toFasta(list_,N):
    with open('./resultMHC/queryNetMHC_{0}.fa'.format(N),'w') as file1:
        for index,item in enumerate(list_):
            file1.write('>{0} mer\n'.format(index+1))
            file1.write('{0}\n'.format(item.strip('\'')))

def inspectGTEx(event,tissue='all',plot=True):
    flag = 0

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


def run(NeoJBaml):
    PID = os.getpid()

    print('retrieve all junction site sequence-Process:{0}\n'.format(PID))
    NeoJBaml.retrieveJunctionSite(dict_exonCoords,dict_fa)
    print('finished retrieving all junction site sequence-Process:{0}\n'.format(PID))
    print('inspecting each event see if they could match up with any existing transcripts-Process:{0}\n'.format(PID))
    NeoJBaml.matchWithExonlist(df_exonlist,dict_exonCoords)
    print('finished inspecting each event see if they could match up with any existing transcripts-Process:{0}\n'.format(PID))
    print('getting most likely ORF and ORFaa for each event that could match up with existing transcript-Process:{0}\n'.format(PID))
    NeoJBaml.getORF()
    NeoJBaml.getORFaa()
    print('finished getting most likely ORF and ORFaa-Process:{0}\n'.format(PID))
    print('checking the phase of each event\'s junction site-Process:{0}\n'.format(PID))
    NeoJBaml.phaseTranslateJunction()
    print('finished checking the phase of each event\'s junction site-Process:{0}\n'.format(PID))
    print('getting Nmer, only for first round-Process:{0}\n'.format(PID))
    NeoJBaml.getNmer()
    print('first round ends-Process:{0}\n'.format(PID))

    print('starting to mannal check,second round-Process:{0}\n'.format(PID))
    NeoJBaml.mannual()
    print('finished mannual check-Process:{0}\n'.format(PID))
    print('starting to deploy netMHCpan-Process:{0}\n'.format(PID))
    if MHCmode == 'MHCI':
        NeoJBaml.netMHCresult(HLA,software,MHCmode)
    elif MHCmode == 'MHCII':
        NeoJBaml.netMHCresult(HLA,software,MHCmode) 
    print('finished binding affinity prediction-Process:{0}\n'.format(PID))
    return NeoJBaml.df


def main(intFile,taskName,k,HLA,software,MHCmode,mode,Core=mp.cpu_count(),checkGTEx=False,EventAnnotationFile='dummy.txt',GroupsFile='proxy.txt'):

    startTime = process_time()
    global df
    global df_exonlist
    global dict_exonCoords
    global dict_fa
    global dictExonList
    print('loading input files\n')
    df = pd.read_csv(intFile,sep='\t') 
    print('finished loading input files\n')
    print('loading all existing transcripts files\n')
    df_exonlist = pd.read_csv('./data/mRNA-ExonIDs.txt',sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    print('finished loading all existing transcripts files\n')
    print('loading all subexon coordinates files\n')
    dict_exonCoords = exonCoords_to_dict('./data/Hs_Ensembl_exon.txt','\t')
    print('finished loading all subexon coordinates files\n')
    print('loading exon sequence fasta files, 2GB\n')
    dict_fa = fasta_to_dict('./data/Hs_gene-seq-2000_flank.fa')
    print('finished loading exon sequence fasta files\n')
    print('converting subexon coordinates to a dictionary-Process\n')
    dictExonList = convertExonList(df_exonlist)
    print('finished converting subexon coordinates to a dictionary-Process\n')
    if mode == 'OncoSplice' and checkGTEx==True:
        dfori = pd.read_csv(EventAnnotationFile,sep='\t')
        dfgroup = pd.read_csv(GroupsFile,sep='\t',header=None,names=['TCGA-ID','group','label'])
        print('loading GTEx dataset, it will take 20 mins, please be patient')

        start = process_time()

        with bz2.BZ2File('dicTissueExp2.pbz2','rb') as f1:
            dicTissueExp = cpickle.load(f1)  
            end = process_time()    
        print('consume {0}'.format(end-start))
        
        metaBaml = Meta(df) #Instantiate Meta object
        metaBaml.getPercent(dfgroup,dfori,'{0}_All'.format(taskName),write=True)
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_newMethod_checkGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')
    
    if mode == 'OncoSplice' and checkGTEx == False:
        dfori = pd.read_csv(EventAnnotationFile,sep='\t')
        dfgroup = pd.read_csv(GroupsFile,sep='\t',header=None,names=['TCGA-ID','group','label'])
        metaBaml = Meta(df) #Instantiate Meta object
        metaBaml.getPercent(dfgroup,dfori,'{0}_All'.format(taskName),write=True)
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_newMethod_noGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')
    
    
        
    if mode == 'TumorAndControl' and checkGTEx == True:
        print('loading GTEx dataset, it will take 20 mins, please be patient')
        start = process_time()
        with bz2.BZ2File('./data/dicTissueExp2.pbz2','rb') as f1:
            dicTissueExp = cpickle.load(f1)  
            end = process_time()    
        print('consume {0}'.format(end-start))
        metaBaml = Meta(df) #Instantiate Meta object
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_oldMethod_checkGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')
        
    if mode=='TumorAndControl' and checkGTEx == False:
        metaBaml = Meta(df) #Instantiate Meta object
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_oldMethod_noGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')
        
    if mode == 'singleSample' and checkGTEx == True:
        print('loading GTEx dataset, it will take 20 mins, please be patient')
        start = process_time()
        with bz2.BZ2File('dicTissueExp2.pbz2','rb') as f1:
            dicTissueExp = cpickle.load(f1)  
            end = process_time()    
        print('consume {0}'.format(end-start))
        metaBaml = Meta(df) #Instantiate Meta object
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_testMethod_checkGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')
        
    if mode == 'singleSample' and checkGTEx == False:
        metaBaml = Meta(df) #Instantiate Meta object
        print('generate NeoJunctions\n')
        dfNeoJunction = neoJunction_testMethod_noGTEx(metaBaml.df)
        print('NeoJunctions are ready\n')
       
    print('start analysis and spawn subprocesses\n')




    df_split = np.array_split(dfNeoJunction, Core, axis=0)    # cut by rows, equally splited 
    obj_split = [NeoJ(df,k) for df in df_split]
    pool = Pool(Core)
    df_out_list = pool.map(run, obj_split)
    pool.close()
    pool.join()
    Crystal = pd.concat(df_out_list)   # default is stacking by rows
    
    Crystal.to_csv('./resultMHC/NeoJunction_{0}_mark.txt'.format(k),sep='\t',header=True,index = False)
    
    endTime = process_time()
    print('Time Usage: {} seconds'.format(endTime-startTime))


def usage():


    print('Usage:')
    print('python3 mhcPresent.py -i ./Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt -t test -k 8 -H HLA-A29:02,HLA-B51:01,HLA-B54:01,HLA-B57:01 -s /data/salomonis2/LabFiles/Frank-Li/python3/netMHCpan-4.1/netMHCpan -M MHCI -m singleSample -C 8 -c False')
    print('Options:')
    print('-i : path of input file')
    print('-t : the name of your task')
    print('-k : Kmer')
    print('-H : queried HLA alleles')
    print('-s : full path for where netMHCpan sit')
    print('-M : MHCmode, either MHCI or MHCII')
    print('-m : mode for generating Neo-Junctions')
    print('-C : how many processes you wanna spawn')
    print('-c: check GTEx data or not')
    print('--EventAnnotationFile: Oncosplice mode, path of EventAnnotation File')
    print('--GroupsFile: Oncosplice mode, path of OncoSplice File')
    print('-h --help : check help information ')
    print('Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati, 2020') 
    



if __name__ == "__main__":

#    log_err = open('queryGTEx.stderr.log','a')
#    log_out = open('queryGTEx.stdout.log','a')
#    sys.stderr = log_err
#    sys.stdout = log_out
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'hi:t:k:H:s:M:m:C:c:',['help','EventAnnotationFile=','GroupsFile='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('-i'):
            intFile = arg
            print('Input file is:', arg)
        elif opt in ('-t'):
            taskName = arg
            print('give your task a name:',arg)
        elif opt in ('-k'):
            k= int(arg)
            print('kmer:', arg)
        elif opt in ('-H'):
            HLA = arg
            print('Queried HLA allele:',arg)
        elif opt in ('-s'):
            software = arg
            print('Full path of netMHCpan standalone version:',arg)
        elif opt in ('-M'):
            MHCmode = arg
            print('MHCI or MHCII:',arg)
        elif opt in ('-m'):
            mode = arg
            print('How to get Neo-Junction:',arg)
        elif opt in ('-C'):
            Core = int(arg)
            print('How many processes to use:',arg)
        elif opt in ('-c'):
            checkGTEx = bool(arg)
            print('check GTEx or not:',arg)
        elif opt in ('--EventAnnotationFile'):
            EventAnnotationFile = arg
            print('Oncosplice mode, you have to specify full path of EventAnnotation:',arg)
        elif opt in ('--GroupsFile'):
            GroupsFile = arg
            print('Oncosplice mode, you have to specify full path of GroupsFile:',arg)
        elif opt in ('--help','-h'):
            usage() 
            sys.exit() 

    counter = 0
    main(intFile,taskName,k,HLA,software,MHCmode,mode,Core=mp.cpu_count(),checkGTEx=False,EventAnnotationFile='dummy.txt',GroupsFile='proxy.txt')
    print(counter)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




