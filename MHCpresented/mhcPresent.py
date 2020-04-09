#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:32:27 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop/project/')
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
import backgroundGTEx as GTEx

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
            exam_seq.append(exam_site_1 + exam_site_2)
            back_site_1 = subexon_tran(back_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
            back_site_2 = subexon_tran(back_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
            back_seq.append(back_site_1 + back_site_2)
            
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
    def __init__(self,df):
        self.df=df

    def phaseTranslateJunction(self):
        col0,col1 = [],[]
        for i in range(self.df.shape[0]):
            wholeTrans = list(self.df['exam_whole_transcripts'])[i]
            ORFtrans = list(self.df['exam_ORF_tran'])[i]
            junction = list(self.df['exam_seq'])[i]
            phaseArray,peptideArray = [],[]
            if not wholeTrans: 
                phaseArray.append('#')
                peptideArray.append('#')
                print('This gene {0} is absent in mRNA-ExonID file\n'.format(list(self.df['UID'])[i]))
            for j in range(len(wholeTrans)):
                if not wholeTrans[j]: pass   # splicing evnet can not match to this transcript
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
                        phase = abs(startJun - startORF) % 3 # phase 0 means junction will translate in 0 base
                        peptide = str(Seq(junction[phase:-1],generic_dna).translate(to_stop=False))
                        phaseArray.append(phase)
                        peptideArray.append(peptide)
            
            col0.append(phaseArray) 
            col1.append(peptideArray)
        self.df['phase'] = col0
        self.df['peptide'] = col1
                
    def get9mer(self):
        col = []
        for i in range(self.df.shape[0]):
            peptides = list(self.df['peptide'])[i]
            merArray = []
            if not peptides: merArray.append('MANNUAL')
            else:
                for peptide in peptides:
                    if peptide == '*': merArray.append('*')
                    else: 
                        tempArray = extractNmer(peptide,9)
                        merArray.append(tempArray)
            truthTable1 = [False if mer == '*' or mer == [] else True for mer in merArray]
            if not any(truthTable1): merArray = ['Match but no translation or ealy stop, MANNUAL']
            col.append(merArray)
        self.df['9mer'] = col
        
    def getFinalMers(self):
        col = []
        for i in range(self.df.shape[0]):
            uid, merArray = list(self.df['UID'])[i], list(self.df['9mer'])[i]
            col.append(GTEx.merPassGTExCheck(dictGTEx,uid,merArray))
        self.df['matchedFinalMer'] = col
                 
             
             
def extractNmer(peptide,N):
    starIndex = peptide.find('*')  
    merArray = []
    if starIndex == -1 and len(peptide) >= N:
        for i in range(0,len(peptide)-N+1,1):
            mer = peptide[i:i+N]
            merArray.append(mer)
    if starIndex >= N:
        peptideTrun = peptide[:starIndex]
        for j in range(0,len(peptideTrun)-N+1,1):
            mer = peptideTrun[j:j+N]
    return merArray
                
        
        

def transcript2peptide(cdna_sequence):
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
    except: full_transcript_store = []
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
    return full_transcript_store

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
    dfNeoJunction = df[((df['dPSI'] >= 0.3) & (df['avg-Healthy__U2AF1-CV'] <= 0.23))]
    return dfNeoJunction  
    
def uid(df, i):
    uid = list(df['UID'])[i]       
    gene = uid.split(':')[0]
    dict = {}
    gene = gene + ':' + uid.split('|')[1].split(':')[0] # slicing the ENSG in background event
    x = uid.split('|')
    dict[gene] = [x[0].split(':')[2]]
    try: x[1].split(':')[2]
    except IndexError: dict[gene].append(x[1].split(':')[1])
    else: 
        fusionExon = str(x[1].split(':')[1])+':'+str(x[1].split(':')[2])
        dict[gene].append(fusionExon)
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
                    exon_seq = 'UUUUUUUUUUUUUUUUUUUUUUUU'   
                    print('{0} observes an UTR event {1}'.format(EnsID,subexon))
                else:
                    if flag == 'site2':           
                        exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],EnsID,attrs[1])  # chr,strand, start,end
                    elif flag == 'site1':
                        exon_seq = query_from_dict_fa(dict_fa,attrs[2],suffix,EnsID,attrs[1])
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

if __name__ == "__main__":
    # load necessary input files
    df = pd.read_csv('PSI.AML__U2AF1-CV_vs_Healthy__U2AF1-CV.txt',sep='\t') 
    df_exonlist = pd.read_csv('mRNA-ExonIDs.txt',sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    dict_exonCoords = exonCoords_to_dict('Hs_Ensembl_exon.txt','\t')
    dict_fa = fasta_to_dict('Hs_gene-seq-2000_flank.fa')
    metaBaml = Meta(df) #Instantiate Meta object
    dfNeoJunction = neoJunctions(metaBaml.df)
    NeoJBaml = NeoJ(dfNeoJunction) #Instantiate NeoJ object
    NeoJBaml.retrieveJunctionSite()
    NeoJBaml.matchWithExonlist(df_exonlist,dict_exonCoords)
    NeoJBaml.getORF()
    NeoJBaml.getORFaa()
    NeoJBaml.phaseTranslateJunction()
    NeoJBaml.get9mer()
    
    dictExonList = GTEx.convertExonList(df_exonlist)  
    dictGTEx = GTEx.backgroundGTEx()
    
    NeoJBaml.getFinalMers()
    
    
    NeoJBaml.df.to_csv('NeoJunction.txt',sep='\t',header=True,index = False)
#    
#    def toFasta(list_):
#    with open('queryNetMHC.fa','w') as file1:
#        for index,item in enumerate(list_):
#            file1.write('>{0} mer\n'.format(index+1))
#            file1.write('{0}\n'.format(item.strip('\'')))
#
#
#list1 = ['LRTSSSWEM', 'ILRTSSSWE', 'RTSSSWEMC', 'CMWILRTSS', 'MWILRTSSS', 'SCMWILRTS', 'WILRTSSSW']
#list2 = ['WLHKIGLVV', 'HIAVGQQMN', 'QQMNLHWLH', 'VGQQMNLHW', 'HKIGLVVIL', 'ALCHIAVGQ', 'LALCHIAVG', 'AVGQQMNLH', 'VLALCHIAV', 'IGLVVILAS', 'GQQMNLHWL', 'VVILASTVV', 'LCHIAVGQQ', 'CHIAVGQQM', 'VILASTVVA', 'NLHWLHKIG', 'MNLHWLHKI', 'HWLHKIGLV', 'ILASTVVAM', 'KIGLVVILA', 'LHWLHKIGL', 'QMNLHWLHK', 'LVVILASTV', 'IAVGQQMNL', 'GLVVILAST', 'LHKIGLVVI']
#list3 = list1 + list2
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

