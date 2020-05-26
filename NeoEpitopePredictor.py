#!/Users/ligk2e/opt/anaconda3/envs/python3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 19:43:43 2020

@author: ligk2e
"""
import sys
#sys.path.append('/Users/ligk2e/opt/anaconda3/envs/python3/lib/python3.7/site-packages')
import os
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
import matplotlib.pyplot as plt
import math
import urllib.parse
import urllib.request
import requests
import xmltodict
import subprocess
import argparse
import getopt
import ast


############################################################################################
# part1: branch2.py, match to existing peptides
#############################################################################################



def GetIncreasedPart(df):
    df['sign'] = df['dPSI'].apply(lambda x: True if x>0 else False) # x = lambda a,b:a+b; x(5,6)   
    df_ori = df[df['sign']==True]
    df_ori = df_ori.drop(columns=['sign'])  # how to drop a column
    return df_ori




def UID(df, i):
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

def match_with_exonlist(df_ori,df_exonlist,dict_exonCoords,dict_fa):
    col1 = []
    col2 = []
    
    for i in range(df_ori.shape[0]):
        temp=UID(df_ori,i)
        EnsID=list(temp.keys())[0].split(':')[1]
        Exons_examined = exon_extract(temp,0,EnsID)
        Exons_back = exon_extract(temp,1,EnsID)
        col1.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_examined,dict_fa))
        col2.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_back,dict_fa))
        
    return col1,col2


def exon_update(temp,pos,EnsID):
    Exons_former = list(temp.values())[0][pos].split('-')[0]
    #print(Exons_former)
    if len(Exons_former) > 7:  # non-canonical/novel splicing sites
        print(Exons_former + ' in ' + EnsID + ' is a non-canonocal case\n')
        Exons_former_update = Exons_former
    elif Exons_former.startswith('E'):
        Exons_former_nume = Exons_former.lstrip('E')
        Exons_former_nume_update = str(D(Exons_former_nume) + D('0.1'))
        Exons_former_update = 'E' + Exons_former_nume_update
    elif Exons_former.startswith('I'):
        Exons_former_nume = Exons_former.lstrip('I')
        Exons_former_nume_update = str(D(Exons_former_nume) + D('0.1'))
        Exons_former_update = 'I' + Exons_former_nume_update
    Exons_latter = list(temp.values())[0][pos].split('-')[1]
    Exons = Exons_former_update + '|' + Exons_latter
    print(Exons)
    return Exons

def exon_extract(temp,pos,EnsID):
    Exons = list(temp.values())[0][pos].split('-')[0] + '|' + list(temp.values())[0][pos].split('-')[1]
    return Exons

def core_match(df_exonlist,dict_exonCoords,EnsID,Exons,dict_fa):
    
    try:
        df_certain = df_exonlist[df_exonlist['EnsGID'] == EnsID]
    except: full_transcript_store = []  # EnsGID is absent in df_exonlist
    full_transcript_store = []
    for item in list(df_certain['Exons']):
        full_transcript=''
        Exons1 = '|' + Exons
        Exons2 = Exons + '|'
        
        if re.search(rf'{re.escape(Exons1)}',item) or re.search(rf'{re.escape(Exons2)}',item) or re.search(rf'{re.escape(Exons)}$',item):
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
    #print(strand)
    if strand == '+':
        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse,generic_dna).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000 # python range/slice doesn't include end point
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1,generic_dna)
        exon_seq = str(s.reverse_complement())
    return exon_seq

##############################################################################################    
# part2: pick_peptide.py, find the most likely ORF

###############################################################################################    

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
    


def transcript2peptide(cdna_sequence):
    # TAA,TGA,TAG
    # tips: find/index and split function && their corresponding re version for multiple stirng
    reading_manners = []
    reading_manners.append(cdna_sequence[0:])
    reading_manners.append(cdna_sequence[1:])
    reading_manners.append(cdna_sequence[2:])
    frag_comp_array = []
    for manner in reading_manners:       
#        frag_array = re.split(r'TAA|TGA|TAG',manner), it is for multiple condition
        pos = []
        for m in re.finditer(r'TAA|TGA|TAG',manner):   # for multiple instances
            if m.start() % 3 == 0:
                pos.append(m.start())
        frag_array = pos_to_frags(pos,manner)
        for frag in frag_array:
            if 'ATG' not in frag or len(frag) == 0:
                continue
            else:
                for n in re.finditer('ATG',frag):
                    if (len(frag) - n.start()) % 3 == 0:
                        frag_comp = frag[n.start():]
                        frag_comp_array.append(frag_comp)
                        break
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
#            else:
#                print('length diffs within 8, add_score is not able to overturn current max one') 
    max_seq_tran = max_seq
    max_seq_aa = str(Seq(max_seq,generic_dna).translate(to_stop=False))
    max_seq_final = [max_seq_tran,max_seq_aa,frag_comp_array]
    return max_seq_final     

def pos_to_frags(pos,sequence):
    frag_array = []
    if pos: #tips: this elegant way to define criteria        
        frag_array.append(sequence[0:pos[0]])
        i = 0
        while i < len(pos)-1:
            frag_array.append(sequence[pos[i]+3:pos[i+1]])
            i += 1
        last_seq = sequence[pos[-1]+3:]
        if not any(codon in last_seq for codon in ['TAA','TAG','TGA']):
            frag_array.append(sequence[pos[-1]+3:])
#    else:
#        pass
    return frag_array
        
    
def final_conversion(col):   # get most likely ORF
    output_array_aa = []
    output_array_tran = []

    for event in col:
        temp_array_aa = []
        temp_array_tran = []

        for transcript in event:
            if transcript == '':
                temp_array_aa.append(transcript)
                temp_array_tran.append(transcript)
            else:
                max_pep = transcript2peptide(transcript)[1]
                max_tran = transcript2peptide(transcript)[0]
                temp_array_aa.append(max_pep)
                temp_array_tran.append(max_tran)
        output_array_aa.append(temp_array_aa)
        output_array_tran.append(temp_array_tran)
    return output_array_tran,output_array_aa           


def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        next(file)
        for line in file:
#            dict_temp={}
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5])
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}
    return dict_exonCoords

            
def find_longest_AA(listAA):
    max=0
    for item in listAA:
        try:
            stop_pos = item.index('*') # return only first occurence
            length = len(item[:stop_pos])
        except ValueError:
            length=len(item)
        if int(length) > max:
            max = int(length)
            max_item = item
    return max_item

#################################################################################################
# part3: following.py   find the junction sites' sequence.   
#################################################################################################

def retrieveJunctionSite(df,dict_exonCoords,dict_fa):
    exam_seq,back_seq = [],[]
    for i in range(df.shape[0]):
        temp = UID(df,i)
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
        
    df['exam_seq'] = exam_seq
    df['back_seq'] = back_seq
    return df
        
        
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
                try:
                    attrs = dict_exonCoords[fusionGeneEnsID][fusionGeneExon]
                except KeyError:
                    exon_seq = '***********************'
                    print('{0} does not include in {1} exonlists'.format(fusionGeneExon,fusionGeneEnsID))
                else:   
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
                    print('{0} observes an UTR event {1}'.format(EnsID,subexon))
                    chrUTR,strandUTR = utrAttrs(EnsID,dict_exonCoords)

                    exon_seq = utrJunction(suffix,EnsID,strandUTR,chrUTR,flag)  

                else:
                    if flag == 'site2':           
                        exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],EnsID,attrs[1])  # chr,strand, start,end
                    elif flag == 'site1':
                        exon_seq = query_from_dict_fa(dict_fa,attrs[2],suffix,EnsID,attrs[1])
    return exon_seq

def utrJunction(site,EnsGID,strand,chr_,flag):  # U0.1_438493849, here 438493849 means the site
    if flag == 'site1' and strand == '+':  # U0.1_438493849 - E2.2
        otherSite = int(site) - 100 + 1   # extract UTR with length = 100
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
    elif flag == 'site1' and strand == '-':    # 438493849 is the coordinates in forward strand
        otherSite = int(site) + 100 - 1 
        #exon_seq = query_from_dict_fa(dict_fa,site,otherSite,EnsGID,strand)    # site, otherSite must be coordinates in forward strand, strand argument will handle the conversion automatically
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
        exon_seq = str(Seq(exon_seq,generic_dna).reverse_complement())
    elif flag == 'site2' and strand == '+':  # E5.3 - U5.4_48374838
        otherSite = int(site) + 100 -1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
    elif flag == 'site2' and strand == '-':
        otherSite = int(site) - 100 + 1
        #print(EnsGID,chr_,site,otherSite)
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
        exon_seq = str(Seq(exon_seq,generic_dna).reverse_complement())
    return exon_seq

def utrAttrs(EnsID,dict_exonCoords):  # try to get U0.1's attribute, but dict_exonCoords doesn't have, so we just wanna get the first entry for its EnsGID
    exonDict = dict_exonCoords[EnsID] 
    attrs = next(iter(exonDict.values()))
    chr_,strand = attrs[0],attrs[1]
    return chr_,strand



def retrieveSeqFromUCSCapi(chr_,start,end):
    url = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment={0}:{1},{2}'.format(chr_,start,end)
    response = requests.get(url)
    status_code = response.status_code
    assert status_code == 200
    try:
        my_dict = xmltodict.parse(response.content)
    except:
        print(chr_,start,end)
        raise Exception('Sorry,please check above printed stuff')
    exon_seq = my_dict['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n','').upper()
    return exon_seq


##################################################################################################
# part4: narrow down to extracellular instances and check if good representative and do seeding alignment   
###################################################################################################
    
def extract_EnsID(df):
    
    UID = list(df['UID'])
    EnsID_array = []
    for item in UID:
        EnsID = item.split('|')[0].split(':')[1]
        EnsID_array.append(EnsID)
    return EnsID_array

def write_list_to_file(list):
    with open('EnsID4query.txt','w') as f1:
        f1.writelines('%s\n' % EnsID_i for EnsID_i in list)
    return None

 


                        
def check_if_good_representative(df):
    outer_condition_array = []
    outer_condition_detail = []
    
    
    for i in range(df.shape[0]):
        all_whole_tran = df['exam_match_whole_tran'].tolist()[i]
        all_ORF_tran = df['exam_match_tran'].tolist()[i]
        junction = df['exam_seq'].tolist()[i].replace(',','')
        condition_array = []
        for idx in range(len(all_whole_tran)):
            whole = all_whole_tran[idx]
            ORF = all_ORF_tran[idx]
            if whole:
                start_ORF = whole.find(ORF)
                end_ORF = start_ORF + len(ORF)
                start_junction = whole
    # we have to apply fuzzy matching, because junction consists of former part and latter part(i.e. E6.3|E8.1)
    # In my way to handle overlapping 1nt, the former one will always stay constant, but latter one in the whole
    # transcript, it might get trimmed but here we don't trim it, so there might be 1 overhang in jucntion seq.
    
                pattern = regex.compile('(%s){d<=1}' % junction) 
                try:
                    start_junction = pattern.search(whole).span()[0]
                except:
                    print(df.iloc[i]['UID'],idx)
                    raise Exception('bug')
                end_junction = pattern.search(whole).span()[1] - 1

                if start_junction <= end_ORF and end_junction >= start_ORF:
                    condition_array.append(True)
                else:
                    condition_array.append(False)
            else:
                condition_array.append(False)
        if sum(condition_array) == 0: outer_condition_array.append(False)
        else: outer_condition_array.append(True)
        outer_condition_detail.append(condition_array)
    
    df['involvement'] = outer_condition_detail
    df['getTranslated'] = outer_condition_array
    df_retained = df[df['getTranslated']==True]
    df_retained = df_retained.drop(columns=['getTranslated']) 
    df_filtered = df[df['getTranslated']==False] 
    df_filtered = df_filtered.drop(columns=['getTranslated'])
     # filter out events that can not match with existing ones, or events that could match but splicing site won't involve in ORF formation
    return df_retained,df_filtered
            
def alignment_to_uniprot(df,dict_uni_fa,Ens2ACC,mode):
    col1 = []
    col2 = []
    col3 = []
    for i in range(df.shape[0]):
        # collect all uniprot-curated isoform protein sequence
        target = {}
        EnsID = list(df['UID'])[i].split('|')[0].split(':')[1]
        ACCID = Ens2ACC[EnsID] 
        isoforms = list(dict_uni_fa.keys())  # ['Q9NR97','Q9NR97-2'...]
        for iso in isoforms:        
            if ACCID in iso: 
                seq = dict_uni_fa[iso]
                target[iso] = seq
        # collect all mine-predicted isoform protein sequence
        involve = df['involvement'].tolist()[i]  # [True,False,False,True] indicate which transcript would be a good representative
        match_aa = df['exam_match_aa'].tolist()[i]
        repre = []
        for idx,aa in enumerate(match_aa):
            if involve[idx] == True:   # only consider aa that is caused by splicing event
                bucket = chop_sequence(aa,10)   # chopping
                subnotes = {}
                for j in range(len(bucket)):   # for each 10mer
                    frag = bucket[j]
                    for key,value in target.items():   # for each curated isoform
                        try: 
                            subnotes[key].append(True) if frag in value else subnotes[key].append(False)
                        except KeyError:
                            subnotes[key] = []
                            subnotes[key].append(True) if frag in value else subnotes[key].append(False)
                for k,m in subnotes.items():   # k is key, m is value
                    if sum(m)==0: subnotes[k].append('notAligned')
                    elif sum(m)==len(m): subnotes[k].append('totallyAligned')
                    else: subnotes[k].append('partiallyAligned')
                repre.append(subnotes)
            elif involve[idx] == False:
                repre.append('Either splicing site is not involved in ORF formation or splicing event doesn\'t occur in this certain transcript')
        col1.append(repre)
        # define what kind of peptide this splicing sites would generate by interogratting each repre list
        identity = []
        for n in repre:
            definition = ''
            if isinstance(n,dict):
                for p in n.values():   #n will be {'P14061':[True,True,False,'partiallyAligned'],'P14061-2':[True,True,False,'partiallyAligned']}
                    if p[-1] == 'totallyAligned':  
                        definition = 'one of already documented isoforms'
                        break
            else: 
                definition = 'Either splicing site is not involved in ORF formation or splicing event doesn\'t occur in this certain transcript'
            if definition == '': definition = 'novel isoform'
            identity.append(definition)
        col2.append(identity)
        # let's see if it is possible to generate any novel isoform


        if mode=='strigent':   # as long as one of possible concatenation results in totallyAligned, then we stop pursuing

            for idx,w in enumerate(identity):
                crystal = True   
                if w == 'one of already documented isoforms': 
                    crystal = False
                    break
            if crystal == True:
                crystal_ = []
                for idx,w in enumerate(identity):
                    if w == 'novel isoform': 
                        query_aa = match_aa[idx]
                        try: result = TMHMM(query_aa,'{0}_{1}_{2}'.format(i,idx,EnsID))  # there's some rare case, it can match up with exising one,
                        # but it doesn't have reasonable ORFs, so the returned ORF prediction is ''
                        except: result = False
                        if result:
                            crys = (True,idx)  # we need look into that
                            crystal_.append(crys)
                if crystal_ == []: col3.append(False) # no need to look into this event anymore
                else: col3.append(crystal_)
            else:
                col3.append(False)
            
        elif mode=='loose':    # consider every possible novel isoform possibilities
            crystal = []
            for idx,w in enumerate(identity):
                if w == 'novel isoform': 
                    query_aa = match_aa[idx]
                    try: result = TMHMM(query_aa,'{0}_{1}_{2}'.format(i,idx,EnsID))  # there's some rare case, it can match up with exising one,
                    # but it doesn't have reasonable ORFs, so the returned ORF prediction is ''
                    except: result = False
                    if result:
                        crys = (True,idx)  # we need look into that
                        crystal.append(crys)
            if crystal == []: col3.append(False) # no need to look into this event anymore
            else: col3.append(crystal)
            
        
    df['alignment'] = col1
    df['identity'] = col2
    try:df['interest'] = col3
    except: 
        print(col3)
        raise Exception('hi')
    return df

    
          

    
def read_uniprot_seq(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            uniID = title.split('|')[1]
            dict_fa[uniID] = seq
    return dict_fa       

def chop_sequence(seq,kmer):   # how to splice sequence, elegant way to use range
    frag_bucket = []
    for i in range(0,len(seq),kmer):
        try:
            frag_bucket.append(seq[i:i+kmer])
        except:
            frag_bucket.append(seq[i:])
    return frag_bucket

##############################################################################################
# part5: interrogating chromosome stataistics
#############################################################################################
    

def ChroDistribution(df):
    chro_array = []
    for i in range(df.shape[0]):
        ensid = list(df['UID'])[i].split('|')[0].split(':')[1]
        chro = dict_fa[ensid][0]
        chro_array.append(chro)
    freq = collections.Counter(chro_array)
    return freq

'''courtesy by user on stackoverflow'''
def Round2Precision(value,precision:int=0,mode:str=''): # default argument with specified type
    assert precision >= 0 # if true, continue, otherwise raise assertError, using for self-check
    value *= 10 ** precision # if you wanna round by precision, have to do that
    method = round   # round will base on >.5 or <.5
    if mode.lower() == 'up': 
        method = math.ceil     # always round up
    elif mode.lower() == 'down':
        method = math.floor   # always round down
    answer = '{0:.{1}f}'.format(method(value)/10**precision,precision)   
    return float(answer) 


def PlotChroScarse(chro_dict,path):
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.75,0.75])  #[left,bottom, width, height]
    scarse = {}
    for chro,attr in chro_dict.items():
        chro_s = re.split(r'chr',chro)[-1]
        scarse[chro_s] = Round2Precision(attr[0]/attr[1],2)  # genes per 1 Millon bp
    x_axis = list(scarse.keys())
    y_axis = list(scarse.values())
    ax.bar(x_axis,y_axis)
    ax.set(xlabel='chromosome',ylabel='genes per 1 Million bp',title='crowdness of human chromosome')
    
    #ax.legend()    
    #plt.show()
    fig.savefig(path)
    plt.close(fig)
    
def IDmappingACC2Ensembl(lis,mannual=False):  
    if os.path.exists(os.path.join(outFolder,'Ens2ACC.p')): 
        with open(os.path.join(outFolder,'Ens2ACC.p'),'rb') as f2:
            Ens2ACC = pickle.load(f2)
    else:    
        query = ' '.join(lis)
        url = 'https://www.uniprot.org/uploadlists/'    
        params = {
        'from': 'ACC+ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': query
        }
        
        data = urllib.parse.urlencode(params)  # convert dict to string
        data = data.encode('utf-8')     # convert string to byte
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()   
        a = response.decode('utf-8')        # convert byte to string
        
        rows = a.split('\n')
        rows = rows[1:]   # the first item is 'from':'to'
        rows = rows[:-1]  # the last item is empty, it is becasue the last '\n' got recognized as newline
        ACC2Ens,Ens2ACC = {},{}
        for item in rows:
            ACC = item.split('\t')[0]
            Ens = item.split('\t')[1]        
            ACC2Ens[ACC] = Ens
            Ens2ACC[Ens] = ACC
        diff = len(set(lis)) - len(set(list(ACC2Ens.keys())))
        if mannual == False:
            print('There are {0} membrane protein missing!!!!'.format(diff))
        elif mannual == True:
            print('There are {0} membrane protein needs to be mannually checked, you can exit anytime'.format(diff))            
            for acc in lis:     # some ACC doesn't map to a EnsGID, double check see if every one in human membrane list has been mapped
                try: ACC2Ens[acc]
                except KeyError:
                    print('Please mannually check if {0} has corresponding EnsGID:'.format(acc))
                    cond = input('Does {0} have corresponding EnsGID? (y|n|e):'.format(acc))
                    if cond == 'y':
                        EnsGID = input('Please enter the corresponding EnsGID:')
                        ACC2Ens[acc] = EnsGID
                        Ens2ACC[EnsGID] = acc
                        print('Great!, {0} has been mannually added to dictionary.'.format(EnsGID))
                    if cond == 'n':
                        print('Caveat: {0} is a humen membrane protein, it will not be considered in following analysis'.format(acc))
                    if cond == 'e':
                        break
           
            with open(os.path.join(outFolder,'Ens2ACC.p'),'wb') as f1:
                pickle.dump(Ens2ACC,f1)
    return Ens2ACC    # as complete as possible

def TMHMM(aa,name):
    # download TMHMM linux version, untar it.
    # change shabang of tmhmm and tmhmmformat.pl to the path of perl 5+ you loaded
    # export the path of tmhmm to $PATH, finsh configuration

    # in my use case, save those for convenience
    # perl: /usr/local/perl/5.20.1/bin/perl
    # tmhmm: /data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin
    if not os.path.exists(os.path.join(outFolder,'TMHMM_temp')): os.makedirs(os.path.join(outFolder,'TMHMM_temp'))
    with open(os.path.join(outFolder,'TMHMM_temp','{}.fasta'.format(name)),'w') as f1:
        f1.write('>peptide_{}\n'.format(name))
        f1.write(aa)
    with open(os.path.join(outFolder,'TMHMM_temp','{}.out'.format(name)),'w') as f2:
        subprocess.run(['tmhmm',os.path.join(outFolder,'TMHMM_temp','{}.fasta'.format(name))],stdout = f2)
    with open(os.path.join(outFolder,'TMHMM_temp','{}.out'.format(name)),'r') as f3:
        next(f3)
        punchline = f3.readline().rstrip('\n').split(' ')
        TMn = int(punchline[-1])
    result = True if TMn > 0 else False
    return result

def diffNovelFromNotInvolved(df):
    col = []

    for i in range(df.shape[0]):
        cond = True
        exam_match_whole_tran = df['exam_match_whole_tran'].iloc[i]   # will be a list already, 1*1
        for item in exam_match_whole_tran:
            if item: 
                cond = False
                break
        col.append(cond)
    df['cond'] = col
    try:df_novel = df[df['cond']]
    except:
        print(len(col),df.shape[0],col)
        df.to_csv(os.path.join(outFolder,'fjdfjd.txt'),sep='\t',index=None)
        raise Exception('jkjk')
    return df_novel




        
def main(intFile,dataFolder,outFolder,mode):
#    intFile = parser.intFile
#    dataFolder = parser.dataFolder
#    outFolder = parser.outFolder
#    print(intFile,dataFolder,outFolder)
    # get increased part
    df = pd.read_csv(intFile,sep='\t')
    
    # load the files
    df_ori = GetIncreasedPart(df)
    df_exonlist = pd.read_csv(os.path.join(dataFolder,'mRNA-ExonIDs.txt'),sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    dict_exonCoords = exonCoords_to_dict(os.path.join(dataFolder,'Hs_Ensembl_exon.txt'),'\t')
    
    dict_fa = fasta_to_dict(os.path.join(dataFolder,'Hs_gene-seq-2000_flank.fa'))
    
    # overlapping with human membrane proteins
    df_membraneProteins = pd.read_csv(os.path.join(dataFolder,'human_membrane_proteins.txt'),sep='\t')  
    ACClist = df_membraneProteins['Entry'].tolist()
    Ens2ACC = IDmappingACC2Ensembl(ACClist,True) 
    EnsID = extract_EnsID(df_ori)
    df_ori['condition'] = [True if item in list(Ens2ACC.keys()) else False for item in EnsID]
    df_ori_narrow = df_ori[df_ori['condition'] == True]
    df_ori_narrow = df_ori_narrow.drop(columns=['condition'])

    # match with all existing transcripts
    col1,col2 = match_with_exonlist(df_ori_narrow,df_exonlist,dict_exonCoords,dict_fa)
    
    df_ori_narrow['exam_match_whole_tran'] = col1
    df_ori_narrow['back_match_whole_tran'] = col2
    
    # derive the most likely ORF for each whole transcript sequence
    output_exam_tran,output_exam_aa = final_conversion(col1)
    output_back_tran,output_back_aa = final_conversion(col2)    
    
    
    df_ori_narrow['exam_match_tran'] = output_exam_tran
    df_ori_narrow['exam_match_aa'] = output_exam_aa
    df_ori_narrow['back_match_tran'] = output_back_tran
    df_ori_narrow['back_match_aa'] = output_back_aa
    
    # derive the junction site sequence and add two columns to df_ori
    new_df_narrow = retrieveJunctionSite(df_ori_narrow,dict_exonCoords,dict_fa)
    
    df_retained,df_filtered = check_if_good_representative(new_df_narrow)
    
    # for retained ones

    
    ### final alignment
    dict_uni_fa = read_uniprot_seq(os.path.join(dataFolder,'uniprot_isoform.fasta'))
        
    df_retained_aligned = alignment_to_uniprot(df_retained,dict_uni_fa,Ens2ACC,mode)   
    df_retained_aligned.to_csv(os.path.join(outFolder,'df_retained.txt'),sep='\t',index=None)
    
    # for filtered ones
    df_novel = diffNovelFromNotInvolved(df_filtered)
    df_novel.to_csv(os.path.join(outFolder,'df_novel.txt'),sep='\t',index=None)    

def usage():
    print('Usage:')
    print('python3 NeoEpitopePredictor.py --intFile /data/salomonis2/LabFiles/Frank-Li/breast_cancer_run/receptor/pseudoPSImatrix_percentage.txt --dataFolder /data/salomonis2/LabFiles/Frank-Li/python3/data --outFolder /data/salomonis2/LabFiles/Frank-Li/breast_cancer_run/receptor --mode strigent')
    print('Options:')
    print('--intFile : path of input file')
    print('--dataFolder : path of data folder')
    print('--outFolder : output folder')    
    print('--mode : using TMHMM or not')
    
    

if __name__ == "__main__":
    #os.chdir('/Users/ligk2e/Desktop/project_LUAD')
    
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'h',['help','intFile=','dataFolder=','outFolder=','mode='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('--intFile'):
            intFile = arg
            print('Input file is:', arg)
        elif opt in ('--dataFolder'):
            dataFolder = arg
            print('Data folder is:',arg)
        elif opt in ('--outFolder'):
            outFolder = arg
            print('output folder:',arg)
        elif opt in ('--mode'):
            mode = arg

            print('Using TMHMM?:',arg)
        elif opt in ('-h','--help'):
            usage()
            sys.exit(1)

    main(intFile,dataFolder,outFolder,mode)
    
    
    
    
#    parser = argparse.ArgumentParser(description='Receptor Protein Prediction') # ArgumentParser object
#    parser.add_argument('--intFile',type=str,default='.',help='input file path')
#    parser.add_argument('--dataFolder',type=str,default='./data',help='data folder path')
#    parser.add_argument('--outFolder',type=str,default='.',help='output folder path')
#    args = parser.parse_args()   # namespace object
#    main(args)


    
    # summarize the distribution of splicing event in df_all and df_increased
#    freq_all = ChroDistribution(df)
#    freq_increased = ChroDistribution(df_ori)
#    print(freq_all,freq_increased)
    
    # continue exploit on chromosomes
#    chro_dict = {
#            'chr1': [1961,248,'Metacentric'],    #[1961genes,248 or so million bp, type of centromere]
#            'chr2': [1194,242,'Submetacentric'],
#            'chr3': [1024,198,'Metacentric'],
#            'chr4': [727,190,'Submetacentric'],
#            'chr5': [839,181,'Submetacentric'],
#            'chr6': [996,170,'Submetacentric'],
#            'chr7': [862,159,'Submetacentric'],
#            'chr8': [646,145,'Submetacentric'],
#            'chr9': [739,138,'Submetacentric'],
#            'chr10': [706,133,'Submetacentric'],
#            'chr11': [1224,135,'Submetacentric'],
#            'chr12': [988,133,'Submetacentric'],
#            'chr13': [308,114,'Acrocentric'],
#            'chr14': [583,107,'Acrocentric'],
#            'chr15': [561,101,'Acrocentric'],
#            'chr16': [795,90,'Metacentric'],
#            'chr17': [1124,83,'Submetacentric'],
#            'chr18': [261,80,'Submetacentric'],
#            'chr19': [1357,58,'Metacentric'],
#            'chr20': [516,64,'Metacentric'],
#            'chr21': [215,46,'Acrocentric'],
#            'chr22': [417,50,'Acrocentric'],
#            'chrX': [804,156,'Submetacentric'],
#            'chrY': [63,57,'Acrocentric']}  
#    PlotChroScarse(chro_dict,'human chromosome genes distribution.pdf')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    