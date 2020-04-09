#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:40:56 2020

@author: ligk2e
"""
import pandas as pd
import pickle
import os


def backgroundGTEx():
    '''
    Run following command in your terminal, get sampleType.txt, wholeBloodID,allsampleID.txt
    
    cut -f 1,7 GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | awk 'BEGIN {FS="\t";} {if (NR>1) {print $1 "\t" $2;}}' > sampleType.txt
    awk '{if ($2 ~/Whole/) {print $1}}' sampleType.txt > wholeBloodID.txt 
    head -n 3 GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct | tail -n 1 | tr '\t' '\n' | awk 'NR>2' > allsampleID.txt
    
    '''
    
    wholeBloodID = list(pd.read_csv('wholeBloodID.txt',sep='\n',header=None,names=['id'])['id'])
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
            if highExpression: merArrayNormal = getNmerNormal(EnsGID,EnsTID)  # merArrayNormal is all 9mer for each transcript that expresses in normal tissue
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
    


def merPassGTExCheck(dictGTEx,uid,merArray):
    EnsGID = uid.split('|')[0].split(':')[1]
    splicingEvent = uid.split('|')[0].split(':')[2:]  # E12.3-E13.4 # E12.4_48334893-E23.4  # E12.4-ENSG00000898765:E3.4
    normalMer = dictGTEx[EnsGID][1]['wholeMerArrayNormal']
    arrangedNormal = flattenNestedList(normalMer,1)
    arrangedSplicing = flattenNestedList(merArray,1,True)
    uniqueInSplicing = list(arrangedSplicing - arrangedNormal)
    return uniqueInSplicing
    
    
    
    

    
    
    
def flattenNestedList(list_,mode=0,clean=False): # mode = 0 means returning flatten list, mode = 1 means converting to set, mode = 2 means converting to list(set())
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




           
def getNmerNormal(EnsID,EnsTID):
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
        merArray = extractNmer(peptide,9)  # all the 9mer the normal
        return merArray   # ['SRTTIFDFJD','FJDJFKDJFKDJ','FJKDFJKDF']
   
  
            
if __name__ == "__main__":
    dictExonList = convertExonList(df_exonlist)  
    dictGTEx = backgroundGTEx()      
    merArray = ['*', ['SCMWILRTS', 'CMWILRTSS', 'MWILRTSSS', 'WILRTSSSW', 'ILRTSSSWE', 'LRTSSSWEM', 'RTSSSWEMC'], ['SCMWILRTS', 'CMWILRTSS', 'MWILRTSSS', 'WILRTSSSW', 'ILRTSSSWE', 'LRTSSSWEM', 'RTSSSWEMC'], ['SCMWILRTS', 'CMWILRTSS', 'MWILRTSSS', 'WILRTSSSW', 'ILRTSSSWE', 'LRTSSSWEM', 'RTSSSWEMC']]    
    
        
    
    
    


