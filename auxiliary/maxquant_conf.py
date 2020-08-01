#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 10:36:00 2020

@author: ligk2e
"""

import xmltodict
import copy


def add_database_file(doc,intFile):
    with open(intFile,'r') as f:
        for line in f:

            line = line.rstrip('\n')
            #print(line)
            try:
                template = copy.deepcopy(doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'][-1]) # non-first time
            except:
                template = copy.deepcopy(doc['MaxQuantParams']['fastaFiles']['FastaFileInfo']) # first time
            finally:
                template['fastaFilePath'] = line
                template['identifierParseRule'] = '>(.*)'
                template['descriptionParseRule'] = '>(.*)'
                #print(template)
                try:
                    doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'].append(template)
                except:
                    doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'] = []
                    doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'].append(template)
                    #print(doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'])
    return doc


def numThreads(doc,thread):
    doc['MaxQuantParams']['numThreads'] = str(thread)
    return doc

def add_input_files(doc,intFile):
    with open(intFile,'r') as f:
        experiment = 1
        for line in f:
            line = line.rstrip('\n')
            # 'filePaths'
            if not type(doc['MaxQuantParams']['filePaths']['string']) == list:
                doc['MaxQuantParams']['filePaths']['string'] = []
                doc['MaxQuantParams']['filePaths']['string'].append(line)
            else:
                doc['MaxQuantParams']['filePaths']['string'].append(line)
            
            # 'experiments'
            if not type(doc['MaxQuantParams']['experiments']['string']) == list:
                doc['MaxQuantParams']['experiments']['string'] = []
                doc['MaxQuantParams']['experiments']['string'].append(experiment)
            else:
                doc['MaxQuantParams']['experiments']['string'].append(experiment)
            experiment += 1
            
            # 'fractions'
            if not type(doc['MaxQuantParams']['fractions']['short']) == list:
                doc['MaxQuantParams']['fractions']['short'] = []
                doc['MaxQuantParams']['fractions']['short'].append('32767')
            else:
                doc['MaxQuantParams']['fractions']['short'].append('32767')

             
            # 'ptms'
            if not type(doc['MaxQuantParams']['ptms']['boolean']) == list:
                doc['MaxQuantParams']['ptms']['boolean'] = []
                doc['MaxQuantParams']['ptms']['boolean'].append('False')
            else:
                doc['MaxQuantParams']['ptms']['boolean'].append('False')
                
            # 'paramGroupIndices'
            if not type(doc['MaxQuantParams']['paramGroupIndices']['int']) == list:
                doc['MaxQuantParams']['paramGroupIndices']['int'] = []
                doc['MaxQuantParams']['paramGroupIndices']['int'].append('0')
            else:
                doc['MaxQuantParams']['paramGroupIndices']['int'].append('0')
                
            # 'referenceChannel'
            if not type(doc['MaxQuantParams']['referenceChannel']['string']) == list:
                doc['MaxQuantParams']['referenceChannel']['string'] = []
                doc['MaxQuantParams']['referenceChannel']['string'].append(None)
            else:
                doc['MaxQuantParams']['referenceChannel']['string'].append(None)
    
    return doc
    
                
            
            
def change_enzymes(doc,enzymes,mode):

    '''
    mode 0: unspecific
    mode 3: semi-specific
    mode 4: unspecific
    mode 5: no digestion

    '''

    doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymeMode'] = str(mode)

    if enzymes == None:
        doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes'] = None
    else:
        for enzyme in enzymes:
            if not type(doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string']) == list:
                doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string'] = []
                doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string'].append(enzyme)
            else:
                doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string'].append(enzyme)
     
    return doc           





def change_fdr(doc,protein_fdr,peptide_fdr,site_fdr):
    doc['MaxQuantParams']['proteinFdr'] = protein_fdr
    doc['MaxQuantParams']['peptideFdr'] = peptide_fdr
    doc['MaxQuantParams']['siteFdr'] = site_fdr
    return doc


def change_contaminants(doc,include):
    doc['MaxQuantParams']['includeContaminants'] = str(include)
    return doc



def change_length(doc,minPepLen,maxPeptideMass,minPeptideLengthForUnspecificSearch,maxPeptideLengthForUnspecificSearch):
    doc['MaxQuantParams']['minPepLen'] = minPepLen
    doc['MaxQuantParams']['maxPeptideMass'] = maxPeptideMass
    doc['MaxQuantParams']['minPeptideLengthForUnspecificSearch'] = minPeptideLengthForUnspecificSearch
    doc['MaxQuantParams']['maxPeptideLengthForUnspecificSearch'] = maxPeptideLengthForUnspecificSearch
    return doc





if __name__ == '__main__':

    with open('/data/salomonis2/LabFiles/Frank-Li/python3/mqpar_iTRAQ_4plex.xml','r') as fd:
        doc = xmltodict.parse(fd.read())   # default mqpar.xml file
    
    doc_1 = add_database_file(doc,'/data/salomonis2/LabFiles/Frank-Li/CPTAC/TCGA_breast/TCGA_E2-A10A-01/raw_nodigest_fdr10/database.txt')
    doc_2 = numThreads(doc_1,10)
    doc_3 = add_input_files(doc_2,'/data/salomonis2/LabFiles/Frank-Li/CPTAC/TCGA_breast/TCGA_E2-A10A-01/raw_nodigest_fdr10/intFile.txt')
    
    
    #doc_4 = change_enzymes(doc_3,['Trypsin','LysC'],mode=3)
    doc_4 = change_enzymes(doc_3,None,mode=5)
    doc_5 = change_fdr(doc_4,protein_fdr=1,peptide_fdr=0.1,site_fdr=1)
    doc_6 = change_contaminants(doc_5,True)
    doc_7 = change_length(doc_6,8,3000,8,15)
    a = xmltodict.unparse(doc_4,pretty=True,encoding='utf-8')
    a = a.replace('&gt;','>')
    with open('/data/salomonis2/LabFiles/Frank-Li/CPTAC/TCGA_breast/TCGA_E2-A10A-01/raw_nodigest_fdr10/mqpar.xml','w') as f1:
        f1.write(a)
  























  