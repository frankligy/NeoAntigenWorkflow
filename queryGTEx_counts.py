#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:41:02 2020

@author: ligk2e
"""

import os
#os.chdir('/Users/ligk2e/Desktop/project_test')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import hickle as hkl




    
class Spliter():
    def __init__(self,good):   # good is seperable based on some index
        self.good = good
        self.queue = []  # to store splited data
    
    def split_df(self,n):
        dim  = self.good.shape[0]
        lis = [i for i in range(dim)]
        size = len(lis)//n + 1
        for j in range(0,len(lis),size):
            part_index = lis[j:j+size]
            part_df = self.good.iloc[part_index]
            
            self.queue.append(part_df)

    def split_ndarray(self,n):
        dim = len(self.good)
        lis = [i for i in range(dim)]
        size = len(lis)//n + 1
        for j in range(0,len(lis),size):
            part_ndarray = self.good[j:j+size]
            self.queue.append(part_ndarray)



def scratchPlusView1():     # directly use pool for multiprocessing

    global dicSRA

    sraTable = pd.read_csv('./data/GTEx_SRARunTable.txt',sep='\t')
    
    
    sraData_c = pd.read_csv('./data/Hs_RNASeq_Group-1_vs_Group-2.txt',sep='\t')

    colname = sraData_c.columns.tolist()[1:]   # first colun is 'UID'
    colname1 = [item.split(':')[1] for item in colname]   # trim the extra 'Group1,2'
    colname2 = ['UID']+colname1
    sraData_c.columns = colname2   # change the colname to simplified version
    
    sraData_c.index = sraData_c['UID']
    UID = sraData_c['UID'].values    
    truth = np.array([True if ':' in uid else False for uid in UID])
    valid = UID[truth]   # ndarray
    sraData_c = sraData_c.loc[valid]
    
    
    
    conversion = sraTable[['Run','body_site']]   # this is a complete table
    
    analyzed = sraData_c.columns.tolist()[1:]
    SRR_ana = [item.split('_')[0] for item in analyzed]    # out data only has 1274 samples
    
    
    
    conversion.index = conversion['Run'].tolist()
    
    new_conversion = conversion.loc[list(set(conversion['Run'].tolist()).intersection(set(SRR_ana)))]  # take intersection
    
    
    
    dicSRA = {}
    for i in range(new_conversion.shape[0]):
        SRR = new_conversion.iloc[i]['Run']
        tissue = new_conversion.iloc[i]['body_site']
        try: 
            dicSRA[tissue].append(SRR)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(SRR) 
            
    

        
    Spliter1 = Spliter(sraData_c) 
    Spliter1.split_df(40)  
    liaisonData = Spliter1.queue 
    
    
    import multiprocessing as mp
    p = mp.Pool(processes = 40)
    dicts = p.map(constructDic,liaisonData)
    
    def merge_dicts(dic_list):
    
        result = {}
        for dictionary in dic_list:
            result.update(dictionary)
        return result
    
    dicTissueExp = merge_dicts(dicts)

    hkl.dump(dicTissueExp, 'dicTissueExp_c_gzip.hkl', mode='w', compression='gzip')
    return dicTissueExp


def constructDic(sraData):

    dicTissueExp = {}
    for i in range(sraData.shape[0]):
        print('this is the {0}th run of process{1}'.format(i,os.getpid()))
        event = sraData['UID'].tolist()[i]
        for tissue,accID in dicSRA.items():
            try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values  # here PSI value will be stored as ndarray
            except KeyError:            
                dicTissueExp[event] = {}
                dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
    return dicTissueExp

    


def loadPlusView():
    from time import process_time
    start = process_time()
    import bz2
    import _pickle as cpickle
    import pickle
    with bz2.BZ2File('./data/dicTissueExp_c.pbz2','rb') as f1:
         dicTissueExp = cpickle.load(f1)  
    end = process_time()
    
    print('consume {0}'.format(end-start))
    return dicTissueExp


def inspectGTEx(dicTissueExp,event,cutoff,tissue,plot):
    flag = 0
    import warnings
    warnings.filterwarnings("ignore")

    if tissue=='all':
        try:
            tissueExp = dicTissueExp[event]
        except KeyError:
            print('Don\'t detect expression of {0} in normal tissue'.format(event))
        else:

            for tis,exp in tissueExp.items():
                exp = exp.astype('int')
                
                exp[np.isnan(exp)] = 0   # because nan just means they don't even have expression

                if exp.size == 0: print('{0} data incomplete'.format(tis))
                elif np.any(exp):   # have non-zero element
                    
                    hits = sum([True if item > cutoff else False for item in exp]) # how many samples has PSI > 0.1
                    size = exp.size   # how many samples in total

                    print('{0}:({1}/{2}) has read counts > {3}'.format(tis,hits,size,cutoff))
                    
                    if plot=='True':
                        fig = plt.figure()
                        plt.bar(np.arange(len(exp)),exp,width=0.2,label=tis)
                        plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                        plt.xlabel('GTEx Samples')
                        plt.ylabel('PSI value')
                        plt.legend()
                        if not os.path.exists('./GTEx'): os.makedirs('./GTEx')
                        plt.savefig('./GTEx/{1}.pdf'.format(event,tis),bbox_inches='tight')
                        plt.close(fig)
                    else: continue
                else: 
                    flag += 1
                    print('No expression in {}'.format(tis))
            print('{0} has no expression in {1} tissue types'.format(event,flag))
            
    else:
        try:
            expression = dicTissueExp[event][tissue]
        except KeyError:
            print('Don\'t detect expression of {0} in normal tissue'.format(event))
        else:
            exp = expression.astype('int')
            exp[np.isnan(exp)] = 0
            if exp.size == 0: print('{0} data incomplete'.format(tissue))
            elif np.any(exp):   # have non-zero element
                
                hits = sum([True if item > cutoff else False for item in exp]) # how many samples has PSI > 0.1
                size = exp.size   # how many samples in total
                print('{0}:({1}/{2}) has read counts > {3}'.format(tissue,hits,size,cutoff))
                
                fig = plt.figure()
                plt.bar(np.arange(len(exp)),exp,width=0.2,label='tissue')
                plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
                plt.xlabel('GTEx Samples')
                plt.ylabel('PSI value')
                plt.legend()
                if not os.path.exists('./GTEx'): os.makedirs('./GTEx')
                plt.savefig('./{}.pdf'.format(tissue),bbox_inches='tight')
                plt.show()

  

def usage():


    print('Usage:')
    print('python3 queryGTEx_counts.py -e ENSG00000137944:E4.1-I4.1_88965413 -c 5 -m savage -t all -p True')
    print('Options:')
    print('-e --event : Splicing event you want to interrogate')
    print('-c --cutoff : cutoff value for read counts')
    print('-m --mode : load off-the-shelf GTEx data or generating it from scratch')
    print('-t --tissue : tissue-specific abundance you want to inspect')
    print('-p --plot : Do you want to plot all tissue-specific expression bar charts?')
    print('-h --help : check help information ')
    print('Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati, 2020') 

def main(event,mode,cutoff,tissue,plot):
    if mode == 'denovo':
        dicTissueExp = scratchPlusView1()
    if mode == 'savage':
        print('Please wait for around 10 min to load the GTEx Data')
        dicTissueExp = loadPlusView()
    flag = inspectGTEx(dicTissueExp,event,cutoff,tissue,plot)
    return flag

if __name__ == '__main__':
    import getopt
    import sys
#    log_err = open('queryGTEx.stderr.log','a')
#    log_out = open('queryGTEx.stdout.log','a')
#    sys.stderr = log_err
#    sys.stdout = log_out
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'he:c:m:t:p:',['help','event=','cutoff=','mode=','tissue=','plot='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('-e','--event'):
            event = arg
            print('Queried examined splicing event:', arg)
        elif opt in ('-c','--cutoff'):
            cutoff = float(arg)
            print('cutoff value is:', arg)
        elif opt in ('-m','--mode'):
            mode = arg
            print('Choose the mode to obtain GTEx data',arg)
        elif opt in ('-t','--tissue'):
            tissue = arg
            print('Tissue I want to inspect:', arg)
        elif opt in ('-p','--plot'):
            plot = arg
            print('Generating plot or not:',arg)
        elif opt in ('--help','-h'):
            usage() 
            sys.exit()       
    main(event,mode,cutoff,tissue,plot)

