#!/Users/ligk2e/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 12:27:22 2020

@author: ligk2e
"""
import os
#os.chdir('/Users/ligk2e/Desktop/project_test')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def scratchPlusView1():     # directly use pool for multiprocessing

    global dicSRA

    sraTable1 = pd.read_csv('./data/SraRunTable-GTEX1.txt',sep=',')
    sraTable2 = pd.read_csv('./data/SraRunTable-GTEX2.txt',sep=',')
    sraTable3 = pd.read_csv('./data/SraRunTable-GTEX3.txt',sep=',')
    sraData = pd.read_csv('./data/GTEx_EventAnnotation.txt',sep='\t')
    
    dicSRA = {}
    for i in range(sraTable1.shape[0]):
        accID, tissue= sraTable1['Run'].tolist()[i],sraTable1['body_site'].tolist()[i]
        try: dicSRA[tissue].append(accID)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(accID)
    
    for i in range(sraTable2.shape[0]):
        accID, tissue= sraTable2['Run'].tolist()[i],sraTable2['body_site'].tolist()[i]
        try: dicSRA[tissue].append(accID)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(accID)
        
    for i in range(sraTable3.shape[0]):
        accID, tissue= sraTable3['Run'].tolist()[i],sraTable3['body_site'].tolist()[i]
        try: dicSRA[tissue].append(accID)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(accID)  
            
    
    
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
        
    Spliter1 = Spliter(sraData) 
    Spliter1.split_df(20)  
    liaisonData = Spliter1.queue 
    
    
    import multiprocessing as mp
    p = mp.Pool(processes = 20)
    dicts = p.map(constructDic,liaisonData)
    
    def merge_dicts(dic_list):
    
        result = {}
        for dictionary in dic_list:
            result.update(dictionary)
        return result
    
    dicTissueExp = merge_dicts(dicts)
    return dicTissueExp

def constructDic(sraData):

    dicTissueExp = {}
    for i in range(sraData.shape[0]):
        print('this is the {0}th run of process{1}'.format(i,os.getpid()))
        event = sraData['UID'].tolist()[i]
        for tissue,accID in dicSRA.items():
            try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
            except KeyError:            
                dicTissueExp[event] = {}
                dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
    return dicTissueExp
    






def scratchPlusView2():     # use manager.queue for shared memory
    sraTable1 = pd.read_csv('./data/SraRunTable-GTEX1.txt',sep=',')
    sraTable2 = pd.read_csv('./data/SraRunTable-GTEX2.txt',sep=',')
    sraTable3 = pd.read_csv('./data/SraRunTable-GTEX3.txt',sep=',')
    sraData = pd.read_csv('./data/GTEx_EventAnnotation.txt',sep='\t')
    
    dicSRA = {}
    for i in range(sraTable1.shape[0]):
        accID, tissue= sraTable1['Run'].tolist()[i],sraTable1['body_site'].tolist()[i]
        try: dicSRA[tissue].append(accID)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(accID)
    
    for i in range(sraTable2.shape[0]):
        accID, tissue= sraTable2['Run'].tolist()[i],sraTable2['body_site'].tolist()[i]
        try: dicSRA[tissue].append(accID)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(accID)
        
    for i in range(sraTable3.shape[0]):
        accID, tissue= sraTable3['Run'].tolist()[i],sraTable3['body_site'].tolist()[i]
        try: dicSRA[tissue].append(accID)
        except KeyError:
            dicSRA[tissue] = []
            dicSRA[tissue].append(accID)  
            
    
    
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
        
    Spliter1 = Spliter(sraData) 
    Spliter1.split_df(4)  
    liaisonData = Spliter1.queue 
      
    import multiprocessing as mp
    
    class ParallelWorker4GeneratingGTExDictionary(object):
        def __init__(self,n):  # n means how many process you specify
            self.n = n
            self.manager = mp.Manager()
            self.queue = self.manager.JoinableQueue()
            self.processes = []
       
        @staticmethod    
        def constructDic(sraData,queue,index):
            dicTissueExp = {}
            for i in range(sraData.shape[0]):
                print('this is the {0}th run of process{1}'.format(i,index))
                event = sraData['UID'].tolist()[i]
                for tissue,accID in dicSRA.items():
                    try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
                    except KeyError:            
                        dicTissueExp[event] = {}
                        dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
            queue.put(dicTissueExp)
    
        def conf(self):
            global liaisonData
    
    
            for i in range(self.n):
                p = mp.Process(target=ParallelWorker4GeneratingGTExDictionary.constructDic, args=(liaisonData[i],self.queue,i+1,))
                self.processes.append(p)
                
            self.queue.join()
            
        def run(self):
            [x.start() for x in self.processes]
            [x.join() for x in self.processes]
            
        def collect(self):
            self.result = []
            for i in range(self.n):
                self.result.append(self.queue.get())
    
    
             
    multiWorker = ParallelWorker4GeneratingGTExDictionary(4)
    multiWorker.conf()
    multiWorker.run()
    multiWorker.collect()
            
        
    def merge_dicts(dic_list):
    
        result = {}
        for dictionary in dic_list:
            result.update(dictionary)
        return result
    
    dicTissueExp = merge_dicts(multiWorker.result)
    return dictissueExp

def loadPlusView():
    from time import process_time
    start = process_time()
    import bz2
    import _pickle as cpickle
    import pickle
    with bz2.BZ2File('./data/dicTissueExp2.pbz2','rb') as f1:
         dicTissueExp = cpickle.load(f1)  
    end = process_time()
    
    print('consume {0}'.format(end-start))
    return dicTissueExp


def inspectGTEx(dicTissueExp,event,tissue='all',plot=True):
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
                exp = exp.astype('float64')
                exp=exp[np.logical_not(np.isnan(exp))]
                if exp.size == 0: print('{0} data incomplete'.format(tis))
                elif np.any(exp):   # have non-zero element
                    if plot==True:
                        fig = plt.figure()
                        plt.bar(np.arange(len(exp)),exp,width=0.2,label=tis)
                        plt.xticks(np.arange(len(exp)),np.arange(len(exp))+1)
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
            exp = expression.astype('float64')
            exp=exp[np.logical_not(np.isnan(exp))]
            if exp.size == 0: print('{0} data incomplete'.format(tissue))
            elif np.any(exp):   # have non-zero element
                plt.bar(np.arange(len(exp)),exp,width=0.2,label='tissue')
                plt.legend()
                if not os.path.exists('./GTEx'): os.makedirs('./GTEx')
                plt.savefig('./{}.pdf'.format(tissue),bbox_inches='tight')
                plt.show()
                print(expression)
  

def usage():


    print('Usage:')
    print('python3 queryGTEx.py -e KYAT3:ENSG00000137944:E4.1-I4.1_88965413 -b ENSG00000137944:E4.1-E5.1 -m savage -t Prostate -p True')
    print('Options:')
    print('-e --event : Splicing event you want to interrogate')
    print('-b --back: background event')
    print('-m --mode : load off-the-shelf GTEx data or generating it from scratch')
    print('-t --tissue : tissue-specific abundance you want to inspect')
    print('-p --plot : Do you want to plot all tissue-specific expression bar charts?')
    print('-h --help : check help information ')
    print('Author: Guangyuan(Frank) Li <li2g2@mail.uc.edu>, PhD Student, University of Cincinnati, 2020') 

def main(event,mode,tissue='all',plot=True):
    if mode == 'denovo':
        dicTissueExp = scratchPlusView1()
    if mode == 'savage':
        print('Please wait for around 10 min to load the GTEx Data')
        dicTissueExp = loadPlusView()
    flag = inspectGTEx(dicTissueExp,event,tissue,plot)
    return flag

if __name__ == '__main__':
    import getopt
    import sys
#    log_err = open('queryGTEx.stderr.log','a')
#    log_out = open('queryGTEx.stdout.log','a')
#    sys.stderr = log_err
#    sys.stdout = log_out
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'he:b:m:t:p:',['help','event=','back=','mode=','tissue=','plot='])
    except getopt.GetoptError as err:
        print('ERROR:', err)
        usage()
        sys.exit(1)
    for opt, arg in options:
        if opt in ('-e','--event'):
            event = arg
            print('Queried examined splicing event:', arg)
        elif opt in ('-b','--back'):
            back = arg
            print('Background event:', arg)
        elif opt in ('-m','--mode'):
            mode = arg
            print('Choose the mode to obtain GTEx data',arg)
        elif opt in ('-t','--tissue'):
            tissue = arg
            print('Tissue I want to inspect:', arg)
        elif opt in ('-p','--plot'):
            plot = bool(arg)
            print('Generating plot or not:',arg)
        elif opt in ('--help','-h'):
            usage() 
            sys.exit() 
    full = event + '|' + back        
    main(full,mode,tissue,plot)
    
    





















