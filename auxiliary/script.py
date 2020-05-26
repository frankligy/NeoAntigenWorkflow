# TASK1: May 9th, try to generate all tissue-specific splicing event expression in GTEx dataset. 
# import pandas as pd

# sraData = pd.read_csv('GTEx_EventAnnotation.txt',sep='\t')

# import pickle
# with open('dicSRA.p','rb') as file1:
#     dicSRA=pickle.load(file1)


# dicTissueExp = {}
# for i in range(sraData.shape[0]):
#     print('this is the {0}th run'.format(i))
#     event = sraData['UID'].tolist()[i]
#     for tissue,accID in dicSRA.items():
#         try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
#         except KeyError:            
#             dicTissueExp[event] = {}
#             dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values

# with open('dicTissueExp.p','wb') as file2:
#     pickle.dump(dicTissueExp,file2)

# TASK2:May 1oth, try to parallelize May 9th code, and use bz2 library to compress pickle object
# import pandas as pd

# sraData = pd.read_csv('GTEx_EventAnnotation.txt',sep='\t')

# import pickle
# with open('dicSRA.p','rb') as file1:
#     dicSRA=pickle.load(file1)

# dicTissueExp = {}
# semaphore = 0

# def search():
#     global semaphore
#     global dicTissueExp
#     while semaphore < sraData.shape[0]:
#         print('this is the {0}th run'.format(semaphore))
#         event = sraData['UID'].tolist()[semaphore]
#         for tissue,accID in dicSRA.items():
#             try: dicTissueExp[event][tissue] = sraData.iloc[semaphore][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
#             except KeyError:          
#                 dicTissueExp[event] = {}
#                 dicTissueExp[event][tissue] = sraData.iloc[semaphore][[accID+'_1.bed' for accID in dicSRA[tissue]]].values 
#         semaphore += 1

# def thread_task(lock): 
#     for _ in range(sraData.shape[0]//2+1): 
#         lock.acquire() 
#         search()
#         lock.release() 


# import threading 
# lock = threading.Lock()
# t1 = threading.Thread(target=thread_task, args=(lock,)) 
# t2 = threading.Thread(target=thread_task, args=(lock,)) 

# t1.start() 
# t2.start() 
# t1.join() 
# t2.join()  

# import bz2
# import _pickle as cpickle
# with bz2.BZ2File('dicTissueExp.pbz2','wb') as f1:
#     cpickle.dump(dicTissueExp,f1)


# TASK3: May 10th, modify TASK2 to multiprocessing version
# import pandas as pd

# sraData = pd.read_csv('GTEx_EventAnnotation.txt',sep='\t')

# import pickle
# with open('dicSRA.p','rb') as file1:
#     dicSRA=pickle.load(file1)

# dicTissueExp = {}
# semaphore = 0

# def search(ns):
#     semaphore = ns.semaphore
#     dicTissueExp = ns.dicTissueExp
#     sraData = ns.sraData
#     dicSRA = ns.dicSRA
#     while semaphore < sraData.shape[0]:
#         print('this is the {0}th run'.format(semaphore))
#         event = sraData['UID'].tolist()[semaphore]
#         for tissue,accID in dicSRA.items():
#             try: dicTissueExp[event][tissue] = sraData.iloc[semaphore][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
#             except KeyError:          
#                 dicTissueExp[event] = {}
#                 dicTissueExp[event][tissue] = sraData.iloc[semaphore][[accID+'_1.bed' for accID in dicSRA[tissue]]].values 
#         semaphore += 1

# def process_task(ns,lock): 
#     for _ in range(sraData.shape[0]//4+1): 
#         lock.acquire() 
#         search(ns)
#         lock.release() 

# import multiprocessing
# with multiprocessing.Manager() as manager:
#     ns = manager.Namespace()
#     ns.sraData = sraData
#     ns.dicSRA = dicSRA
#     ns.dicTissueExp = dicTissueExp
#     ns.semaphore = semaphore



# lock = multiprocessing.Lock()
# p1 = multiprocessing.Process(target=process_task,args=(ns,lock,))
# p2 = multiprocessing.Process(target=process_task,args=(ns,lock,))
# p3 = multiprocessing.Process(target=process_task,args=(ns,lock,))
# p4 = multiprocessing.Process(target=process_task,args=(ns,lock,))

# p1.start()
# p2.start()
# p3.start()
# p4.start()

# p1.join()
# p2.join()
# p3.join()
# p4.join()



# import bz2
# import _pickle as cpickle
# with bz2.BZ2File('dicTissueExp1.pbz2','wb') as f1:
#     cpickle.dump(dicTissueExp,f1)   

# TASK4: May 17th, take up the challenge to multi-process TASK3 again, this time use manager.JoinableQueue()

# import pandas as pd

# sraTable1 = pd.read_csv('./SraRunTable-GTEX1.txt',sep=',')
# sraTable2 = pd.read_csv('./SraRunTable-GTEX2.txt',sep=',')
# sraTable3 = pd.read_csv('./SraRunTable-GTEX3.txt',sep=',')
# sraData = pd.read_csv('./GTEx_EventAnnotation.txt',sep='\t')

# dicSRA = {}
# for i in range(sraTable1.shape[0]):
#     accID, tissue= sraTable1['Run'].tolist()[i],sraTable1['body_site'].tolist()[i]
#     try: dicSRA[tissue].append(accID)
#     except KeyError:
#         dicSRA[tissue] = []
#         dicSRA[tissue].append(accID)

# for i in range(sraTable2.shape[0]):
#     accID, tissue= sraTable2['Run'].tolist()[i],sraTable2['body_site'].tolist()[i]
#     try: dicSRA[tissue].append(accID)
#     except KeyError:
#         dicSRA[tissue] = []
#         dicSRA[tissue].append(accID)
    
# for i in range(sraTable3.shape[0]):
#     accID, tissue= sraTable3['Run'].tolist()[i],sraTable3['body_site'].tolist()[i]
#     try: dicSRA[tissue].append(accID)
#     except KeyError:
#         dicSRA[tissue] = []
#         dicSRA[tissue].append(accID)  
        


# class Spliter():
#     def __init__(self,good):   # good is seperable based on some index
#         self.good = good
#         self.queue = []  # to store splited data
    
#     def split_df(self,n):
#         dim  = self.good.shape[0]
#         lis = [i for i in range(dim)]
#         size = len(lis)//n + 1
#         for j in range(0,len(lis),size):
#             part_index = lis[j:j+size]
#             part_df = self.good.iloc[part_index]
            
#             self.queue.append(part_df)
    
# Spliter1 = Spliter(sraData) 
# Spliter1.split_df(20)  
# liaisonData = Spliter1.queue 
  
# import multiprocessing as mp

# class ParallelWorker4GeneratingGTExDictionary(object):
#     def __init__(self,n):  # n means how many process you specify
#         self.n = n
#         self.queue = mp.Queue()
#         self.processes = []
#         self.result=[]
#         self.lock = mp.Lock()
     
#     def constructDic(self,sraData,index):
#         dicTissueExp = {}
#         for i in range(sraData.shape[0]):
#             print('this is the {0}th run of process{1}'.format(i,index))
#             event = sraData['UID'].tolist()[i]
#             for tissue,accID in dicSRA.items():
#                 try: dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
#                 except KeyError:            
#                     dicTissueExp[event] = {}
#                     dicTissueExp[event][tissue] = sraData.iloc[i][[accID+'_1.bed' for accID in dicSRA[tissue]]].values
#         self.lock.acquire()
#         self.queue.put(dicTissueExp)
#         self.result.append(self.queue.get())
#         self.lock.release()

#     def run(self):
#         global liaisonData


#         for i in range(self.n):
#             p = mp.Process(target=self.constructDic, args=(liaisonData[i],i+1,))
#             self.processes.append(p)
#         for x in self.processes:
#             x.start()
#         for x in self.processes:
#             x.join() 
            



         
# multiWorker = ParallelWorker4GeneratingGTExDictionary(20)

# multiWorker.run()

        
    
# def merge_dicts(dic_list):

#     result = {}
#     for dictionary in dic_list:
#         result.update(dictionary)
#     return result

# dicTissueExp = merge_dicts(multiWorker.result)

# print(len(multiWorker.result))

# import pickle
# with open('./imme.p','wb') as file3:
#     pickle.dump(multiWorker.result,file3)
# print('succeed 1')

# import pickle
# with open('./dicTissueExp2.p','wb') as file1:
#     pickle.dump(dicTissueExp,file1)
# print('succeed 2')


# import bz2
# import _pickle as cpickle
# with bz2.BZ2File('dicTissueExp2.pbz2','wb') as f1:
#     cpickle.dump(dicTissueExp,f1)   
# print('succeed 3')


# TASK 5 try again with Pool method
import pandas as pd
import os

sraTable1 = pd.read_csv('./SraRunTable-GTEX1.txt',sep=',')
sraTable2 = pd.read_csv('./SraRunTable-GTEX2.txt',sep=',')
sraTable3 = pd.read_csv('./SraRunTable-GTEX3.txt',sep=',')
sraData = pd.read_csv('./GTEx_EventAnnotation.txt',sep='\t')

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
Spliter1.split_df(20)  
liaisonData = Spliter1.queue 

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



import multiprocessing as mp
p = mp.Pool(processes = 20)
dicts = p.map(constructDic,liaisonData)

def merge_dicts(dic_list):

    result = {}
    for dictionary in dic_list:
        result.update(dictionary)
    return result

dicTissueExp = merge_dicts(dicts)

print(len(dicts))

import pickle
with open('./imme.p','wb') as file3:
    pickle.dump(dicts,file3)
print('succeed 1')

import pickle
with open('./dicTissueExp2.p','wb') as file1:
    pickle.dump(dicTissueExp,file1)
print('succeed 2')


import bz2
import _pickle as cpickle
with bz2.BZ2File('dicTissueExp2.pbz2','wb') as f1:
    cpickle.dump(dicTissueExp,f1)   
print('succeed 3')





















