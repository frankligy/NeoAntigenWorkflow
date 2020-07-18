'''
Author: Frank Li <li2g2@mail.uc.edu>
July 10th 2020 10:36PM

'''

import pandas as pd
import numpy as np
import os
import ast



def extract(k,HLA,folder,taskName,MHC):
    start = k[0]  # assume it will be 8
    start_filePath = os.path.join(folder,'Neoantigen_{0}_{1}.txt'.format(start,taskName))
    start_data = pd.read_csv(start_filePath,sep='\t')
    start_colName = MHC + 'result'
    start_data_base = start_data.drop([start_colName],axis=1)

    for i in k:
        filePath = os.path.join(folder,'Neoantigen_{0}_{1}.txt'.format(i,taskName))
        data = pd.read_csv(filePath,sep='\t')
        colName = MHC + 'result'
        targets = data[colName]   # Series object
        start_data_base = start_data_base.join(pd.Series(targets,name='{0}result_{1}mer'.format(MHC,i)))
    
    start_data_base.to_csv(os.path.join(folder,'merged_result_{0}.txt'.format(MHC)),sep='\t',index=None)






if __name__ == '__main__':
    taskName = 'PSI.Normal-like_vs_Others'   # same as -t in mhcPresent.py
    folder_prefix = '/data/salomonis2/LabFiles/Frank-Li/Anu/first'   # same as -o in mhcPresent.py
    folder = os.path.join(folder_prefix, 'resultMHC_{0}'.format(taskName))
    k = [8,9,10,11]
    HLA = ['HLA-A02:01','HLA-A01:01','HLA-B07:02','HLA-B08:01','HLA-C07:01','HLA-C04:01']   # same as -H in mhcPresent.py
    MHC = 'MHCI'
    extract(k,HLA,folder,taskName,MHC)