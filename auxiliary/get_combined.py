'''
Author: Frank Li <li2g2@mail.uc.edu>
July 10th 2020 10:36PM

'''

import pandas as pd
import numpy as np
import os
import ast
import argparse



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



def change_HLA_format(HLA):   #'HLA-A32:01,HLA-A68:01,HLA-B40:01,HLA-B51:01,HLA-C02:02,HLA-C16:01'
    result = HLA.split(',')
    return result



def main(args):

    taskName = args.task
    folder_prefix = args.outdir
    HLA = change_HLA_format(args.HLA)


    folder = os.path.join(folder_prefix, 'resultMHC_{0}'.format(taskName))
    k = [8,9,10,11]
    HLA = ['HLA-A02:06','HLA-A01:01','HLA-B35:01','HLA-B35:03','HLA-C04:04','HLA-C04:01']   # same as -H in mhcPresent.py
    MHC = 'MHCI'
    extract(k,HLA,folder,taskName,MHC)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='combine 8mer, 9mer, 10mer, 11mer')
    parser.add_argument('--task',type=str,default=None,help='task name, same as task in mhcPresent.py')
    parser.add_argument('--outdir',type=str,default='.',help='output dicrectory, same as outFolder in mhcPresent.py')
    parser.add_argument('--HLA',type=str,default=None,help='HLA type we want to inspect, same as HLA in mhcPresent.py')
    args=parser.parse_args()
    main(args)

