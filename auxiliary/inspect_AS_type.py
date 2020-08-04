#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 18:46:10 2020

@author: ligk2e
"""

import os
import pandas as pd
import numpy as np
import argparse
import ast



def main(args):
    peptide = args.peptide
    merge = args.merge
    outdir = args.outdir
    task = args.task

    df_peptide = pd.read_csv(peptide,sep='\t')
    df_merge = pd.read_csv(merge,sep='\t')

    # process df_merge
    col = []
    for row in df_merge.itertuples():
        bucket = []
        mer8 = row.MHCIresult_8mer
        mer9 = row.MHCIresult_9mer
        mer10 = row.MHCIresult_10mer
        mer11 = row.MHCIresult_11mer
        for item in [mer8,mer9,mer10,mer11]:
            if item == 'No candidates': continue
            else:
                item = ast.literal_eval(item)
                collect = list(item.values())   # [([],[]),([],[]),([],[])..]
                for each in collect:
                    bucket.extend(each[0])
                    bucket.extend(each[1])
        bucket = list(set(bucket))
        col.append(bucket)
    df_merge=df_merge.join(pd.Series(col,name='collection'))

    # start to match
    ms = set(df_peptide.iloc[:,0].tolist())
    cond = []
    for row in df_merge.itertuples():
        lis = row.collection   # all the neoantigen predicted for each event
        if lis:   # if not empty
            comm = set(lis).intersection(ms)   # compare if two list have common item, use set.intersection
            if comm:
                cond.append(True)
            else:
                cond.append(False)
        else:
            cond.append(False)
    #print(len(cond),df_merge.shape[0])
    df_merge = df_merge.loc[pd.Series(cond)]
    df_merge = df_merge.set_index(pd.Index(np.arange(df_merge.shape[0])))

    df_merge.to_csv(os.path.join(outdir,'{0}.AS_type.txt'.format(task)),sep='\t',index=None)
                    
            

        
    






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='inspect the AS types each MS-supported Neoangigens belong to')
    parser.add_argument('--peptide',type=str,default='.',help='path to the peptide list we identified in MaxQuant')
    parser.add_argument('--merge',type=str,default='.',help='path to the merge file')
    parser.add_argument('--outdir',type=str,default='.',help='path to the output folder')
    parser.add_argument('--task',type=str,default='.',help='give your task a name')
    args = parser.parse_args()
    main(args)
    



