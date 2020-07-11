'''
Author: Frank Li <li2g2@mail.uc.edu>
July 10th 2020 10:36PM

'''

import pandas as pd
import numpy as np
import os
import ast

def extract(k,HLA,folder,taskName,MHC):
    with open(os.path.join(folder,'augment_{0}.fasta'.format(taskName)),'w') as f:
        for i in k:
            filePath = os.path.join(folder,'Neoantigen_{0}_{1}.txt'.format(i,taskName))
            data = pd.read_csv(filePath,sep='\t')
            colName = MHC + 'result'
            targets = data[colName]   # Series object
            for target in targets:
                if target == 'No candidates': continue
                else:
                    target = ast.literal_eval(target)   # a dict
                    for hla in HLA:
                        strongBinding = target[hla][0]
                        for sb in strongBinding:
                            f.write('>{0}mer_{1}_strongBinding\n'.format(i,hla))
                            f.write(sb+'\n')
                        weakBinding = target[hla][1]
                        for wb in weakBinding:
                            f.write('>{0}mer_{1}_weakBinding\n'.format(i,hla))
                            f.write(wb+'\n')




if __name__ == '__main__':
    taskName = 'TCGA-E2-A10A-01'
    folder_prefix = '/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast'
    folder = os.path.join(folder_prefix, 'resultMHC_{0}'.format(taskName))
    k = [8,9,10,11]
    HLA = ['HLA-A02:06','HLA-A01:01','HLA-B35:01','HLA-B35:03','HLA-C04:04','HLA-C04:01']
    MHC = 'MHCI'
    extract(k,HLA,folder,taskName,MHC)
