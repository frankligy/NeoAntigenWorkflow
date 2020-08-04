import pandas as pd
import numpy as np
import os
from collections import Counter




def read_hla_fasta(path):
    dic = {}
    with open(path,'r') as f:
        last = None
        for line in f:
            line=line.rstrip('\n')
            if line.startswith('>'):
                dic[line[1:]] = []
                last = line[1:]
            else:
                dic[last].extend(line.split(' '))
    return dic


if __name__ == '__main__':
    ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/OptiType/PAM50ClassifiedTCGAPatients.txt',sep='\t')
    hla = read_hla_fasta('/data/salomonis2/LabFiles/Frank-Li/OptiType/result.txt')


    basal = []
    luminalA = []
    luminalB = []
    her2 = []
    normal = []
    na = []

    col = []
    for i in range(ori.shape[0]):
        id_ = ori.iloc[i]['Complete TCGA ID']
        type_ = ori.iloc[i]['PAM50 mRNA']

        col.append(hla[id_])




        if type_ == 'Basal-like':
            basal.extend(hla[id_])
        elif type_ == 'HER2-enriched':
            her2.extend(hla[id_])
        elif type_ == 'Luminal A':
            luminalA.extend(hla[id_])
        elif type_ == 'Luminal B':
            luminalB.extend(hla[id_])
        elif type_ == 'Normal-like':
            normal.extend(hla[id_])
        else:
            na.extend(hla[id_])

    new = ori.join(pd.Series(col,name='hla genotypes'))
    new.to_csv('/data/salomonis2/LabFiles/Frank-Li/OptiType/with_hla_types.txt',sep='\t',index=None)

    with open('counts_hla_type.txt','w') as f:
        print('basal subtype:\n',Counter(basal),file=f)
        print('luminalA subtype:\n',Counter(luminalA),file=f)
        print('luminalB subtype:\n',Counter(luminalB),file=f)
        print('normal-like subtype:\n',Counter(normal),file=f)
        print('na subtype:\n',Counter(na),file=f)
    
        

