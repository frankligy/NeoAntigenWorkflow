#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 20:39:42 2020

@author: ligk2e
"""





import pandas as pd
import os
import numpy as np
os.chdir('/Users/ligk2e/Desktop/project_LUAD')



df_event = pd.read_csv('LUAD_Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-filtered-75p.txt',sep='\t')
df_event = df_event.drop('UID',axis=1)
df_event = df_event.rename(columns={'UID.1':'UID'})




df_groups = pd.read_csv('groups.txt',sep='\t',header=None)

dic = {}
for i in range(df_groups.shape[0]):
    id_, label = df_groups[0].tolist()[i],df_groups[2].tolist()[i]
    if label == 'R1-V7':
        try: dic['R1-V7'].append(id_)
        except KeyError:
            dic['R1-V7'] = []
            dic['R1-V7'].append(id_)
    elif label == 'Healthy':
        try: dic['Healthy'].append(id_)
        except KeyError:
            dic['Healthy'] = []
            dic['Healthy'].append(id_)

col = []
for i in range(df_event.shape[0]):
    mean_h = np.mean(df_event.iloc[i][dic['Healthy']].tolist())
    cond = True if mean_h==0 else False
    col.append(cond)
df_guide = df_event[pd.Series(col)]
df_guide.to_csv('feature_EventAnnotation.txt',sep='\t',index=None)
    
    
    









dicTissueExp['TSPAN6:ENSG00000000003:I3.1-E4.1|ENSG00000000003:E3.4-E4.1']['Cells - Cultured fibroblasts']
inspectGTEx('KYAT3:ENSG00000137944:E4.1-I4.1_88965413|ENSG00000137944:E4.1-E5.1')
inspectGTEx('NSUN5P2:ENSG00000106133:I2.1-E3.1|ENSG00000106133:E2.1-E3.1')   # inconsistent between GTEx and TCGA

tissue = list(dicSRA.keys())

def inspectGTEx(event,tissue='all',plot=True):
    flag = 0
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    import numpy as np
    global dicTissueExp
    if tissue=='all':
        tissueExp = dicTissueExp[event]
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
                    plt.savefig('./figures/{1}.pdf'.format(event,tis),bbox_inches='tight')
                    plt.close(fig)
                else: continue
            else: 
                flag += 1
                print('No expression in {}'.format(tis))
            
    else:
        expression = dicTissueExp[event][tissue]
        exp = expression.astype('float64')
        exp=exp[np.logical_not(np.isnan(exp))]
        if exp.size == 0: print('{0} data incomplete'.format(tissue))
        elif np.any(exp):   # have non-zero element
            plt.bar(np.arange(len(exp)),exp,width=0.2,label='tissue')
            plt.legend()
            plt.savefig('./{}.pdf'.format(tissue),bbox_inches='tight')
            plt.show()
            print(expression)
    return flag  
        
"""        
['Artery - Tibial',
 'Adipose - Subcutaneous',
 'Breast - Mammary Tissue',
 'Heart - Atrial Appendage',
 'Uterus',
 'Thyroid',
 'Adrenal Gland',
 'Lung',
 'Spleen',
 'Brain - Cortex',
 'Adipose - Visceral (Omentum)',
 'Colon - Sigmoid',
 'Skin - Not Sun Exposed (Suprapubic)',
 'Small Intestine - Terminal Ileum',
 'Whole Blood',
 'Esophagus - Gastroesophageal Junction',
 'Muscle - Skeletal',
 'Cells - Cultured fibroblasts',
 'Liver',
 'Minor Salivary Gland',
 'Skin - Sun Exposed (Lower leg)',
 'Stomach',
 'Ovary',
 'Heart - Left Ventricle',
 'Esophagus - Muscularis',
 'Nerve - Tibial',
 'Colon - Transverse',
 'Brain - Hypothalamus',
 'Brain - Nucleus accumbens (basal ganglia)',
 'Artery - Aorta',
 'Pituitary',
 'Testis',
 'Artery - Coronary',
 'Pancreas',
 'Brain - Cerebellum',
 'Kidney - Cortex',
 'Vagina',
 'Esophagus - Mucosa',
 'Prostate',
 'Brain - Caudate (basal ganglia)',
 'Brain - Putamen (basal ganglia)',
 'Cells - EBV-transformed lymphocytes',
 'Brain - Cerebellar Hemisphere',
 'Brain - Frontal Cortex (BA9)',
 'Brain - Spinal cord (cervical c-1)',
 'Brain - Amygdala',
 'Brain - Substantia nigra',
 'Brain - Hippocampus',
 'Brain - Anterior cingulate cortex (BA24)']        
"""        
        
        
import pandas as pd
import ast
pep = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/peptide_no_digestion.txt',header=None)[0].tolist()
idx8,idx9,idx10,idx11=[],[],[],[]
for item in pep:
    length = len(item)
    if length==8:
        mer8_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_8_mark.txt',sep='\t')
        mer8 = mer8_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer8):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-A29:02']
                if [item] in candidate: idx8.append(idx)

    if length==9:
        mer9_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_9_mark.txt',sep='\t')
        mer9 = mer9_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer9):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-A29:02']
                if [item] in candidate: idx9.append(idx)
       
    if length==10:
        mer10_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_10_mark.txt',sep='\t')
        mer10 = mer10_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer10):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-A29:02']
                if [item] in candidate: idx10.append(idx)
   
    if length==11:
        mer11_df = pd.read_csv('/Users/ligk2e/Desktop/project_test/resultMHC/NeoJunction_11_mark.txt',sep='\t')
        mer11 = mer11_df['MHCIresult'].tolist()
        for idx,candidate in enumerate(mer11):
            if not candidate == 'No candidates':
                candidate = ast.literal_eval(candidate)
                candidate = candidate['HLA-A29:02']
                if [item] in candidate: idx11.append(idx)

        
new8_df = mer8_df.iloc[idx8]
new9_df = mer9_df.iloc[idx9] 
new10_df = mer10_df.iloc[idx10] 
new11_df = mer11_df.iloc[idx11] 
new = pd.concat([new8_df,new9_df,new10_df,new11_df])



    
        
        
        
        
        
        
        
            
            