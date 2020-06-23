import os
import pandas as pd
import numpy as np
import h5py
import _pickle as cpickle
import bz2






def check_GTEx_PSI_hdf5(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio):
    col = []
    for i in range(df.shape[0]):
        print(i)
        UID = df.iloc[i]['UID']
        event = UID.split('|')[0]       # foreground event
        with h5py.File(os.path.join(dataFolder,'dicTissueExp.hdf5'),'r') as f:
            try:                
                tissueExp = f[event]
            except:
                cond = 'Unknown'
                col.append(cond)
            else:
                tissueCounter = 0
                for tis,exp in tissueExp.items():
                    exp = np.array(exp)
        
                    exp = exp.astype('float64')
                    exp[np.isnan(exp)] = 0.0   # nan means can not detect the gene expression
                    hits = sum([True if i > cutoff_PSI else False for i in exp])   # in a tissue, how many samples have PSI > cutoff value
                    total = exp.size    # how many samples for each tissue type
                    sampleRatio = hits/total    # percentage of sampels that are expressing this event
                    if sampleRatio > cutoff_sampleRatio: tissueCounter += 1   # this tissue is expressiing this event
                tissueRatio = tissueCounter/54    # 54 tissue types in total
                if tissueRatio > cutoff_tissueRatio:
                    cond = 'False'
                    col.append(cond)
                else: 
                    cond = 'True'
                    col.append(cond)
    df['cond'] = col
    df_True = df[df['cond']=='True']
    df_Unknown = df[df['cond']=='Unknown']
    df_False = df[df['cond']=='False']
    df_True.to_csv(os.path.join(outFolder,'df_True.txt'),sep='\t',index=None)
    df_Unknown.to_csv(os.path.join(outFolder,'df_Unknown.txt'),sep='\t',index=None)
    df_False.to_csv(os.path.join(outFolder,'df_False.txt'),sep='\t',index=None)



def check_GTEx_PSI_pbz2(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio):
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
        dicTissueExp = cpickle.load(f1)  
    col = []
    for i in range(df.shape[0]):
        print(i)
        UID = df.iloc[i]['UID']
        event = UID.split('|')[0]       # foreground event
        try:
            tissueExp = dicTissueExp[event]  # {heart:[],brain:[]}   # values are ndarray
        except KeyError:
            cond = 'Unknown'
            col.append(cond)
        else:
            tissueCounter = 0
            for tis,exp in tissueExp.items():
    
                exp = exp.astype('float64')
                exp[np.isnan(exp)] = 0.0   # nan means can not detect the gene expression
                hits = sum([True if i > cutoff_PSI else False for i in exp])   # in a tissue, how many samples have PSI > cutoff value
                total = exp.size    # how many samples for each tissue type
                sampleRatio = hits/total    # percentage of sampels that are expressing this event
                if sampleRatio > cutoff_sampleRatio: tissueCounter += 1   # this tissue is expressiing this event
            tissueRatio = tissueCounter/54    # 54 tissue types in total
            if tissueRatio > cutoff_tissueRatio:
                cond = 'False'
                col.append(cond)
            else: 
                cond = 'True'
                col.append(cond)
    df['cond'] = col
    df_True = df[df['cond']=='True']
    df_Unknown = df[df['cond']=='Unknown']
    df_False = df[df['cond']=='False']
    #df_True.to_csv(os.path.join(outFolder,'df_True.txt'),sep='\t',index=None)
    #df_Unknown.to_csv(os.path.join(outFolder,'df_Unknown.txt'),sep='\t',index=None)
    #df_False.to_csv(os.path.join(outFolder,'df_False.txt'),sep='\t',index=None)

    df_train = df[df['cond']!='Unknown']
    df_train.to_csv(os.path.join(outFolder,'df_train_LUAD.txt'),sep='\t',index=None)


dataFolder = '/data/salomonis2/LabFiles/Frank-Li/python3/data'
outFolder = '/data/salomonis2/LabFiles/Frank-Li/specificity-score'
cutoff_PSI = 0.1
cutoff_sampleRatio = 0.1
cutoff_tissueRatio = 0.1
df = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/specificity-score/test.txt',sep=' ')
check_GTEx_PSI_pbz2(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio)