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
                    if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                        continue
                    else:

                        exp = np.array(exp)
            
                        exp = exp.astype('float64')
                        exp[np.isnan(exp)] = 0.0   # nan means can not detect the gene expression
                        hits = sum([True if i > cutoff_PSI else False for i in exp])   # in a tissue, how many samples have PSI > cutoff value
                        total = exp.size    # how many samples for each tissue type
                        sampleRatio = hits/total    # percentage of sampels that are expressing this event
                        if sampleRatio > cutoff_sampleRatio: tissueCounter += 1   # this tissue is expressiing this event
                tissueRatio = tissueCounter/51    # 51 tissue types in total,excluding 3 cancer cell lines
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

    



def check_GTEx_PSI_pbz2(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio,taskName):
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
        dicTissueExp = cpickle.load(f1)  
    col = []
    for i in range(df.shape[0]):
        #print(i)
        UID = df.iloc[i]['UID']
        event = UID.split('|')[0]       # 0 means foreground event, 1 means background event
        try:
            tissueExp = dicTissueExp[event]  # {heart:[],brain:[]}   # values are ndarray
        except KeyError:
            cond = 'Unknown'
            col.append(cond)
        else:
            tissueCounter = 0
            for tis,exp in tissueExp.items():

                if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                    continue
                else:
                    exp = exp.astype('float64')
                    exp[np.isnan(exp)] = 0.0   # nan means can not detect the gene expression
                    hits = sum([True if i > cutoff_PSI else False for i in exp])   # in a tissue, how many samples have PSI > cutoff value
                    total = exp.size    # how many samples for each tissue type
                    sampleRatio = hits/total    # percentage of sampels that are expressing this event
                    if sampleRatio > cutoff_sampleRatio: tissueCounter += 1   # this tissue is expressiing this event
            tissueRatio = tissueCounter/51    # 51 tissue types in total, excluding 3 cancer cell lines
            if tissueRatio > cutoff_tissueRatio:
                cond = 'False'
                col.append(cond)
            else: 
                cond = 'True'
                col.append(cond)
    df['check'] = col
    df_true = df[df['check']=='True']
    df_unknown = df[df['check']=='Unknown']
    df_false = df[df['check']=='False']
    #df_True.to_csv(os.path.join(outFolder,'df_True.txt'),sep='\t',index=None)
    #df_Unknown.to_csv(os.path.join(outFolder,'df_Unknown.txt'),sep='\t',index=None)
    #df_False.to_csv(os.path.join(outFolder,'df_False.txt'),sep='\t',index=None)

    df_train = df.loc[df['check']!='Unknown']
    df_train = df_train.set_index(pd.Index(np.arange(df_train.shape[0])))
    df_train.to_csv(os.path.join(outFolder,'df_train_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_true = df.loc[df['check'] == 'True']
    df_true = df_true.set_index(pd.Index(np.arange(df_true.shape[0])))
    df_true.to_csv(os.path.join(outFolder,'df_true_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_unknown = df.loc[df['check'] == 'Unknown']
    df_unknown = df_unknown.set_index(pd.Index(np.arange(df_unknown.shape[0])))
    df_unknown.to_csv(os.path.join(outFolder,'df_unknown_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_false = df.loc[df['check'] == 'True']
    df_false = df_false.set_index(pd.Index(np.arange(df_false.shape[0])))
    df_false.to_csv(os.path.join(outFolder,'df_false_{0}.txt'.format(taskName)),sep='\t',index=None)   

    df.to_csv(os.path.join(outFolder,'df_whole_{0}.txt'.format(taskName)),sep='\t',index=None)


if __name__ == '__main__':

    dataFolder = '/data/salomonis2/LabFiles/Frank-Li/python3/data'
    outFolder = '/data/salomonis2/LabFiles/Frank-Li/specificity-score/now'
    taskName = 'breast'
    cutoff_PSI = 0.1
    cutoff_sampleRatio = 0.1
    cutoff_tissueRatio = 0.1
    df = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/specificity-score/before/TCGA-E2-A10A-01.ori.txt',sep='\t')
    check_GTEx_PSI_pbz2(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio,taskName)