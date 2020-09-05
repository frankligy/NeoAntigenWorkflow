'''
this script is for analyzing TCGA Breast cancer patients AS-derived neoantigen burden
'''


'''
cat counts.TCGA-BRCA.txt | head -n 1| awk '{for(i=0;i<=NF;i++){print $i}}' | tail -n +2 > all_counts_column_id.txt
grep "TCGA-[A-Za-z0-9]*-[A-Za-z0-9]*-11" all_counts_column_id.txt > all_matched_control.txt

'''

import pandas as pd
import numpy as np
import _pickle as cpickle
import bz2
import os
import ast
import argparse
import shelve
import statistics
from sklearn import preprocessing
import pickle
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis,QuadraticDiscriminantAnalysis
from joblib import dump,load


def convert_nan_to_zero(col2,tcga_tumor):  # col2 is a Series
    array = col2.values   # ndarray
    array[np.isnan(array)] = 0.0
    col2 = pd.Series(array,name=tcga_tumor)
    return col2     # col2 is a series


def extract_foreground_event(col1):
    '''
    AC233968.1:ENSG00000274615:E12.1_38124630-ENSG00000278599:E2.1|ENSG00000274615:E12.1_38124630-ENSG00000278599:E3.1
    AC233968.1:ENSG00000274615:E5.1_38117740-E6.1|ENSG00000274615:E5.1-E6.1

    part will be : 
    ENSG00000274615:E12.1_38124630-ENSG00000278599:E2.1
    ENSG00000274615:E5.1_38117740-E6.1
    '''
    part = []
    UID = col1.tolist()
    for i in range(len(UID)):   
        uid = UID[i]
        x = uid.split('|')    # seperate foreground and background
        try: x[0].split(':')[3]         # if index 3 exist, foreground event will be a trans-splicing event
        except IndexError: event = x[0].split(':')[2]   # index 3 doesn't exsit, foreground won't be a trans-splicing event
        else: event = str(x[0].split(':')[2])+':'+str(x[0].split(':')[3])
        query = uid.split('|')[1].split(':')[0] + ':' + event   # ensg in background part + ':' + foreground event
        part.append(query)
    return part

def find_PSI_col(tcga,events):
    # tcga: TCGA-A1-A0SB (Patient ID)
    # In events, column either TCGA-A1-A0SB-01(tumor sample) or TCGA-A1-A0SB-11 (matched control)
    tcga_tumor = tcga + '-01'
    col1 = events['UID']
    part = extract_foreground_event(col1)
    col2 = events[tcga_tumor]
    col2 = convert_nan_to_zero(col2,tcga_tumor)
    df_start = pd.concat([col1,col2],axis=1)
    '''
            UID    TCGA-A1-A0SB-01
    0       ...       ...
    1       ...       ...

    '''
    return df_start,part


def convert_controls_to_list(controls):
    return controls[0].tolist()


def ave_counts_control(counts_control):
    mat = counts_control.iloc[:,1:].to_numpy()
    mat[np.isnan(mat)] = 0
    row_mean = np.mean(mat,axis=1)   # sames like np.mean, row-wise must be axis=1
    col1 = counts_control['AltAnalyze_ID']
    col2 = pd.Series(row_mean,name='mean')
    return pd.concat([col1,col2],axis=1)


def construct_count_dic(counts_id_useful,query_id_counts,mean_counts):   # these three are all list
    dic = {}
    for i in range(len(counts_id_useful)):
        ensg = counts_id_useful[i].split(':')[0]
        event = counts_id_useful[i].split(':')[1:]   
        event = ''.join(event)
        info = query_id_counts[i]
        mean_ = mean_counts[i]
        try:
            dic[ensg].append((event,info,mean_))
        except KeyError:
            dic[ensg] = []
            dic[ensg].append((event,info,mean_))
    return dic

def dissect_counts_extract(counts_extract):
    '''
    counts_extract:
    [
        (counts,mean_)
        (counts,mean_)
    ]

    '''

    a = pd.Series([item[0] for item in counts_extract],name='counts')
    b = pd.Series([item[1] for item in counts_extract],name='mean_in_control')

    return a,b


def match_countsBYdic(part,dic):
    counts_extract = []
    for i in range(len(part)):
        flag = False
        part_ensg = part[i].split(':')[0]
        searchSpace = dic[part_ensg]
        for j in searchSpace:
            match = part_ensg+':'+j[0]
            if part[i] == match:
                counts_extract.append(j[1:])
                flag = True
                break
        if flag == False:
            counts_extract.append(('no matched counts','no matched counts'))
    
    a,b = dissect_counts_extract(counts_extract)
    return a,b
    


def add_counts(tcga,df_start,part,counts,mean_df):
    '''
    ENSG00000000419:E1.4-E2.1_50955256=chr20:50958363-50955256
    '''
    tcga = tcga + '-01'
    counts_id = counts['AltAnalyze_ID'].tolist()
    counts_id_useful = [item.split('=')[0] for item in counts_id]   # remove chr coordination, retain ENSG

    new = [name[:15] for name in counts.columns]   # TCGA-AR-A2LJ-01A-12R-A19W-07.bed to TCGA-AR-A2LJ-01
    counts.columns = new   # change all column names to truncated version
    query_id_counts = counts[tcga].tolist()

    mean_counts = mean_df['mean'].tolist()

    dic = construct_count_dic(counts_id_useful,query_id_counts,mean_counts)
    a,b = match_countsBYdic(part,dic)
    #print(a,b,len(a),len(b))
    df_new = df_start.join(a).join(b)
    return df_new

def filter_PSI(df,tcga):
    psi = df[tcga]
    df = df.loc[psi!=0.0,:]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    return df

def filter_counts(df):
    counts = df['counts']
    df = df.loc[counts!='no matched counts',:]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    return df

def filter_mean(df):    # control counts should be <= 3
    mean_ = df['mean_in_control']
    df = df.loc[mean_ <= 3.0,:]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    return df

def filter_further(df):   # tumor counts should be > 3
    col = df['counts']
    df = df.loc[col > 3,:]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    return df

def filter_effect_size(df):    # tumor counts / mean_in_control should be >= 3
    col1 = df['counts'].values   # ndarray
    col2 = df['mean_in_control'].values    # ndarray
    col2 = [0.0001 if item== 0 else item for item in col2]
    divide = col1 / col2 
    #divide[np.isinf(divide)] = 10000
    df = df.loc[pd.Series(divide) > 3, :]
    df = df.set_index(pd.Index(np.arange(df.shape[0])))
    return df



def filtering(df_new,tcga):
    tcga = tcga + '-01'
    df_new = filter_PSI(df_new,tcga)
    df_new = filter_counts(df_new)
    df_new = filter_mean(df_new)
    df_new = filter_further(df_new)
    df_new = filter_effect_size(df_new)
    return df_new




def convert_hla(HLA):
    HLA = HLA.split(',')
    new = []
    for item in HLA:
        item = item.replace('*','')
        item = 'HLA-' + item
        new.append(item)
    return new


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


    df_train = df.loc[df['check']!='Unknown']
    df_train = df_train.set_index(pd.Index(np.arange(df_train.shape[0])))
    df_train.to_csv(os.path.join(outFolder,'df_train_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_final = df.loc[df['check']!='False']
    df_final = df_final.set_index(pd.Index(np.arange(df_final.shape[0])))
    df_final.to_csv(os.path.join(outFolder,'df_final_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_true = df.loc[df['check'] == 'True']
    df_true = df_true.set_index(pd.Index(np.arange(df_true.shape[0])))
    df_true.to_csv(os.path.join(outFolder,'df_true_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_unknown = df.loc[df['check'] == 'Unknown']
    df_unknown = df_unknown.set_index(pd.Index(np.arange(df_unknown.shape[0])))
    df_unknown.to_csv(os.path.join(outFolder,'df_unknown_{0}.txt'.format(taskName)),sep='\t',index=None)

    df_false = df.loc[df['check'] == 'False']
    df_false = df_false.set_index(pd.Index(np.arange(df_false.shape[0])))
    df_false.to_csv(os.path.join(outFolder,'df_false_{0}.txt'.format(taskName)),sep='\t',index=None)   

    df.to_csv(os.path.join(outFolder,'df_whole_{0}.txt'.format(taskName)),sep='\t',index=None)
    return df_final

def PSI(df,cutoff):
    ave_PSI, median_PSI, percentage_PSI = [],[],[]
    print('loading dicTissueExp_psi file')
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp.pbz2'),'rb') as f1:
        dicTissueExp_psi = cpickle.load(f1)
    print('successfully load dicTissueExp_psi file into RAM')
    for i in range(df.shape[0]):
        cumulative = 0   # cumulative PSI values, with the intent of calculating average PSI
        lis = []        # record all PSI value for calcating median
        hits = 0        # how many samples are above threshold
        n = 0             # total amounts of samples
        event = df.iloc[i]['UID'].split('|')[0]
        tissueDic = dicTissueExp_psi[event]
        for tis,exp in tissueDic.items():
            exp = np.array(exp)   # turn to ndarray, and dtype is object
            exp = exp.astype('float64')
            exp[np.isnan(exp)] = 0.0   # because nan just means they don't even have expression
            if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                continue
            else: 
                if np.any(exp):   # exist expresson in this tissue type
                    cumulative += sum(exp)   # add all PSI value in this tissue to cumulative
                    lis.extend(exp.tolist())       # inject all samples' PSI to lis
                    hits += sum([True if item>cutoff else False for item in exp])   # calculate how many hits and add it to hits
                    n += exp.size           # calculate how many tissues and add to n
                else:    # all zero
                    lis.extend(exp.tolist())
                    n += exp.size
        # calculate three metrics for each event
        ave = cumulative / n
        median = statistics.median(lis)
        percentage = hits / n 
        # append them to final lists
        ave_PSI.append(ave)
        median_PSI.append(median)
        percentage_PSI.append(percentage)
    return ave_PSI,median_PSI,percentage_PSI

def ReadCounts(df,cutoff):
    ave_counts, median_counts, percentage_counts = [],[],[]
    print('loading dicTissueExp_counts file')
    with bz2.BZ2File(os.path.join(dataFolder,'dicTissueExp_counts.pbz2'),'rb') as f1:
        dicTissueExp_counts = cpickle.load(f1)
    print('successfully load dicTissueExp_counts file into RAM')
    for i in range(df.shape[0]):
        cumulative = 0   # cumulative read counts values, with the intent of calculating average PSI
        lis = []        # record all read counts value for calcating median
        hits = 0        # how many samples are above threshold
        n = 0             # total amounts of samples
        event = df.iloc[i]['UID'].split('|')[0]
        event = ':'.join(event.split(':')[1:])   # trim out the gene symbol, only keep ENSG:E3.4-E5.6
        tissueDic = dicTissueExp_counts[event]
        for tis,exp in tissueDic.items():
            exp = np.array(exp)   # turn to ndarray, and dtype is object
            exp = exp.astype('int')
            if tis == 'Cells - Cultured fibroblasts' or tis == 'Cells - Leukemia cell line (CML)' or tis == 'Cells - EBV-transformed lymphocyte':  # these tissue are tumor tissue, should be excluded
                continue
            else: 
                if np.any(exp):   # exist expresson in this tissue type
                    cumulative += sum(exp)   # add all PSI value in this tissue to cumulative
                    lis.extend(exp.tolist())       # inject all samples' PSI to lis
                    hits += sum([True if item>cutoff else False for item in exp])   # calculate how many hits and add it to hits
                    n += exp.size           # calculate how many tissues and add to n
                else:    # all zero
                    lis.extend(exp.tolist())
                    n += exp.size
        # calculate three metrics for each event
        ave = cumulative / n
        median = statistics.median(lis)
        percentage = hits / n 
        # append them to final lists
        ave_counts.append(ave)
        median_counts.append(median)
        percentage_counts.append(percentage)
    return ave_counts,median_counts,percentage_counts



def scoring_process(df_evaluation,cutoff_PSI,cutoff_readcounts,max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector):
    ave_PSI,median_PSI,percentage_PSI = PSI(df_evaluation,cutoff_PSI)
    ave_counts,median_counts,percentage_counts = ReadCounts(df_evaluation,cutoff_readcounts)
    
    mat_eval = np.column_stack((ave_PSI,median_PSI,percentage_PSI,ave_counts,median_counts,percentage_counts))
    print('shape of mat_eval:',mat_eval.shape)
    print(mat_eval)
    mat_eval_new = np.zeros(mat_eval.shape)

    for j in range(mat_eval.shape[1]):
        for i in range(mat_eval.shape[0]):
            new_ij = core_function(max_mat_ori[j],min_mat_ori[j],max_mat_rescale[j],min_mat_rescale[j],mat_eval[i,j])
            mat_eval_new[i,j] = new_ij
    print('shape of mat_eval_new:',mat_eval_new.shape)
    print(mat_eval_new)
    IWscore = []
    for m in range(df_evaluation.shape[0]):
        score = core_IW(mat_eval_new[m,:],leading_eigenvector)
        inverse_score = (-1.0) * float(score)
        sigmoid_score = sigmoid(inverse_score)
        IWscore.append(sigmoid_score)
    df_evaluation['IWscore'] = IWscore
    return df_evaluation



def wrapper_scoring_process(scoreFile,dataFolder):
    scoring = scoreFile
    # max_mat_ori = [9.98564855e-01,1.00000000e+00,1.00000000e+00,1.24263241e+04,1.02370000e+04,1.00000000e+00]
    # min_mat_ori = [0.00092801,0.0,0.00326531,0.01387755,0.0,0.0]
    # max_mat_rescale = [1.6203541,1.246249,0.98267483,72.70393268,80.23512846,0.77875848]
    # min_mat_rescale = [-1.69933277,-1.43360037,-2.65506607,-0.30237015,-0.30285234,-2.89327914]
    # leading_eigenvector = [0.48264742,0.47174347,0.47383551,0.21692184,0.22945297,0.46934607]
    s = shelve.open(os.path.join(dataFolder,'training_parameters_IW'))
    max_mat_ori = s['max_mat_ori']
    min_mat_ori = s['min_mat_ori']
    max_mat_rescale = s['max_mat_rescale']
    min_mat_rescale = s['min_mat_rescale']
    leading_eigenvector = s['leading_eigenvector']
    s.close()
    df_new = scoring_process(scoring,0.1,3,max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector)
    return df_new

def core_function(max_ori,min_ori,max_re,min_re,old):
    new = (max_re * (old - min_ori) + min_re * (max_ori - old)) / (max_ori - min_ori)
    return new

def core_IW(array,weight):
    IW = np.dot(array,weight)     # vectorization
    return IW

def sigmoid(x):
    x = float(x)
    y = 1/(1+np.exp(-x))
    return y

def HLA_seperator(hla):
    HLA1='HLA-A02:01,HLA-A01:01,HLA-C07:02,HLA-A03:01,HLA-C07:01,HLA-B07:02,HLA-A24:02,HLA-C03:04,HLA-C06:02,HLA-C04:01'     # netMHCpan4.1 format
    HLA2='HLA-C05:01,HLA-B44:02,HLA-B40:01,HLA-B08:01,HLA-A11:01,HLA-B58:01,HLA-B51:01,HLA-B15:01,HLA-A23:01,HLA-B35:01'
    HLA3='HLA-B18:01,HLA-C03:03,HLA-C01:02,HLA-B44:03,HLA-B13:02,HLA-C08:02,HLA-A68:01,HLA-C02:02,HLA-B14:02,HLA-A31:01'
    HLA4='HLA-C12:03,HLA-C15:02,HLA-C03:02,HLA-A26:01,HLA-B49:01,HLA-C16:01,HLA-B50:01,HLA-A33:03,HLA-C14:02,HLA-A33:01'
    HLA5='HLA-B35:03,HLA-A33:01,HLA-B15:17,HLA-C07:04,HLA-B57:01,HLA-B38:01,HLA-B53:01,HLA-B27:05,HLA-A32:01,HLA-B15:03'
    HLA6='HLA-B58:02,HLA-C02:10,HLA-A25:01,HLA-A02:05,HLA-B52:01,HLA-C12:02,HLA-B55:01,HLA-A68:02,HLA-A29:01,HLA-C16:02'
    HLA7='HLA-B39:06,HLA-B07:05,HLA-C15:05,HLA-A02:07,HLA-B46:01,HLA-A02:06,HLA-A24:23,HLA-B39:01,HLA-A74:01,HLA-B53:02'
    HLA8='HLA-A36:01,HLA-B35:02,HLA-B57:03,HLA-B27:13,HLA-A02:17,HLA-B35:08,HLA-B37:01,HLA-A66:01,HLA-B41:02,HLA-C17:01'
    HLA9='HLA-A69:01,HLA-B15:16,HLA-A29:02,HLA-B44:05,HLA-B27:09,HLA-B27:14,HLA-B27:02,HLA-B47:01,HLA-B15:18,HLA-B51:09'
    HLA10='HLA-A24:07,HLA-B40:06,HLA-B27:06,HLA-C01:03,HLA-B54:01,HLA-A34:02,HLA-B14:01,HLA-B40:02'
    HLA11='HLA-B08:09,HLA-B35:05,HLA-B27:03,HLA-C04:03,HLA-A02:11,HLA-B27:04,HLA-B38:02,HLA-A03:02,HLA-A24:03,HLA-B55:02,HLA-B39:24,HLA-B27:10,HLA-B78:01,HLA-A68:24'
    HLA12='HLA-A02:03,HLA-B56:01,HLA-A02:20,HLA-C16:04,HLA-B48:01,HLA-A66:02,HLA-B15:12,HLA-A02:02,HLA-B35:77,HLA-B51:08,HLA-C03:05,HLA-A30:01,HLA-B44:10,HLA-C03:17'
    HLA13='HLA-C15:04,HLA-B82:01,HLA-A30:04,HLA-B27:07,HLA-B18:03,HLA-C08:01,HLA-B56:04,HLA-B15:10,HLA-B41:01,HLA-B15:21,HLA-B42:01,HLA-B39:10,HLA-B45:01,HLA-B15:02'
    HLA14='HLA-A30:02,HLA-C18:01,HLA-C14:03,HLA-B13:01,HLA-A01:02,HLA-C04:04'
    array = [HLA1,HLA2,HLA3,HLA4,HLA5,HLA6,HLA7,HLA8,HLA9,HLA10,HLA11,HLA12,HLA13,HLA14]
    for i in range(len(array)):
        if hla in array[i]:
            index_ = i + 1  # No.index_ folder
            break
    return index_

def single(row,mer,hla):
    # mer is a Series, row is a int
    if mer[row] == 'No candidates': return 'No candidates'
    zoom = ast.literal_eval(mer[row])   # will be a dict
    value = zoom[hla]  # [[],[]]
    return {hla:value}





def from_merge_file(uid,hla,index_,row):
    merge = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/resultMHC_all_{0}/merged_result_all_{0}.txt'.format(str(index_)),sep='\t')
    if row is None:
        try:
            row = pd.Index(merge['UID']).get_loc(uid)
        except KeyError:
            return 'No candidates','No candidates','No candidates','No candidates',None   # this is a bug condition but haven't figured out
    mer8 = merge['HLAIpresented_dict_8mer']
    mer9 = merge['HLAIpresented_dict_9mer']
    mer10 = merge['HLAIpresented_dict_10mer']
    mer11 = merge['HLAIpresented_dict_11mer']
    #print(row,hla)
    dic8 = single(row,mer8,hla)
    dic9 = single(row,mer9,hla)
    dic10 = single(row,mer10,hla)
    dic11 = single(row,mer11,hla)
    #print(dic8,dic9,dic10,dic11)
    return dic8,dic9,dic10,dic11,row



def single_mer(DIC,dic):
    if DIC != 'No candidates':
        if dic != 'No candidates':
            try: 
                DIC.update(dic)
            except:
                DIC = dic
        else:
            DIC = 'No candidates'
    return DIC

def organize(DIC,HLA):
    for i in range(len(HLA)):
        if DIC[i] == 'No candidates':
            DIC[i] = {HLA[i]:[[],[]]}
    #print(DIC)
    base = DIC[0]
    #print(base)
    for i in range(1,6):
        base.update(DIC[i])
    #print(base)
    return base


def retrieve(df,HLA):
    col8,col9,col10,col11 = [],[],[],[]
    for i in range(df.shape[0]):
        print(i)
        uid = df['UID'].iloc[i]
        DIC8,DIC9,DIC10,DIC11 = [],[],[],[]
        row = None
        for hla in HLA:
            index_ = HLA_seperator(hla)
            dic8,dic9,dic10,dic11,row = from_merge_file(uid,hla,index_,row)
            DIC8.append(dic8)
            DIC9.append(dic9)
            DIC10.append(dic10)
            DIC11.append(dic11)

        DIC8 = organize(DIC8,HLA)
        DIC9 = organize(DIC9,HLA)
        DIC10 = organize(DIC10,HLA)
        DIC11 = organize(DIC11,HLA)

        col8.append(DIC8)
        col9.append(DIC9)
        col10.append(DIC10)
        col11.append(DIC11)
    df['HLAIpresented_dict_8mer'] = col8
    df['HLAIpresented_dict_9mer'] = col9
    df['HLAIpresented_dict_10mer'] = col10       
    df['HLAIpresented_dict_11mer'] = col11
    return df

def combine_all_mers(df,HLA):
    combine_col = []
    for i in range(df.shape[0]):
        combine = []
        col8 = df.iloc[i]['HLAIpresented_dict_8mer']
        col9 = df.iloc[i]['HLAIpresented_dict_9mer']
        col10 = df.iloc[i]['HLAIpresented_dict_10mer']
        col11 = df.iloc[i]['HLAIpresented_dict_11mer']
        for col in [col8,col9,col10,col11]:
            target = col
            #print(type(target))
            if target == 'No candidates': continue
            else:
                target = target  # a dict

                #print(HLA,target)
                for hla in HLA:
                    strongBinding = target[hla][0]
                    weakBinding = target[hla][1]
                    combine.extend(strongBinding)
                    combine.extend(weakBinding)
        combine = list(set(combine))
        combine_col.append(combine)
    df['neoantigens'] = combine_col
    return df

def count_all_mers(df):
    count_col = []
    for i in range(df.shape[0]):
        neoantigens = df['neoantigens'].iloc[i]
        count = len(neoantigens)
        count_col.append(count)
    df['neoantigen_counts'] = count_col
    return df

def construct_df(uid):
    return pd.DataFrame({'UID':uid},index=[0])
            
def run_specificity_score(df_final):
    dataFolder = '/data/salomonis2/LabFiles/Frank-Li/python3/data'
    col = []
    for i in range(df_final.shape[0]):
        uid = df_final.iloc[i]['UID']
        check = df_final.iloc[i]['check']
        if check == 'Unknown':
            new = 1.0
        elif check == 'True':
            scoreFile = construct_df(uid)
            df_new = wrapper_scoring_process(scoreFile,dataFolder)
            new = df_new['IWscore']
        col.append(new)
    df_final['specificity_score'] = col
    return df_final



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='breast cancer neoantigen analysis')
    parser.add_argument('--tcga',type=str,default=None,help='task name')
    parser.add_argument('--HLA',type=str,default=None,help='HLA type we want to inspect')
    args=parser.parse_args()

    tcga = args.tcga
    HLA = args.HLA

    global dataFolder
    #tcga = 'TCGA-A1-A0SB'
    #HLA = 'A*11:01 A*25:01 B*35:01 B*44:02 C*04:01 C*05:01'

    if not os.path.exists('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}'.format(tcga)):
        os.makedirs('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}'.format(tcga))

    if not os.path.exists('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}/{0}.ori.txt'.format(tcga)):
        print('loading counts and events')
        counts = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/counts.TCGA-BRCA.txt',sep='\t')
        events = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-filtered-names-75p.txt',sep='\t')
        controls = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/all_matched_control.txt',sep='\t',header=None)
        controls = convert_controls_to_list(controls)   # will be a list
        print('loading only healthy control')
        col_to_use = ['AltAnalyze_ID']
        col_to_use.extend(controls)
        counts_control = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/counts.TCGA-BRCA.txt',sep='\t',usecols=col_to_use)
        #counts_control.to_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/healthy.txt',sep='\t',index=None)
        print('computing mean value in controls')
        mean_df = ave_counts_control(counts_control) 
        #mean_df.to_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/inter.txt',sep='\t',index=None)
        '''
            AltAnalyze_ID   mean
        0
        1   

        '''
        print('obtaining base')
        df_start,part = find_PSI_col(tcga,events)

        print('add counts')
        df_new = add_counts(tcga,df_start,part,counts,mean_df)
        '''
                UID     tcga_tumor    counts     mean_in_control
        0
        1
        '''
        print('filtering')
        df_final = filtering(df_new,tcga)
        df_final.to_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}/{0}.ori.txt'.format(tcga),sep='\t',index=None)


    # let's start to retrieve neoantigens, first, check GTEx as well
    HLA = convert_hla(HLA)  # ['HLA-A11:01','HLA-A25:01','HLA-B35:01','HLA-B44:02','HLA-C04:01','HLA-C05:01']
    outFolder = '/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{}'.format(tcga)
    taskName = tcga
    if not os.path.exists(os.path.join(outFolder,'df_final_{0}.txt'.format(taskName))):
        dataFolder = '/data/salomonis2/LabFiles/Frank-Li/python3/data'
        cutoff_PSI = 0.1
        cutoff_sampleRatio = 0.1
        cutoff_tissueRatio = 0.1
        df = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}/{0}.ori.txt'.format(tcga),sep='\t')
        df_final = check_GTEx_PSI_pbz2(df,cutoff_PSI,cutoff_sampleRatio,cutoff_tissueRatio,taskName)
    else:
        df_final = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}/df_final_{0}.txt'.format(tcga),sep='\t')

    # let's run specificity score
    dataFolder = '/data/salomonis2/LabFiles/Frank-Li/python3/data'
    df_final = run_specificity_score(df_final)


    # then we retrieve the neoantigens
    #HLA = convert_hla(HLA)  # ['HLA-A11:01','HLA-A25:01','HLA-B35:01','HLA-B44:02','HLA-C04:01','HLA-C05:01']
    #df_final = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}/df_final_{0}.txt'.format(tcga),sep='\t')
    #df_final = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/TCGA-A1-A0SB/down.txt',sep='\t')
    df = retrieve(df_final,HLA)
    df = combine_all_mers(df,HLA)
    df = count_all_mers(df)
    df.to_csv('/data/salomonis2/LabFiles/Frank-Li/mhc/TCGA-breast/collated_result/{0}/merged.{0}.txt'.format(tcga),sep='\t',index=None)
