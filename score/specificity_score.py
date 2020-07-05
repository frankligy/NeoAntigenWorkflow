import os
import pandas as pd
import numpy as np
import bz2
import _pickle as cpickle
import statistics
from sklearn import preprocessing
import pickle







# PSI portion, average PSI, median PSI, percentage above user-defined cutoff(0.1)
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

def training_process(training,cutoff_PSI,cutoff_readcounts):   
    # get value matrix
    ave_PSI,median_PSI,percentage_PSI = PSI(training,cutoff_PSI)
    ave_counts,median_counts,percentage_counts = ReadCounts(training,cutoff_readcounts)

    mat_ori = np.column_stack((ave_PSI,median_PSI,percentage_PSI,ave_counts,median_counts,percentage_counts))  # original giant matrix
    max_mat_ori = np.amax(mat_ori,axis=0)   # store max value for each metric
    min_mat_ori = np.amin(mat_ori,axis=0)   # store min value for each metric
    print('shape of original matrix is:',mat_ori.shape)
    print('original matrix is:\n',mat_ori)
    with open('mat_ori.p','wb') as f1:
        pickle.dump(mat_ori,f1)

    ## L2 normalize ave_counts and median_counts
    old = np.column_stack((ave_counts,median_counts))  # (n*2) array
    new = preprocessing.normalize(old,norm='l2',axis=0)   # column-wise,axis=0
    ave_counts_new = new[:,0]
    median_counts_new = new[:,1]

    mat_rescale_before = np.column_stack((ave_PSI,median_PSI,percentage_PSI,ave_counts_new,median_counts_new,percentage_counts))   # giant matrix storing value for all 6 metrics
    mat_rescale = preprocessing.scale(mat_rescale_before,axis=0)  # to mean=0,var = 1
    max_mat_rescale = np.amax(mat_rescale,axis=0)   # store max value for each metric
    min_mat_rescale = np.amin(mat_rescale,axis=0)   # store min value for each metric
    print('shape of rescaled matrix is:',mat_rescale.shape)
    print('rescaled matrix is:\n',mat_rescale)
    with open('mat_rescale.p','wb') as f2:
        pickle.dump(mat_rescale,f2)

    # get covariance matrix
    covariance_matrix = np.cov(np.transpose(mat_rescale))   # np.cov recognize random variable by row, so need to transpose the matrix
    print('shape of covariance matrix is:',covariance_matrix.shape)
    print('covariance_matrix is:\n', covariance_matrix)
    with open('mat_cov.p','wb') as f3:
        pickle.dump(covariance_matrix,f3)

    # get leading eigenvector (largest eigenvalue)
    eigen = np.linalg.eig(covariance_matrix)
    leading_eigenvector = eigen[1][:,np.argmax(eigen[0])]   # index of largest eigenvalue corresponds to the column of second array
    print('shape of leading_eigenvector:',leading_eigenvector.size)
    # order: ave_PSI,median_PSI,percentage_PSI,ave_counts,median_counts,percentage_counts
    return max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector



def scoring_process(df_evaluation,cutoff_PSI,cutoff_readcounts,max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector):
    ave_PSI,median_PSI,percentage_PSI = PSI(df_evaluation,cutoff_PSI)
    ave_counts,median_counts,percentage_counts = ReadCounts(df_evaluation,cutoff_readcounts)
    
    mat_eval = np.column_stack((ave_PSI,median_PSI,percentage_PSI,ave_counts,median_counts,percentage_counts))
    print('shape of mat_eval:',mat_eval.shape)
    print(mat_eval)
    mat_eval_new = np.zeros(mat_eval.shape)
    print('shape of mat_eval_new:',mat_eval_new.shape)
    print(mat_eval_new)
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



if __name__ == '__main__':
    # load training set 
    dataFolder = '/data/salomonis2/LabFiles/Frank-Li/python3/data'
    training = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/specificity-score/df_train_breast.txt',sep='\t')
    # max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector = training_process(training,0.1,3)
    # print(max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector)

    max_mat_ori = [9.98564855e-01,1.00000000e+00,1.00000000e+00,1.24263241e+04,1.02370000e+04,1.00000000e+00]
    min_mat_ori = [0.00092801,0.0,0.00326531,0.01387755,0.0,0.0]
    max_mat_rescale = [1.6203541,1.246249,0.98267483,72.70393268,80.23512846,0.77875848]
    min_mat_rescale = [-1.69933277,-1.43360037,-2.65506607,-0.30237015,-0.30285234,-2.89327914]
    leading_eigenvector = [0.48264742,0.47174347,0.47383551,0.21692184,0.22945297,0.46934607]
    df_new = scoring_process(training,0.1,3,max_mat_ori,min_mat_ori,max_mat_rescale,min_mat_rescale,leading_eigenvector)
    df_new.to_csv('/data/salomonis2/LabFiles/Frank-Li/specificity-score/df_new_breast.txt',sep='\t',index=None)
    

