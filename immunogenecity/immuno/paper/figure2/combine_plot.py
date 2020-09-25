import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import pandas as pd
import numpy as np
from seperateCNN import *
from utils import *
import shelve
from ResLikeCNN import *
from aaindex_encoding_ResLikeCNN import model_aaindex,construct_aaindex,pull_hla_aaindex,pull_label_aaindex,pull_peptide_aaindex

def draw_combined_ROC(arrayP,arrayT):
    fig = plt.figure()
    colormap = ['red','green','blue','cyan','magenta','orange']
    legend = ['iedb','logistic regression','deephlapan','PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)']
    # colormap = ['cyan','magenta','orange']
    # legend = ['PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)']
    for i in range(len(arrayP)):
        fpr,tpr,_ = roc_curve(arrayT[i],arrayP[i],pos_label=1)
        area = auc(fpr, tpr)
        lw = 2
        plt.plot(fpr, tpr, color=colormap[i],
                 lw=lw, label='{0} (area = {1:0.2f})'.format(legend[i],area))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.show()
    plt.savefig('paper/figure2/testing_ROC.pdf')

def draw_combined_PR(arrayP,arrayT):
    fig = plt.figure()
    colormap = ['red', 'green', 'blue', 'cyan', 'magenta', 'orange']
    legend = ['iedb', 'logistic regression', 'deephlapan', 'PIER(Blosum Encoding)', 'PIER(AAindex Encoding)',
              'PIER(Ensemble)']
    # colormap = ['cyan','magenta','orange']
    # legend = ['PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)']
    for i in range(len(arrayP)):
        precision,recall,_ = precision_recall_curve(arrayT[i],arrayP[i],pos_label=1)
        area = auc(recall, precision)
        lw = 2
        plt.plot(recall, precision, color=colormap[i],
                 lw=lw, label='{0} (area = {1:0.2f})'.format(legend[i],area))
    plt.plot([0, 1], [0.12, 0.12], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc="best")
    plt.show()
    plt.savefig('paper/figure2/testing_pr.pdf')

if __name__ == '__main__':
    # load testing dataset
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)




    # result from deepimmuno
    seperateCNNmodel = seperateCNN()
    seperateCNNmodel.load_weights('secondFilter32_epoch150/')
    result = seperateCNNmodel.predict(x=[input1_test,input2_test])

    # result from baseline logistic regression
    s = shelve.open('logistic')
    y_pred = s['y_pred']
    Y_test = s['Y_test']
    s.close()
    y_pred = y_pred[:,1]
    Y_test = Y_test.astype(np.int32)

    # result from deephlapan
    df1 = pd.read_csv('deephlapan/ineo_testing_new.txt', sep='\t')
    df2 = pd.read_csv('deephlapan/ineo_testing_new_final_predicted_result.csv')
    y_deephlapan = df1['immunogenecity'].values
    y_pred_deephlapan = df2['immunogenic score'].values

    # result from IEDB
    iedb = pd.read_csv('IEDB.csv')
    y_iedb = iedb['label'].tolist()
    y_pred_iedb = iedb['score'].tolist()

    # result from ResLikeCNN
    ResLikeCNN = model()
    ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')
    res = ResLikeCNN.predict([input1_test,input2_test])

    # result from aaindex12_ResLikeCNN
    ResLikeCNN_index = model_aaindex()
    ResLikeCNN_index.load_weights('aaindex12_encoding_ReslikeCNN_reproduce/')
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',
                           sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct_aaindex(ori_test, hla, dic_inventory,after_pca)
    input1_test_a = pull_peptide_aaindex(dataset_test)
    input2_test_a = pull_hla_aaindex(dataset_test)
    label_test_a = pull_label_aaindex(dataset_test)
    res_aaindex = ResLikeCNN_index.predict(x=[input1_test_a, input2_test_a])

    # do combined ensemble, first run enselbl.py
    #df = pd.read_csv('ensemble_result.txt',sep='\t')
    result_ensemble = np.mean(np.concatenate([res,res_aaindex],axis=1),axis=1)

    # # result from transfer_learning, well, if you count this best-performed one
    # ResLikeCNN = model()
    # ResLikeCNN.load_weights('transfer_learning_best_performance_achieved')
    # tran = ResLikeCNN.predict([input1_test,input2_test])

    # let's draw the figure
    arrayP = [y_pred_iedb,y_pred,y_pred_deephlapan,res,res_aaindex,result_ensemble]
    arrayT = [y_iedb,Y_test,y_deephlapan,label_test,label_test,label_test]
    draw_combined_ROC(arrayP,arrayT)
    draw_combined_PR(arrayP,arrayT)



