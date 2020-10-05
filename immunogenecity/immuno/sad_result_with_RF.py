import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import pandas as pd
import numpy as np
import shelve


def draw_combined_ROC(arrayP,arrayT):
    fig = plt.figure()
    colormap = ['green','black','cyan','magenta','orange','red']
    legend = ['random_forest(blosum)','random_forest(aaindex)',
              'PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)','all']
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
    plt.savefig('paper/supplementary figure2/bench_ROC.pdf')

def draw_combined_PR(arrayP,arrayT):
    fig = plt.figure()
    # colormap = ['red','green','blue','black','cyan','magenta','orange','yellow']
    # legend = ['ANN','random_forest(max=2)','random_forest(max=10)','random_forest(no depth restriction)',
    #           'PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)','all']
    colormap = ['green','black','cyan','magenta','orange','red']
    legend = ['random_forest(blosum)','random_forest(aaindex)',
              'PIER(Blosum Encoding)','PIER(AAindex Encoding)','PIER(Ensemble)','all']
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
    plt.savefig('paper/supplementary figure2/bench_pr.pdf')


def draw_ROC(y_true,y_pred):

    fpr,tpr,_ = roc_curve(y_true,y_pred,pos_label=1)
    area_mine = auc(fpr,tpr)
    fig = plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
            lw=lw, label='ROC curve (area = %0.2f)' % area_mine)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

def draw_PR(y_true,y_pred):
    precision,recall,cutoff = precision_recall_curve(y_true,y_pred,pos_label=1)
    area_PR = auc(recall,precision)
    baseline = np.sum(np.array(y_true) == 1) / len(y_true)

    plt.figure()
    lw = 2
    plt.plot(recall,precision, color='darkorange',
            lw=lw, label='PR curve (area = %0.2f)' % area_PR)
    plt.plot([0, 1], [baseline, baseline], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve example')
    plt.legend(loc="lower right")
    plt.show()
    return cutoff

if __name__ == '__main__':
    rf_df = pd.read_csv('rf_no_depth.txt',sep='\t')
    pier_df = pd.read_csv('ensemble_result.txt',sep='\t')

    y_label = rf_df['label']
    y_pred_rf_blosum = rf_df['rf_blosum']
    y_pred_rf_aaindex = rf_df['rf_aaindex']
    y_pred_rf_max2 = rf_df['aaindex_max2']
    y_pred_rf_max10 = rf_df['aaindex_max10']
    y_pred_rf_max50 = rf_df['aaindex_max50']

    y_pred_pier_blosum = pier_df['blosum']
    y_pred_pier_aaindex = pier_df['aaindex']
    y_pred_pier_ensemble = pier_df['ensemble']

    # compare
    df_tmp = pd.concat([y_label,y_pred_rf_aaindex,y_pred_rf_blosum,y_pred_pier_blosum,y_pred_pier_aaindex],axis=1)
    df_tmp.to_csv('/Users/ligk2e/Desktop/tmp/rf_pier.txt',sep='\t',index=None)

    # let's simply ensemble them together
    mat = np.empty([len(y_label),4])
    mat[:,0] = y_pred_rf_blosum
    mat[:,1] = y_pred_rf_aaindex
    mat[:,2] = y_pred_pier_blosum
    mat[:,3] = y_pred_pier_aaindex
    result = np.mean(mat,axis=1)
    draw_ROC(y_label,result)
    draw_PR(y_label,result)

    s = shelve.open('ANN')
    y_pred = s['y_pred']
    Y_test = s['Y_test']
    s.close()
    y_pred_ANN = y_pred[:,1]
    Y_test_ANN = Y_test.astype(np.int32)

    # let's draw the figure
    arrayP = [y_pred_rf_blosum,y_pred_rf_aaindex,y_pred_pier_blosum,y_pred_pier_aaindex,y_pred_pier_ensemble,result]
    arrayT = [y_label,y_label,y_label,y_label,y_label,y_label]
    draw_combined_ROC(arrayP,arrayT)
    draw_combined_PR(arrayP,arrayT)
