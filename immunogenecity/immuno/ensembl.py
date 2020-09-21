import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from aaindex_encoding_ResLikeCNN import model_aaindex,construct_aaindex,pull_label_aaindex,pull_hla_aaindex,pull_peptide_aaindex
from seperateCNN import *
from ResLikeCNN import *
from utils import *
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score


if __name__ == '__main__':
    # aaindex model result
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',
                           sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct_aaindex(ori_test, hla, dic_inventory,after_pca)
    input1_test = pull_peptide_aaindex(dataset_test)
    input2_test = pull_hla_aaindex(dataset_test)
    label_test = pull_label_aaindex(dataset_test)
    aaindex_model = model_aaindex()
    aaindex_model.load_weights('aaindex12_encoding_ReslikeCNN')
    result_a = aaindex_model.predict(x=[input1_test, input2_test])
    hard_a = [1 if i > 0.5 else 0 for i in result_a]

    # ResLIkecnn RESULT
    ResLikeCNN = model()
    ResLikeCNN.load_weights('ResLikeCNN_epoch8_sigmoid/')
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    #seperateCNNmodel.evaluate(x=[input1_test,input2_test],y=label_test,batch_size=512)
    result = ResLikeCNN.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result]


    # compose them together
    df = pd.DataFrame({'label':label_test[:,0],'aaindex':result_a[:,0],'reslike':result[:,0]})

    # get gan
    import pickle
    with open('/Users/ligk2e/Desktop/tmp/GAN_result.p','rb') as f1:
        gan = pickle.load(f1)
    df['gan'] = gan


    combine = np.mean(np.array([df['aaindex'].values,df['neg'].values]),axis=0)
    df['combined'] = combine
    #df.to_csv('ensemble_result.txt',sep='\t',index=None)
    ensemble = [1 if i > 0.5 else 0 for i in combine]
    confusion_matrix(label_test,ensemble)
    draw_ROC(label_test,combine)
    draw_PR(label_test,combine)





