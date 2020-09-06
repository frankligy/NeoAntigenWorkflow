import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
from utils import *
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import matplotlib.pyplot as plt
import numpy as np

from seperateCNN import *
from ResLikeCNN import *




if __name__ == '__main__':
    # load training binding data
    ori = pd.read_csv('data/transfer_training.txt',sep='\t')
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    dataset = construct(ori, hla, dic_inventory)   # [ (10,21,1),(46,21,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    label = pull_label(dataset)

    # load validation binding data
    ori_val = pd.read_csv('data/transfer_validation.txt',sep='\t')
    dataset_val = construct(ori_val, hla, dic_inventory)
    input1_val = pull_peptide(dataset_val)
    input2_val = pull_hla(dataset_val)
    label_val = pull_label(dataset_val)

    # do training
    ResLikeCNN = model()
    ResLikeCNN.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )   # if you want to do model.summary, just run a toy case to let the model pick up the input dimension

    history = ResLikeCNN.fit(
        x=[input1,input2],   # feed a list into
        y=label,
        validation_data = ([input1_val,input2_val],label_val),
        batch_size=512,
        epochs=11,
        #class_weight = {0:0.2,1:0.8}   # I have 20% positive and 80% negative in my training data
    )
    draw_history(history)

    ResLikeCNN.save_weights('trained_on_binding')


    # let's test it
    ori_test = pd.read_csv('data/transfer_testing.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)

    result = ResLikeCNN.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result]
    confusion_matrix(label_test,hard)
    f1_score(label_test,hard)
    accuracy_score(label_test,hard)
    draw_ROC(label_test,result)
    draw_PR(label_test,result)
    '''
    The performance is pretty bad, not sure whether due to testing-dataset or anything else.
    '''

    # let's do transfer learning and fine-tuning
    print(ResLikeCNN.summary())
    '''
    Model: "model_2"
    _________________________________________________________________
    Layer (type)                 Output Shape              Param #   
    =================================================================
    cnn_peptide_2 (CNN_peptide)  multiple                  46496     
    _________________________________________________________________
    cnn_mhc_2 (CNN_MHC)          multiple                  116704    
    _________________________________________________________________
    flatten_2 (Flatten)          multiple                  0         
    _________________________________________________________________
    dense_4 (Dense)              multiple                  32896     
    _________________________________________________________________
    dense_5 (Dense)              multiple                  129       
    =================================================================
    Total params: 196,225
    Trainable params: 194,177
    Non-trainable params: 2,048
    _________________________________________________________________  
    '''
    ## first experiment, don't freeze any layer, just use binding data as a start point
    ## second experiemnt, freeze all conv layer, fine-tune dense layers
    ResLikeCNN = model()

    #ResLikeCNN.load_weights('trained_on_binding')
    ResLikeCNN.load_weights('transfer_learning_best_performance_achieved')
    ResLikeCNN.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )  # if you want to do model.summary, just run a toy case to let the model pick up the input dimension
    ResLikeCNN.layers[0].trainable = False
    ResLikeCNN.layers[1].trainable = False
    ResLikeCNN.layers[2].trainable = False
    ResLikeCNN.layers[3].trainable = False

    ori_retrain = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_retrain = construct(ori_retrain, hla, dic_inventory)
    input1_retrain = pull_peptide(dataset_retrain)
    input2_retrain = pull_hla(dataset_retrain)
    label_retrain = pull_label(dataset_retrain)

    ori_retrain_reval = pd.read_csv('data/shuffle_validation_filter910.txt',sep='\t')
    dataset_retrain_reval = construct(ori_retrain_reval, hla, dic_inventory)
    input1_retrain_reval = pull_peptide(dataset_retrain_reval)
    input2_retrain_reval = pull_hla(dataset_retrain_reval)
    label_retrain_reval = pull_label(dataset_retrain_reval)


    history_retrain = ResLikeCNN.fit(
        x=[input1_retrain,input2_retrain],   # feed a list into
        y=label_retrain,
        validation_data = ([input1_retrain_reval,input2_retrain_reval],label_retrain_reval),
        batch_size=512,
        epochs=10,
        class_weight = {0:0.2,1:0.8}   # I have 20% positive and 80% negative in my training data
    )

    # let's test it
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    #seperateCNNmodel.evaluate(x=[input1_test,input2_test],y=label_test,batch_size=512)
    result = ResLikeCNN.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result]
    confusion_matrix(label_test,hard)
    f1_score(label_test,hard)
    accuracy_score(label_test,hard)
    draw_ROC(label_test,result)
    draw_PR(label_test,result)





