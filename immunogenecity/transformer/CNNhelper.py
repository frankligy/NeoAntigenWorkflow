#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 18:40:59 2020

@author: ligk2e
"""

import os
import torch
import torch.nn as nn
from utils import *
import pandas as pd
from torch.utils.data import Dataset, DataLoader, random_split,Subset
import shelve
from skorch import NeuralNet,NeuralNetClassifier
import skorch
import numpy as np
from sklearn.model_selection import GridSearchCV



class CNN(nn.Module):
    def __init__(self,height=14,width=46,channel=25,hidden=128,p=0.4):
        super(CNN,self).__init__()
        self.height = height
        self.width = width
        self.channel = channel
        self.hidden = hidden
        self.p = p
        
        self.layer1 = nn.Sequential(
            nn.Conv2d(self.channel,32,kernel_size=(5,8)),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d((1,3)),
            )
        
        layer1_conv2d_output_H = CNN.calculate_conv2d_dimension(self.height,5)
        layer1_conv2d_output_W = CNN.calculate_conv2d_dimension(self.width,8)
        layer1_maxpool2d_output_H = CNN.calculate_conv2d_dimension(layer1_conv2d_output_H,1,stride=1)
        layer1_maxpool2d_output_W = CNN.calculate_conv2d_dimension(layer1_conv2d_output_W,3,stride=3)
        
        self.layer2 = nn.Sequential(
            nn.Conv2d(32,64,kernel_size=(5,4)),
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d((1,2)),
            )
        
        layer2_conv2d_output_H = CNN.calculate_conv2d_dimension(layer1_maxpool2d_output_H,5)
        layer2_conv2d_output_W = CNN.calculate_conv2d_dimension(layer1_maxpool2d_output_W,4)
        layer2_maxpool2d_output_H = CNN.calculate_conv2d_dimension(layer2_conv2d_output_H,1,stride=1)
        layer2_maxpool2d_output_W = CNN.calculate_conv2d_dimension(layer2_conv2d_output_W,2,stride=2)
        
        self.layer3 = nn.Sequential(
            nn.Dropout(self.p),
            nn.Linear(64 * layer2_maxpool2d_output_H * layer2_maxpool2d_output_W, self.hidden),
            nn.Dropout(self.p),
            nn.ReLU(),
            nn.Linear(self.hidden,2),
            )
        
        self.sigmoid = nn.Sigmoid()
        
    def forward(self,x):  # x [batch,channel,max_H,max_W]
        #print(x.size())
        out = self.layer1(x)
        #print(out.size())
        out = self.layer2(out)
        #print(out.size())
        out = out.reshape(out.shape[0],-1)
        #print(out.size())
        out = self.layer3(out)
        out = self.sigmoid(out)
        return out    # [batch,2]
        
            
        
        
        
    
    @staticmethod
    def calculate_conv2d_dimension(dim_in,kernel,padding=0,dilation=1,stride=1):
        dim_out = (dim_in + 2*padding - dilation*(kernel-1) -1) / stride +1
        return int(dim_out)







if __name__ == '__main__':
    
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


    ori = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/shuffle_training.txt',sep='\t')
    hla = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    
    index = wrapper_preprocess()
    training_dataset = CNN_dataset(ori, hla, dic_inventory, index)
    
    # modelObj,training_dataset,optimizer,criterion,batch_size,num_epochs,outdir
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.CrossEntropyLoss()
    batch_size = 512
    num_epochs = 5
    outdir = '/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/transformer/test.pth'
    pytorch_training(CNN, training_dataset, optimizer, criterion, batch_size, num_epochs, outdir)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    