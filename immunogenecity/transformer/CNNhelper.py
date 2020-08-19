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
    def __init__(self,height=10,width=46,channel=25,hidden=128,p=0.4):
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
            nn.Conv2d(32,32,kernel_size=(2,4)),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d((1,2)),
            )
        
        layer2_conv2d_output_H = CNN.calculate_conv2d_dimension(layer1_maxpool2d_output_H,2)
        layer2_conv2d_output_W = CNN.calculate_conv2d_dimension(layer1_maxpool2d_output_W,4)
        layer2_maxpool2d_output_H = CNN.calculate_conv2d_dimension(layer2_conv2d_output_H,1,stride=1)
        layer2_maxpool2d_output_W = CNN.calculate_conv2d_dimension(layer2_conv2d_output_W,2,stride=2)
        
        self.layer3 = nn.Sequential(
            nn.Dropout(self.p),
            nn.Linear(32 * layer2_maxpool2d_output_H * layer2_maxpool2d_output_W,self.hidden),
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


    ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_training_test.txt',sep='\t')
    hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    
    index = wrapper_preprocess()
    training_dataset = CNN_dataset(ori, hla, dic_inventory, index)
    print(training_dataset[1][0].size())
    
    model = CNN().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.CrossEntropyLoss()


    training_loader = DataLoader(training_dataset,batch_size=512,shuffle=True,drop_last=False)

    num_epochs = 50
    for epoch in range(num_epochs):
        loss_list = []
        acc_list = []
        for i in training_loader:
            X = i[0].to(device)
            y = i[1].to(device)
            optimizer.zero_grad()
            
            y_pred = model(X)
            print(y_pred)
            loss = criterion(y_pred,y)
            loss.backward()
            optimizer.step()
            loss_list.append(loss.item())
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)
            print(predictions)
            print(y)
    
            num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it
    
            num_samples  += predictions.size(0)
    
            acc_list.append(float(num_correct)/float(num_samples)*100)
            
        loss,acc = sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    

        print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(epoch+1,num_epochs,loss,acc))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    