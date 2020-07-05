#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 10:08:07 2020

@author: ligk2e
"""

import torch
import pickle
import numpy as np
from sklearn.model_selection import train_test_split
import sklearn

N, D_in, H, D_out = 20,320,30,1

with open('/Users/ligk2e/Desktop/immunogenecity/X.p','rb') as f:
    X = pickle.load(f)
    
with open('/Users/ligk2e/Desktop/immunogenecity/Y.p','rb') as f:
    Y = pickle.load(f)
    
X_train,X_test,Y_train,Y_test = train_test_split(X, Y,test_size = 0.2,random_state = 5)


X_feed = np.array_split(X_train,9)
Y_feed = np.array_split(Y_train,9)

X_feed, Y_feed= np.array(X_feed),np.array(Y_feed)


class TwoLayerNet(torch.nn.Module):
    def __init__(self,D_in,H,D_out):
        super(TwoLayerNet,self).__init__()
        self.linear1 = torch.nn.Linear(D_in,H)
        self.linear2 = torch.nn.Linear(H,D_out)
        
    def forward(self,x):   # accept a Tensor of input data and return a Tensor of output data
        h_relu = self.linear1(x).clamp(min=0)    # clamp(min=0) is exactly ReLu
        y_pred = self.linear2(h_relu)
        sigmoid = torch.nn.Sigmoid()
        y_pred_sigmoid = sigmoid(y_pred)
        return y_pred_sigmoid
    
model = TwoLayerNet(D_in, H, D_out)    

# model = torch.nn.Sequential(
#     torch.nn.Linear(D_in,H),
#     torch.nn.ReLU(),
#     torch.nn.Linear(H,D_out),
#     )


loss_fn = torch.nn.MSELoss(reduction='sum')
learning_rate = 1e-4  
optimizer = torch.optim.Adam(model.parameters(),lr=learning_rate)                  

for t in range(50):    # epoch
    for i in range(9):   # batch size
        x = torch.from_numpy(X_feed[i]).float()
        y = torch.from_numpy(Y_feed[i]).float().view(-1,1)
        
        y_pred = model(x)
        loss = loss_fn(y_pred,y)
        if t% 100 == 99: print(t,loss.item())    # every 100 epoches, print the loss
        optimizer.zero_grad()
        loss.backward()   # calculate gradients
        optimizer.step()


        
pred = model(torch.from_numpy(X_test).float()).detach().numpy().reshape([42])

sklearn.metrics.roc_auc_score(Y_test,pred)

precision,recall,threshoulds = sklearn.metrics.precision_recall_curve(Y_test,pred)

import matplotlib.pyplot as plt
plt.plot(recall,precision)

np.trapz(precision,dx=0.05)

















