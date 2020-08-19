#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 21:17:34 2020

@author: ligk2e
"""
import os
import torch
import torch.nn as nn
from utils import dict_inventory,rescue_unknown_hla,dataset,balancedBinaryLoader
import pandas as pd
from torch.utils.data import Dataset, DataLoader, random_split,Subset
import shelve
from skorch import NeuralNet,NeuralNetClassifier
import skorch
import numpy as np
from sklearn.model_selection import GridSearchCV



class SelfAttention(nn.Module):
    def __init__(self,embed_size,heads):
        super(SelfAttention,self).__init__()
        self.embed_size = embed_size
        self.heads = heads
        self.head_dim = embed_size // heads
        
        assert(self.head_dim * self.heads == self.embed_size), "Emdedding size is not divisible by heads"
        
        self.values = nn.Linear(self.head_dim,self.head_dim,bias=False)
        self.keys = nn.Linear(self.head_dim,self.head_dim,bias=False)
        self.queries = nn.Linear(self.head_dim,self.head_dim,bias=False)
        self.fc_out = nn.Linear(heads * self.head_dim,embed_size)
        
    def forward(self,values,keys,queries,mask):
        # values, keys, queries are all [batch,seq_len,embedding_size]
        batch = values.shape[0]
        value_len =values.shape[1]
        key_len = keys.shape[1]
        query_len = queries.shape[1]
        
        # reshape them to sub-chunk, [batch,seq_len,heads,head_dim]
        values = values.reshape(batch,value_len,self.heads,self.head_dim)
        keys = keys.reshape(batch,key_len,self.heads,self.head_dim)
        queries = queries.reshape(batch,query_len,self.heads,self.head_dim)
       
        # these three linear layer doesn't change dimension, [batch,seq_len,heads,head_dim]
        values = self.values(values)
        keys = self.keys(keys)
        queries = self.queries(queries)
        
        # energy [batch, heads, seq_len,seq_len], energy is the product of query * keys.T

        energy = torch.einsum('bqhd,bkhd -> bhqk', [queries,keys])  # q means query_len, k means key_len
        
        if mask is not None:
            energy = energy.masked_fill(mask == 0, float("-1e20"))
        
        # attention [batch,heads,seq_len,seq_len]
        attention = torch.softmax(energy / (self.embed_size ** (1/2)), dim=3)
        
        # out is the product of attention * v, [batch,seq_len,heads,head_dim], then flatten [batch,seq_len,embed_size]
        out = torch.einsum('bhqk,bqhd -> bqhd',[attention,values])  # q is value_len, also query_len since they are the same
        out = out.reshape(batch,value_len,self.heads*self.head_dim)
        # this linear doesn't change dimemsion,[batch,seq_len,embed_size]
        out = self.fc_out(out)
        return out
    
    
class TransformerBlock(nn.Module):
    # this block doesn't change dimension
    def __init__(self,embed_size,heads,dropout,forward_expansion):
        super(TransformerBlock,self).__init__()
        self.attention = SelfAttention(embed_size,heads)
        self.norm1 = nn.LayerNorm(embed_size)
        self.norm2 = nn.LayerNorm(embed_size)
        
        self.feed_forward = nn.Sequential(
            nn.Linear(embed_size,forward_expansion*embed_size),
            nn.ReLU(),
            nn.Linear(forward_expansion*embed_size,embed_size),
            )
        
        self.dropout = nn.Dropout(dropout)
        
    
    def forward(self,value,key,query,mask):
        attention = self.attention(value,key,query,mask)
        
        # add skip connection
        x = self.dropout(self.norm1(attention+query))
        forward = self.feed_forward(x)
        # add skip conenction
        out= self.dropout(self.norm2(forward+x))
        
        return out

class Encoder(nn.Module):
    def __init__(self,seq_len,embed_size,num_layers,heads,forward_expansion,dropout,max_length):
        super(Encoder,self).__init__()
        self.embed_size = embed_size
        self.position_embedding = nn.Embedding(max_length,embed_size)
        self.layers = nn.ModuleList(
            [ TransformerBlock(embed_size, heads, dropout, forward_expansion) for _ in range(num_layers)]
            )

        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(embed_size*seq_len,2)

        self.sigmoid = nn.Sigmoid()
        
    def forward(self,x):          # x [batch,seq_len,embed_size]
 
        #print(x.size())
        mask = None
        batch,seq_len = x.shape[0],x.shape[1]
        positions = torch.arange(0,seq_len).expand(batch,seq_len).to(device)  # [batch,seq_len], entry is the indices
        out = self.dropout(x + self.position_embedding(positions))
        #print(out.size())
        for layer in self.layers:
            out = layer(out,out,out,mask)
        
        #print(out.size())
        out = out.reshape(out.shape[0],-1)
        #print(out.size())
        out = self.fc(out)
        out = self.sigmoid(out)
        
        return out   # [batch,2]
            
        

                



if __name__ == '__main__':

    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_training.txt',sep='\t')
    # hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable.txt',sep='\t',header=None,names=['hla','paratope'])
    # inventory = hla['hla']
    # dic_inventory = dict_inventory(inventory)
    # training_dataset = dataset(ori,hla,dic_inventory)
    # max_length = training_dataset.max_length
    
    
    # net = NeuralNetClassifier(
    #     # dataset=uninitialized torch.Dataset
    #     # iterator_train=uninitialized torch.DataLoader
    #     # iterator_valid=uninitialized torch.DataLoader
    #     module=Encoder,
    #     module__seq_len=max_length,
    #     module__embed_size=20,
    #     module__num_layers=6,
    #     module__heads=1,
    #     module__forward_expansion=4,
    #     module__dropout=0.1,
    #     module__max_length=max_length,
    #     iterator_train__batch_size=512,
    #     iterator_train__shuffle=True,
    #     train_split=skorch.dataset.CVSplit(10),
    #     device = device,
    #     optimizer=torch.optim.Adam,
    #     optimizer__lr=0.0001,
    #     criterion=nn.CrossEntropyLoss,
    #     max_epochs=20,

    #     )    
    
    

    
    # X = torch.stack([item[0] for item in training_dataset],dim=0)  # x [N,seq_len,20]
        
    # y = torch.stack([item[1] for item in training_dataset],dim=0)   # y: [N], 1d
    # # net.fit(X,y)
    # # predictions = net.predict(X)
    
    # # GridSearchCV
    # params = {
    #     'module__forward_expansion':[4,8,12],
    #     'module__num_layers':[6,10,14],
    #     'module_dropout': [0,0.1],
    #     'iterator_train__batch_size':[32,128,512,1024],
    #     'optimizer__lr':[0.01,0.001,0.0001],
    #     }
    # gs = GridSearchCV(net, params, refit=False, scoring='accuracy')
    # gs.fit(X,y)
    # print(gs.best_score_, gs.best_params_)








    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


    ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_validation_filter910.txt',sep='\t')
    hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    testing_dataset = dataset(ori,hla,dic_inventory)
    max_length = 56
    testing_loader = DataLoader(testing_dataset,batch_size=3279,shuffle=True,drop_last=True)
    model = Encoder(max_length,21,12,1,10,0,max_length).to(device)
    model.load_state_dict(torch.load('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/model/restrict10.pth'))

    model.eval()
    with torch.no_grad():
        for i,(X,y) in enumerate(testing_loader):

            x = X.to(device)
            y = y.to(device)
            y_pred = model(x)
            print(y_pred)
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)
            print('Predicted:',predictions)
            print('True:',y)

            num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it

            num_samples  += predictions.size(0)

            print('Accuracy:',float(num_correct)/float(num_samples)*100)

            # y_pred_post = y_pred.detach().cpu().numpy()
            # y_post = y.detach().cpu().numpy()


            s = shelve.open('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/testing/testing14')
            s['y_pred'] = y_pred.detach().cpu().numpy()
            s['predictions'] = predictions.detach().cpu().numpy()
            s['y'] = y.detach().cpu().numpy()
            s.close()

    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_training_test.txt',sep='\t')
    # hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    # inventory = hla['hla']
    # dic_inventory = dict_inventory(inventory)
    
    # training_dataset = dataset(ori,hla,dic_inventory)
    # max_length = training_dataset.max_length
    
    # training_loader = DataLoader(training_dataset,batch_size=512,shuffle=True,drop_last=False)
    # #training_loader = balancedBinaryLoader(training_dataset,batch_size=512)
    
    # model = Encoder(max_length,21,12,1,10,0,max_length).to(device)
    # optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    # loss_f=nn.CrossEntropyLoss()
    
    
    # num_epochs = 150
    # for epoch in range(num_epochs):
    #     loss_list = []
    #     acc_list = []
    #     for i in training_loader:
    #         X = i[0].to(device)
    #         y = i[1].to(device)
    #         optimizer.zero_grad()
            
    #         y_pred = model(X)
    #         loss = loss_f(y_pred,y)
    #         loss.backward()
    #         optimizer.step()
    #         loss_list.append(loss.item())
            
    #         num_correct = 0
    #         num_samples = 0
    #         _,predictions = y_pred.max(1)
    
    #         num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it
    
    #         num_samples  += predictions.size(0)
    
    #         acc_list.append(float(num_correct)/float(num_samples)*100)
            
    #     loss,acc = sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    

    #     print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(epoch+1,num_epochs,loss,acc))

    # torch.save(model.state_dict(),'/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/model/restrict10.pth')
    




























    
        
        
        
        