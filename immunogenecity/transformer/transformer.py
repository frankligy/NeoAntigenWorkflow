#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 21:17:34 2020

@author: ligk2e
"""
import os
import torch
import torch.nn as nn
from utils import dict_inventory,rescue_unknown_hla,dataset
import pandas as pd
from torch.utils.data import Dataset, DataLoader, random_split,Subset


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
        self.fc_final = nn.Linear(embed_size*seq_len,2)
        
    def forward(self,x,mask):
        # x [batch,seq_len,embed_size]
        batch,seq_len = x.shape[0],x.shape[1]
        positions = torch.arange(0,seq_len).expand(batch,seq_len).to(device)  # [batch,seq_len], entry is the indices
        out = self.dropout(x + self.position_embedding(positions))
        for layer in self.layers:
            out = layer(out,out,out,mask)
        
        out = out.reshape(out.shape[0],-1)
        out = self.fc_final(out)
        
        return out   # [batch,2]
            
        
        


        
# model = SelfAttention(21,1) 
# model_tran = TransformerBlock(21, 1, 0.1, 4)
# y = model(a,a,a,None)   
# y1 = model_tran(y,y,y,None)

# model_all = Encoder(79,21,6,1,4,0.1,79).to(device)
# a = torch.rand([64,79,21]).to(device)
# y2 = model_all(a,None)

if __name__ == '__main__':
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    ori = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data_new.txt',sep='\t')
    hla = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla2paratopeTable.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    
    training_dataset = dataset(ori,hla,dic_inventory)
    max_length = training_dataset.max_length
    
    training_loader = DataLoader(training_dataset,batch_size=512,shuffle=True,drop_last=True)
    
    model = Encoder(max_length,20,6,1,4,0.1,max_length).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    scheduler = scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, factor=0.1, patience=2, verbose=True)
    # if it observe a non-decreasing loss, give you another (patience-1) more chances, if still not decrease, will reduce 
    # learning rate to factor*learning rate
    loss_f=nn.CrossEntropyLoss()
    
    
    num_epochs = 3
    for epoch in range(num_epochs):
        loss_list = []
        acc_list = []
        for _,(X,y) in enumerate(training_loader):
    
    
            X = X.to(device)   # [batch,seq_len,20]
            y = y.to(device)   # [batch]
            #print(X.size(),y,size())
            optimizer.zero_grad()
            
            y_pred = model(X,None)
            loss = loss_f(y_pred,y)
            loss.backward()
            optimizer.step()
            loss_list.append(loss.item())
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)
    
            num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it
    
            num_samples  += predictions.size(0)
    
            acc_list.append(float(num_correct)/float(num_samples)*100)
            
        loss,acc = sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    
        scheduler.step(loss)
        print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(epoch+1,num_epochs,loss,acc))
    




























    
        
        
        
        