import os
import torch
import torch.nn as nn
from utils import *
import pandas as pd
from torch.utils.data import Dataset, DataLoader, random_split,Subset
import numpy as np




class FC(nn.Module):
    def __init__(self):
        super(FC,self).__init__()
        self.layer1 = nn.Linear(11500,1024)
        self.layer2 = nn.Linear(1024,128)
        self.layer3 = nn.Linear(128,2)
        self.sigmoid = nn.Sigmoid()
        self.dropout = nn.Dropout(0.3)
        self.relu = nn.ReLU()

    def forward(self,x):   # x [batch,11500]
        out = self.relu(self.dropout(self.layer1(x)))
        out = self.relu(self.dropout(self.layer2(out)))
        out = self.relu(self.dropout(self.layer3(out)))
        out = self.sigmoid(out)
        return out    # [batch,2]

if __name__ == '__main__':
    
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


    ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_training_test.txt',sep='\t')
    hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    
    index = wrapper_preprocess()
    training_dataset = FC_dataset(ori, hla, dic_inventory, index)
    
    model = FC().to(device)
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