from torch.utils.data import Dataset, DataLoader, random_split,Subset,ConcatDataset
import torch
from torch.nn.utils.rnn import pad_sequence
import torch.nn as nn
import torch.nn.functional as F

import pandas as pd
import numpy as np
from utils import *
import shelve
import itertools




class dilatedCNN(nn.Module):
    def __init__(self,height=21,width=56,channel=1,hidden=128,p=0.2):
        super(dilatedCNN,self).__init__()
        self.height = height
        self.width = width
        self.channel = channel
        self.hidden = hidden
        self.p = p
        
        self.layer1 = nn.Sequential(
            nn.Conv2d(self.channel,16,kernel_size=(21,5),dilation=(1,3)),   
            nn.BatchNorm2d(16),
            nn.ReLU(),
            nn.MaxPool2d((1,3)),
            )
        
        layer1_conv2d_output_H = dilatedCNN.calculate_conv2d_dimension(self.height,21,dilation=1)
        layer1_conv2d_output_W = dilatedCNN.calculate_conv2d_dimension(self.width,5,dilation=3)
        layer1_maxpool2d_output_H = dilatedCNN.calculate_conv2d_dimension(layer1_conv2d_output_H,1,stride=1)
        layer1_maxpool2d_output_W = dilatedCNN.calculate_conv2d_dimension(layer1_conv2d_output_W,3,stride=3)
        
        self.layer2 = nn.Sequential(
            nn.Conv2d(16,32,kernel_size=(1,4),dilation=(1,2)),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d((1,1)),
            )
        
        layer2_conv2d_output_H = dilatedCNN.calculate_conv2d_dimension(layer1_maxpool2d_output_H,1,dilation=1)
        layer2_conv2d_output_W = dilatedCNN.calculate_conv2d_dimension(layer1_maxpool2d_output_W,4,dilation=2)
        layer2_maxpool2d_output_H = dilatedCNN.calculate_conv2d_dimension(layer2_conv2d_output_H,1,stride=1)
        layer2_maxpool2d_output_W = dilatedCNN.calculate_conv2d_dimension(layer2_conv2d_output_W,1,stride=1)
        
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
    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_training_test.txt',sep='\t')
    # hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    # inventory = hla['hla']
    # dic_inventory = dict_inventory(inventory)
    # training_dataset = dataset(ori,hla,dic_inventory)
    # max_length = training_dataset.max_length
    # #training_loader = DataLoader(training_dataset,batch_size=512,shuffle=True,drop_last=False)
    # training_loader = balancedBinaryLoader(training_dataset,batch_size=512)
    # model = dilatedCNN().to(device)
    # optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    # loss_f=nn.CrossEntropyLoss()
    # num_epochs = 50
    # for epoch in range(num_epochs):
    #     loss_list = []
    #     acc_list = []
    #     for i in training_loader:
    #         X = i[0].unsqueeze(1).permute(0,1,3,2).to(device)
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
    #         print(y_pred)
    #         print(predictions)
    #         print(y)
    
    #         num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it
    
    #         num_samples  += predictions.size(0)
    
    #         acc_list.append(float(num_correct)/float(num_samples)*100)
            
    #     loss,acc = sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    

    #     print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(epoch+1,num_epochs,loss,acc))
    #     torch.save(model.state_dict(),'/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/model/dilatedCNN_balance.pth')


    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


    # ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_validation_filter910.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    # hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    # inventory = hla['hla']
    # dic_inventory = dict_inventory(inventory)
    # testing_dataset = dataset(ori,hla,dic_inventory)
    # max_length = 56
    # testing_loader = DataLoader(testing_dataset,batch_size=3279,shuffle=True,drop_last=True)   # 3279 # 652
    # model = dilatedCNN().to(device)
    # model.load_state_dict(torch.load('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/model/dilatedCNN_balance.pth'))

    # model.eval()
    # with torch.no_grad():
    #     for i,(X,y) in enumerate(testing_loader):

    #         x = X.unsqueeze(1).permute(0,1,3,2).to(device)
    #         y = y.to(device)
    #         y_pred = model(x)
    #         print(y_pred)
            
    #         num_correct = 0
    #         num_samples = 0
    #         _,predictions = y_pred.max(1)
    #         print('Predicted:',predictions)
    #         print('True:',y)

    #         num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it

    #         num_samples  += predictions.size(0)

    #         print('Accuracy:',float(num_correct)/float(num_samples)*100)




    #         s = shelve.open('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/testing/testing19')
    #         s['y_pred'] = y_pred.detach().cpu().numpy()
    #         s['predictions'] = predictions.detach().cpu().numpy()
    #         s['y'] = y.detach().cpu().numpy()
    #         s.close()

    # scoring

    ori = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/data/shuffle_validation_filter910.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = dilatedCNN().to(device)
    model.load_state_dict(torch.load('/data/salomonis2/LabFiles/Frank-Li/immunogenecity/transformer/model/dilatedCNN_balance.pth'))
    merList = ['YPALPHDTAI','IPAASQLDL','MPEAMTIVML']
    HLA = ['HLA-B*3501']
    result = construct_df4deeplearningmodel(merList,HLA,model,device)
    print(result)
