#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 18:10:01 2020

@author: ligk2e
"""

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
barwidth = 0.9

r = np.arange(12)
d = [48,268,83,44,90,64,31,46,41,36,146,185]

label = ['OvCa48','OvCa53','OvCa58','OvCa64','OvCa65','OvCa70','OvCa80','OvCa84','OvCa99',
         'OvCa104','OvCa105','OvCa109']

plt.bar(r,d,width=barwidth)
for i in range(12):
    plt.text(x=i-0.3,y=d[i]+5,s=d[i],fontsize=8)
plt.xticks(r,label,rotation=60)
plt.title('Number of 12 Ovarian Cancer patients\' AS-derived Neoantigens')
plt.hlines(6,0-1,12,linestyle='dashed')
plt.text(x=13,y=6-0.5,s='ASNEO:min#:6',fontsize=8,color='r')
plt.hlines(194,0-1,12,linestyle='dashed')
plt.text(x=13,y=194-0.5,s='ASNEO:max#:194',fontsize=8,color='r')
plt.hlines(69.8,0-1,12,linestyle='dashed')
plt.text(x=13,y=69.8-0.5,s='ASNEO:mean#:69.8',fontsize=8,color='r')
plt.xlabel('Patients#')
plt.ylabel('AS-derived Neoantigens#')    
           
plt.savefig('/Users/ligk2e/Desktop/Ovarian_result.pdf',bbox_inches='tight')
