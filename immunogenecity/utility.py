#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 13:15:33 2020

@author: ligk2e
"""
import numpy as np

properties = {  
        # 'AA':[Hydrophobicity(Kyte-Doolittle),Bulkiness(Zimmerman),'Polarity(Grahtham)']
        # source: TCR contact residue hydrophobicity is a hallmark of immunogenic CD8+ T cell epitopes
        'A':[1.8,11.5,8],    # Alanine              #1
        'C':[2.5,13.46,5.5],  # Cysteine            #2
        'D':[-3.5,11.68,13],   # Aspartic acid      #3
        'E':[-3.5,13.57,12.3], # Glutamic acid      #4
        'F':[2.8,19.8,5.2],    # Phenylalanine      #5
        'G':[-0.4,3.4,9],      # Glycine            #6
        'H':[-3.2,13.69,10.4],  # histidine         #7
        'I':[4.5,21.4,5.2],    # Isoleicine         #8
        'K':[-3.9,15.71,11.3],  # Lysine            #9
        'L':[3.8,21.4,4.9],     # Leucine           #10
        'M':[1.9,16.25,5.7],    # Methionine        #11
        'N':[-3.5,12.82,11.6],  # Asparagine        #12
        'P':[-1.6,17.43,8],     # Proline           #13
        'Q':[-3.5,14.45,10.5],  # Glutamine         #14
        'R':[-4.5,14.28,10.5],  # Arginine          #15
        'S':[-0.8,9.47,9.2],    # Serine            #16
        'T':[-0.7,15.77,8.6],   # Threonine         #17
        'V':[4.2,21.57,5.9],    # Valine            #18
        'W':[-0.9,21.67,5.4],   # Tryptophan        #19
        'Y':[-1.3,18.03,6.2]}   # Tyrosine          #20

def get_hydrophobicity(aa):
    hydrophobicity = 0
    for i in aa:
        a = i.upper()
        hydrophobicity += properties[a][0]
    return hydrophobicity

def KL_divergence(p,q):
    p = np.array(p)
    q = np.array(q)
    pq = np.array([p!=0.0,q!=0.0])
    truth = np.split(pq,len(p),axis=1)  # column-wise
    truth_table = [all(col) for col in truth]
    #print(truth_table)
    return np.sum(np.where(truth_table, p * np.log(p/q), 0))