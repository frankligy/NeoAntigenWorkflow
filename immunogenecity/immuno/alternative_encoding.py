import pandas as pd
import numpy as np
from Bio.SubsMat import MatrixInfo



def dict_to_matrix():
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    dic = MatrixInfo.blosum40
    matrix = np.zeros([21,21])

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            try:
                matrix[i, j] = dic[(amino[i], amino[j])]
            except KeyError:
                try:
                    matrix[i, j] = dic[(amino[j], amino[i])]
                except:
                    matrix[i, j] = -1
    return matrix


def blosum(peptide):
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = dict_to_matrix()
    encoded = np.empty([len(peptide), 21])  # (seq_len,21)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded
