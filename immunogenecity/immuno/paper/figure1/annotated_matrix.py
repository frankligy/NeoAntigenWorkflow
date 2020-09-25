import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Bio.SubsMat import MatrixInfo


def dict_to_matrix():
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    dic = MatrixInfo.blosum50
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

if __name__ == '__main__':

    rows = list('ARNDCQEGHILKMFPSTWYVX')
    columns = list('ARNDCQEGHILKMFPSTWYVX')
    values = dict_to_matrix().astype(np.int)

    fig, ax = plt.subplots()
    im = ax.imshow(values)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(columns)))
    ax.set_yticks(np.arange(len(rows)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(columns)
    ax.set_yticklabels(rows)


    # Loop over data dimensions and create text annotations.
    for i in range(len(rows)):
        for j in range(len(columns)):
            text = ax.text(j, i, values[i, j],
                           ha="center", va="center", color="w")
    plt.savefig('paper/figure1/blosum50.pdf')


    ###### aaindex1
    import pickle
    with open('paper/figure1/after_pca_matrix_aaindex.p','rb') as f:
        after_pca = pickle.load(f)

    after_pca = np.around(after_pca,decimals=0).astype(np.int)
    rows = list('ARNDCQEGHILKMFPSTWYVX')
    columns = list(['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12'])

    values = after_pca

    fig, ax = plt.subplots()
    im = ax.imshow(values)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(columns)))
    ax.set_yticks(np.arange(len(rows)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(columns)
    ax.set_yticklabels(rows)

    plt.setp(ax.get_xticklabels(),rotation=45,ha='right',rotation_mode='anchor')
    # Loop over data dimensions and create text annotations.
    for i in range(len(rows)):
        for j in range(len(columns)):
            text = ax.text(j, i, values[i, j],
                           ha="center", va="center", color="w")
    fig.tight_layout()
    plt.savefig('paper/figure1/aaindex_pca.pdf')




