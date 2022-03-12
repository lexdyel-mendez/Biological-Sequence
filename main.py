import numpy as np
import pandas as pd

def data_extraction(filepath:str,type='np'):
    """
    This function takes a filepath and returns the data in the file as a numpy array

    :param filepath: the path to the file you want to read in
    :type filepath: str
    :param type: The type of data you want to return, defaults to np (optional)
    :return: a numpy array of the data.
    """
    file = open(filepath)
    data = pd.read_csv(filepath)
    return data.to_numpy()


s1 = "atgct"
s2 = "agct"

def matrixInit(matrix,damping):
    """
    Given a matrix and a damping factor, the function will initialize the first row and first column of the matrix with the
    damping factor

    :param matrix: the matrix to be initialized
    :param damping: The damping factor is the probability that the user will continue clicking
    :return: The matrix with the initial values
    """

    for i in range(matrix.shape[0]):
        matrix[i][0] = i * damping
    for j in range(matrix.shape[1]):
        matrix[0][j] = j * damping

    return matrix


def needlemanWunsch(seqA:str, seqB:str):
    """
    The Needleman-Wunsch algorithm is a dynamic programming algorithm for
    finding the longest path in a matrix

    :param seqA: The first sequence to be compared
    :type seqA: str
    :param seqB: The sequence that is being aligned to the sequence in seqA
    :type seqB: str
    :return: The matrix of scores
    """
    m = len(seqA) + 1       # width of matrix (j)
    n = len(seqB) + 1       # height of matrix  (i)
    d = -2                  # Damping constant

    sc_mtrx = matrixInit(matrix=np.zeros([n,m]) , damping=d)

    for i in range(1,n):
        for j in range(1,m):
            # This is a lambda function that is used to calculate the score of the alignment of the two sequences.
            # If the two characters are the same, then the score is 1, otherwise it is -1.
            c = lambda a,b: 1 if (a == b) else -1

            # This is calculating the score of the alignment of the two sequences. If the two characters are the same,
            # then the score is 1, otherwise it is -1.
            S = c(seqA[j-1], seqB[i-1])


            # A is the sum of the value Diagnaly-right upwards with S
            print(f"Matrix Location [{i-i}][{j-1}]",sc_mtrx[i-i][j-1])
            A = sc_mtrx[i-1][j-1] + S

            # The sum of the
            print(f'Matrix Location [{i}][{j-1}]',sc_mtrx[i][j-1])
            B = sc_mtrx[i][j-1]+d

            print(f'Matrix Location [{i-1}][{j}]', sc_mtrx[i-1][j], '\n')
            C = sc_mtrx[i-1][j]+d

            # Setting the value of the cell in the matrix to the maximum of the three values that are calculated.
            sc_mtrx[i][j] = max(A,B,C)


    print(sc_mtrx)
    return sc_mtrx


if __name__ == '__main__':
    needlemanWunsch(s1,s2)