import numpy as np

def model(weight_vec, scores_mat):
    scores_mat = scores_mat.transpose()
    output = np.matmul(scores_mat, weight_vec)
    return output
