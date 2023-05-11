import numpy as np

def model(weight_vec, scores_mat):
    weight_vec = np.array(weight_vec)
    scores_mat = np.array(scores_mat)
    scores_mat = scores_mat.transpose()
    output = np.matmul(scores_mat, weight_vec)
    return output

import pandas as pd
d = {'col1': [1, 2], 'col2': [3, 4]}
df = pd.DataFrame(data=d)
data = pd.read_excel("Trade Off.xlsx", sheet_name="General Trade-off table", usecols=(np.arange(0,12,1)), nrows=17)

print(4)