import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_excel("Trade Off.xlsx", sheet_name="General Trade-off table", usecols=(np.arange(0,12,1)), nrows=17)

scores = data.head(15)
columns = np.array(scores.columns)
for i in range(2, (len(columns))):
    scores[columns[i]] = scores[columns[i]].str.split(";").str[0]


