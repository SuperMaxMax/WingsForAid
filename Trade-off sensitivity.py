import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_excel("Trade Off.xlsx", sheet_name="General Trade-off table", usecols=(np.arange(0,12,1)), nrows=17)

