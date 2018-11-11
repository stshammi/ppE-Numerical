import pandas as pd
import numpy as np

# Python code to extract ligo data from .h5 files

with pd.HDFStore('ligo-data.h5', 'r') as d:
    for key in d.keys():
        df=d.get(key)
        df.to_csv(key.replace('/', '') + '.csv')




