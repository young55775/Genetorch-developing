import pandas as pd
import os
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

def main():
    matplotlib.use('TkAgg')
    file = os.listdir('/home/data3/cyp/full_data_2')
    all = pd.DataFrame()
    for i in file:
        n = pd.read_csv('/home/data3/cyp/full_data_2/{}'.format(i))
        all = pd.concat([n, all], axis=1)
