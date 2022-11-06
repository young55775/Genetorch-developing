import pandas as pd
import os
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('TkAgg')
file = os.listdir('/home/data3/cyp/Random_forest')
all = pd.DataFrame()
for i in file:
    n = pd.read_csv('/home/data3/cyp/Random_forest/{}'.format(i))
    all = pd.concat([n, all], axis=1)
train = all.drop(['mutation_factor','code','HCP3_EEMB','HCP3_LTEMB','LEM2_MXEMB'], axis=1)
X_train, X_test, y_train, y_test = train_test_split(train, all["mutation_factor"], test_size=0.2, random_state=1)
rfr = RandomForestRegressor(n_estimators=30,n_jobs=50)
rfr.fit(X_train, y_train)
score_train = rfr.score(X_train, y_train)
score_test = rfr.score(X_test, y_test)