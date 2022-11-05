import pandas as pd
import os
from sklearn.datasets import load_boston
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
file = os.listdir('MLdataset')
all = pd.DataFrame()
for i in file:
    n = pd.read_csv('MLdataset\\{}'.format(i))
    all = pd.concat([n, all], axis=1)
train = all.drop(['mutation_factor'], axis=1)
X_train, X_test, y_train, y_test = train_test_split(train, whole["mutation_factor"], test_size=0.2, random_state=1)
rfr = RandomForestRegressor(n_estimator=30)
rfr = RandomForestRegressor(n_estimators=30)
rfr.fit(X_train, y_train)
score_train = rfr.score(X_train, y_train)
score_test = rfr.score(X_test, y_test)
rfr = RandomForestRegressor(n_estimators=300)
rfr.fit(X_train, y_train)
score_train2 = rfr.score(X_train, y_train)
score_test2 = rfr.score(X_test, y_test)