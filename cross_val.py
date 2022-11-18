#100272607均分为20份
#5000000
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
    factor = []
    score = []
    for i in range(0,100272607,5000000):
        try:
            test = all.iloc[i:i+5000000]
            train = all.drop(list(range(i,i+len(test))))
            y_test = test['mutation_factor']
            X_test = test.drop(['mutation_factor'],axis=1)
            X_train = train.drop(['mutation_factor'],axis=1)
            y_train = train['mutation_factor']
            rfr = RandomForestRegressor(n_estimators=64 , n_jobs=40,max_depth=15,max_features=0.3)
            rfr.fit(X_train,y_train)
            score = rfr.score(X_test,y_test)
            y_pred = rfr.predict(X_test)
            factor.extend(y_pred)
            score.append(score)
        except:
            factor.append('some error in {}'.format(i))
            score.append('some error in {}'.format(i))
    pd.DataFrame(factor,columns=['mutation_factor']).to_csv('/home/data3/cyp/factor_cross.csv',index=False)
    pd.DataFrame(score,columns=['R2']).to_csv('/home/data3/cyp/r2_score.csv',index=False)

