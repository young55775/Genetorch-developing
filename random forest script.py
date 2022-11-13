import pandas as pd
import os
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import matplotlib
from matplotlib import pyplot as plt
def main():
    matplotlib.use('TkAgg')
    file = os.listdir('/home/data3/cyp/Random_forest')
    all = pd.DataFrame()
    for i in file:
        n = pd.read_csv('/home/data3/cyp/Random_forest/{}'.format(i))
        all = pd.concat([n, all], axis=1)
    train = all.drop(['mutation_factor'], axis=1)
    X_train, X_test, y_train, y_test = train_test_split(train, all["mutation_factor"], test_size=0.2, random_state=23)
    name = train.columns()
    rfr = RandomForestRegressor(n_estimators=64,n_jobs=-1)
    rfr.fit(X_train, y_train)
    score_test = rfr.score(X_test, y_test)
    y_pred = rfr.predict(X_test)
    plt.scatter(y_test,y_pred,s=4,label = score_test)
    plt.title('score {}'.format(n,score_test),loc='right')
    plt.savefig('/home/data3/cyp/64_extimator_with_{}_destroyed.jpg'.format(n))
    plt.close()
# score_train = rfr.score(X_train, y_train)
# score_test = rfr.score(X_test, y_test)

if __name__ == '__main__':
    main()