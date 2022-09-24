import os
import pandas as pd
os.mkdir('genetorch/apf')
file = os.listdir('genetorch\\protein')

for i in file:
    try:
        f = pd.read_csv('genetorch\\protein\\{}'.format(i))
        f['x'] = f['x'].round(2)
        f['y'] = f['y'].round(2)
        f['z'] = f['z'].round(2)
        f.to_csv('apf\\{}'.format(i.replace('.pdb', '.apf')),index=False)
    except:
        print(i)