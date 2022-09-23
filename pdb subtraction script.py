import genetorch.pcluster as pc
import os

os.mkdir('genetorch\\protein')
filelist = os.listdir('genetorch\\pdb')
for i in filelist:
    try:
        model = pc.pdb_model('genetorch\\pdb\\{}'.format(i))
        matirx = pc.get_model_mat(model)
        matirx.to_csv('genetorch\\protein\\{}'.format(i.replace('pdb', 'apf')),index=False)
    except:
        print(i)
        continue
