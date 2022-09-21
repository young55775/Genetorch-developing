import numpy as np
import pandas as pd
# import plotly.plotly as py
from plotly.offline import init_notebook_mode, iplot

init_notebook_mode(connected=True)
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import genetorch.pcluster as pc

model = pc.pdb_model('genetorch\\pdb\\tbb-4.pdb')
data = pc.get_model_mat(model)
trace1 = go.Scatter3d(x=data.x, y=data.y, z=data.z, mode='lines', name="tbb-4",
                      marker=dict(color='rgba(128,128,128, 0.4)'), text=data.index)
import genetorch as gt

tbb = gt.reader.readfile(r"E:\g\tbb-4")
gt.finder.find(tbb)
num_lst = pc.get_num(tbb.result, 'tbb-4')
min_num = int(model['num'].to_list()[0])
max_num = int(model['num'].to_list()[-1])
num_lst = [n for n in num_lst if max_num > n > min_num]
import scipy.cluster.hierarchy as htc

aa_list = pc.aa_list(model, num_lst)
matrix = pc.mat(aa_list)
link = htc.linkage(matrix, method='complete', metric='euclidean')
group = htc.fcluster(link, t=25, criterion='distance')
group = list(group)
matrix.insert(loc=matrix.shape[1], column='group', value=group)
data2 = matrix
trace2 = go.Scatter3d(x=data2.x, y=data2.y, z=data2.z, name='mutation', mode='markers', marker=dict(
    size=5,
    color=data2.group,  # set color to an array/list of desired values
    colorscale='Turbo',  # choose a colorscale
    opacity=0.8
), text=data2.index)
fig = go.Figure(dict(data=[trace1, trace2]))
from IPython.display import HTML

HTML(fig.to_html())
