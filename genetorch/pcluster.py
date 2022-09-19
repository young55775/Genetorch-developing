import pandas as pd
import scipy.cluster.hierarchy as hct
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


class AminoAcid:
    def __init__(self, pro, num):
        self.data = pro[pro['num'] == str(num)]
        self.num = num
        self.name = self.data['aa'].to_list()[0]
        self.x = self.data['x'].to_list()
        self.y = self.data['y'].to_list()
        self.z = self.data['z'].to_list()
        self.x = [float(n) for n in self.x]
        self.y = [float(n) for n in self.y]
        self.z = [float(n) for n in self.z]
        self.median = self.get_grav()

    def get_grav(self):
        n = len(self.x)
        return [sum(self.x) / n, sum(self.y) / n, sum(self.z) / n]


def mat(aa_list):
    matrix = [n.median for n in aa_list]
    ind = [n.name + str(n.num) for n in aa_list]
    return pd.DataFrame(matrix, index=ind, columns=['x', 'y', 'z'])


def cluster_plot(matrix):
    link = hct.linkage(matrix, method='complete', metric='euclidean')
    hct.dendrogram(link, leaf_font_size=10, labels=matrix.index)


def get_num(result, target):
    data = result[result['gene'] == target]['variation'].to_list()[0]
    num_lst = []
    keys = list(data.keys())
    for i in keys:
        if i[-1] != '*':
            num_lst.append("".join(list(filter(str.isdigit, i))))
    while '' in num_lst:
        num_lst.remove('')
    num_lst = [int(n) for n in num_lst]
    return list(set(num_lst))


def pdb_model(model_path):
    model = []
    with open(model_path) as f:
        for i in f:
            if i[0:4] == 'ATOM' or i[0:6] == 'HETATM':
                model.append([i[17:20].strip(), i[22:27].strip(), i[30:38], i[38:46], i[46:54]])

    model = pd.DataFrame(model,
                         columns=['aa', 'num', 'x', 'y', 'z'])
    return model


def aa_list(model, num_lst):
    return [AminoAcid(model, n) for n in num_lst]


def cluster(model_path, result, gene):
    model = pdb_model(model_path)
    min_num = int(model['num'].to_list()[0])
    max_num = int(model['num'].to_list()[-1])
    num_lst = get_num(result, gene)
    num_lst = [str(n) for n in num_lst if max_num > n > min_num]
    aa_lst = aa_list(model, num_lst)
    matrix = mat(aa_lst)
    cluster_plot(matrix)


def get_model_mat(model):
    num_max = model['num'].to_list()[-1]
    num_lst = list(range(1, int(num_max)))
    num_lst = [str(n) for n in num_lst]
    aa = aa_list(model, num_lst)
    return mat(aa)


# PCA
# 3D_plot
import matplotlib.style as mplstyle


def td_plot(model_path, result, gene, dist=25):
    mplstyle.use('fast')
    color = ['orange', 'crimson', 'violet', 'navy', 'y', 'indigo', 'green', 'maroon', 'goldenrod', 'forestgreen',
             'darkslategray', 'darkorange']
    model = pdb_model(model_path)
    min_num = int(model['num'].to_list()[0])
    max_num = int(model['num'].to_list()[-1])
    num_lst = get_num(result, gene)
    num_lst = [str(n) for n in num_lst if max_num > n > min_num]
    aa = aa_list(model, num_lst)
    matrix = mat(aa)
    link = hct.linkage(matrix, method='complete', metric='euclidean')
    group = hct.fcluster(link, t=25, criterion='distance')
    group = list(group)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = matrix['x'].to_list()
    y = matrix['y'].to_list()
    z = matrix['z'].to_list()
    names = [n.name + str(n.num) for n in aa]
    for i in range(len(group)):
        ax.scatter3D(x[i], y[i], z[i], c=color[group[i]], alpha=1)
        ax.text(x[i] + 0.4, y[i] + 0.4, z[i] + 0, names[i], c='k', fontsize=7)
    full_mat = get_model_mat(model)
    xa = full_mat['x'].to_list()
    ya = full_mat['y'].to_list()
    za = full_mat['z'].to_list()
    ax.plot(xs=xa, ys=ya, zs=za, c='dimgray', alpha=0.2)
    ax.set_title('mutation clustering of ' + gene)
    plt.show()



def html_plot(model_path, result, gene, dist=25):
    color = ['orange', 'crimson', 'violet', 'pink', 'y', 'indigo', 'green', 'maroon', 'goldenrod', 'forestgreen',
             'darkslategray', 'darkorange']
    model = pdb_model(model_path)
    min_num = int(model['num'].to_list()[0])
    max_num = int(model['num'].to_list()[-1])
    num_lst = get_num(result, gene)
    num_lst = [str(n) for n in num_lst if max_num > n > min_num]
    aa = aa_list(model, num_lst)
    matrix = mat(aa)
    link = hct.linkage(matrix, method='complete', metric='euclidean')
    group = hct.fcluster(link, t=25, criterion='distance')
    group = list(group)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x = matrix['x'].to_list()
    y = matrix['y'].to_list()
    z = matrix['z'].to_list()
    names = [n.name + str(n.num) for n in aa]
    # for i in range(len(group)):
    #     ax.scatter3D(x[i], y[i], z[i], c=color[group[i]], alpha=1)
    #     ax.text(x[i] + 0.4, y[i] + 0.4, z[i] + 0, names[i], c='k', fontsize=7)
    full_mat = get_model_mat(model)
    xa = full_mat['x'].to_list()
    ya = full_mat['y'].to_list()
    za = full_mat['z'].to_list()
    data = []
    for i in range(len(xa)):
        data.append([xa[i], ya[i], za[i]])
    c = Line3D()
    c.add(series_name=gene, data=data)
    c.set_global_opts(Visualmap_opts = opts.VisualMapOpts(range_color='gray'))
    c.render('cluster of' + gene)
