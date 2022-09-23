import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import pandas as pd
import plotly.graph_objs as go
import scipy.cluster.hierarchy as hct


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
    if model_path[-3:] == 'pdb':
        model = []
        with open(model_path) as f:
            for i in f:
                if i[0:4] == 'ATOM' or i[0:6] == 'HETATM':
                    model.append([i[17:20].strip(), i[22:27].strip(), i[30:38], i[38:46], i[46:54]])

        model = pd.DataFrame(model,
                             columns=['aa', 'num', 'x', 'y', 'z'])
        return model


def aa_list(model, num_lst):
    return [AminoAcid(model, n) for n in num_lst if n in model['num'].to_list()]


def cluster(model_path, result, gene):
    num_lst = get_num(result, gene)
    if model_path[-3:] == "pdb":
        model = pdb_model(model_path)
        min_num = int(model['num'].to_list()[0])
        max_num = int(model['num'].to_list()[-1])
        num_lst = [str(n) for n in num_lst if max_num > n > min_num]
        aa_lst = aa_list(model, num_lst)
        matrix = mat(aa_lst)
        cluster_plot(matrix)
    if model_path[-3:] == "apf":
        model = pd.read_csv(model_path)
        min_num = int(model['num'].to_list()[0])
        max_num = int(model['num'].to_list()[-1])
        num_lst = [str(n) for n in num_lst if max_num > n > min_num]
        string = ','.join(num_lst)
        query = "num in [{}]".format(string)
        matrix = model.query(query)
        matrix = matrix.iloc[:, 1:].set_index('aa')
        cluster_plot(matrix)


def get_model_mat(model):
    num_max = model['num'].to_list()[-1]
    num_lst = list(range(1, int(num_max)))
    num_lst = [str(n) for n in num_lst]
    aa = aa_list(model, num_lst)
    return mat(aa)


# PCA
# 3D_plot
def get_link(model_path, result, gene):
    model = pdb_model(model_path)
    min_num = int(model['num'].to_list()[0])
    max_num = int(model['num'].to_list()[-1])
    num_lst = get_num(result, gene)
    num_lst = [str(n) for n in num_lst if max_num > n > min_num]
    aa = aa_list(model, num_lst)
    matrix = mat(aa)
    link = hct.linkage(matrix, method='complete', metric='euclidean')
    return model, link, matrix,aa


def td_plot(model_path, readfile, gene, dist=25):
    if model_path[-3:] == "pdb":
        mplstyle.use('fast')
        color = ['orange', 'crimson', 'violet', 'navy', 'y', 'indigo', 'green', 'maroon', 'goldenrod', 'forestgreen',
                 'darkslategray', 'darkorange']
        model, link, matrix,aa = get_link(model_path, readfile.result, gene)
        group = hct.fcluster(link, t=dist, criterion='distance')
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


# Try to use plotly to generate html form
def html_plot(model_path, readfile, gene, dist=25):
    if model_path[-3:] == "pdb":
        model, link, matrix,aa = get_link(model_path, readfile.result, gene)
        group = hct.fcluster(link, t=25, criterion='distance')
        group = list(group)
        full_mat = get_model_mat(model)
        protein = go.Scatter3d(x=full_mat.x, y=full_mat.y, z=full_mat.z, mode='lines',
                               marker=dict(color='rgba(128,128,128, 0.4)'), text=full_mat.index)
        matrix.insert(loc=matrix.shape[1], column='group', value=group)
        mutation = go.Scatter3d(x=matrix.x, y=matrix.y, z=matrix.z, name='mutation', mode='markers', marker=dict(
            size=5,
            color=matrix.group,
            colorscale='Turbo',
            opacity=0.8
        ), text=matrix.index, textposition="top center")
        fig = go.Figure(dict(data=[protein, mutation],
                             layout=dict(plot_bgcolor='rgba(233,233,233,1)', paper_bgcolor='rgb(233,233,233)',
                                         title='mutation of {}'.format(gene))))
        with open(readfile.path + '\\' + '{}.html'.format(gene), 'w') as f:
            f.write(fig.to_html())
        print("{}.html is in {}".format(gene, readfile.path))
