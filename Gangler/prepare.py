import genetorch as gt
import os
import pandas as pd
import numpy as np


def is_syn(description):
    """
    to judge a description is synonymous variation or not
    """
    if description[2:5] == description[-3:]:
        return 0
    elif '.' in description:
        if description.split('.')[1].isdigit():
            return 0
    else:
        return 1


def co_data(taglist):
    """

    :param taglist:
    :return: co-data
    """
    co = pd.DataFrame()
    for i in taglist:
        co = pd.concat([co, i])
    return co


def fill_co(co):
    """
    fill the protein column with other variation types.
    """
    for i in range(co.shape[0]):
        if 'intron' in co.iloc[i, 2] and co.iloc[i, 4] == '':
            co.iloc[i, 4] = 'intron_variation'
        if 'splice' in co.iloc[i, 2] and co.iloc[i, 4] == '':
            co.iloc[i, 4] = 'splice_variation'
        if 'UTR' in co.iloc[i, 2] and co.iloc[i, 4] == '':
            co.iloc[i, 4] = 'UTR_variation'
        elif 'upstream' in co.iloc[i, 2] and co.iloc[i, 4] == '':
            co.iloc[i, 4] = 'delete'
        elif 'downstream' in co.iloc[i, 2] and co.iloc[i, 4] == '':
            co.iloc[i, 4] = 'delete'


def fit(co, threshold):
    """
    same as filter
    """
    return co.groupby(['gene', 'base']).filter(lambda x: len(x) <= threshold)


def get_dict(co):
    """
    taglist -> co_data -> fill_co -> fit -> get_dict
    """
    grouped_m = co.groupby('gene')
    dic = {}
    for n, g in grouped_m:
        lst = g['protein'].to_list()
        lst = [n for n in lst if n != '' and n != 'delete']
        mut = len(lst)  # include syn
        var = [n for n in lst if is_syn(n) != 0 and n != 'intron_variantion' and n != 'UTR_variantion']
        vari = len(var)
        dic[n] = [mut, vari]
    return dic


def get_impact(taglist, threshold=5):
    """
    summary of all above
    """
    co = co_data(taglist)
    fill_co(co)
    f = fit(co, threshold)
    impact = get_dict(f)
    # norm(impact)
    return impact


class NormPool:
    def __init__(self, path, threshold):
        data = gt.reader.readfile(path)
        self.taglist = data.taglist
        self.impact_dict = get_impact(data.taglist, threshold)
        self.result = pd.DataFrame()
        self.co_data = []


def to_genelist(data):
    genelist = []
    for i in data:
        genelist.extend(i.keys())
    return list(set(genelist))


def fill_dic(genelist, data):
    for i in genelist:
        for j in data:
            if i not in j.keys():
                j[i] = [0, 0]


def norm(impact_dic):
    sum = 0
    for v in impact_dic.values():
        sum += v[0]
    for v in impact_dic.values():
        v[0] += 0.1
        v[1] += 0.1
        v[0] /= sum
        v[1] /= sum


def impact_rate(impact_array):
    lst = []
    n = len(impact_array[0])
    for i in range(n):
        lst.append(impact_array[1][i] / impact_array[0][i])
    return lst


def accum(num_lst):
    res = 1
    for i in num_lst:
        res *= i
    return res ** (1 / len(num_lst))


def subtract_gene(data, gene):
    lst = []
    for i in data:
        lst.append(i[gene])
    lst = np.asarray(lst)
    lst = lst.T
    return list(lst)

