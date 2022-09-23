# to compare date with mutiple aim
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Gangler import prepare
from matplotlib import pyplot
from collections import Counter
import time


def initiation(a):
    gene = a.result['gene'].tolist()
    variation_num = a.result['variation_number'].tolist()
    # normalization
    # n = a.result[a.result['gene'] == 'ttn-1']['variation_number'].tolist()[0]
    n = sum(variation_num)
    variation_num = [a * 100 / n for a in variation_num]
    pool_file = pd.DataFrame({'gene': gene, 'variation_number': variation_num})
    return pool_file


class snpool:
    def __init__(self, readfile_list, target_list):
        self.target = target_list
        self.initiation_list = []
        self.result = []
        self.readfile_list = readfile_list
        self.data = []
        time1 = time.time()
        for i in self.readfile_list:
            self.data.append(gtfilter(i, 4, rid=[]))
        time2 = time.time()
        print('seperate data' + str(time2-time1))
        self.variation_list = pd.DataFrame()
        self.working()
        time3 = time.time()
        print('get result' + str(time3 - time2))

    def ext_gene(self):
        genelist = []
        for i in self.data:
            genelist.extend(i['gene'].to_list())
        return list(set(genelist))

    def working(self):
        genelist = self.ext_gene()
        data = []
        for i in self.data:
            local_gene = i['gene'].to_list()
            local_var = i['variation_number'].to_list()
            local_tab = dict(zip(local_gene, local_var))
            for j in genelist:
                if not j in local_gene:
                    local_tab[j] = 0
            sum_var = sum(local_tab.values())
            for k, v in local_tab.items():
                local_tab[k] = (v + 0.1) * 100 / sum_var
            data.append(local_tab)
        sum_dict = {}
        for _ in data:
            for k, v in _.items():
                sum_dict.setdefault(k, []).append(v)
        cal_dict = {}
        for k, v in sum_dict.items():
            cal_dict[k] = stat(v)
        self.result = pd.DataFrame.from_dict(cal_dict, orient='index', columns=self.target)
        self.variation_list = pd.DataFrame.from_dict(sum_dict, orient='index', columns=self.target)

    def plot(self, gene):
        x = self.target
        y = list(self.variation_list.loc[gene].values)
        plt.scatter(x, y)
        plt.title('Frequency of ' + gene)


# pick gene and variation_number

def stat(list_ele):
    total = 0
    iter_num = len(list_ele)
    new_list = []
    j = 1
    for i in range(iter_num):
        total = total + list_ele[i]
    for i in range(iter_num):
        if total != 0:
            m = list_ele[i] / total
        else:
            m = 1
        j = j * m
        new_list.append(m)
    alt = []
    for i in list_ele:
        alt.append(i)
    alt.sort()
    sec = alt[-2]
    if max(list_ele) != min(list_ele):
        j = (j / (max(new_list) * (max(list_ele) - sec) ** (len(list_ele) - 1))) ** (1 / (len(list_ele) - 1))
    else:
        j = 'NA'
    m_list = ['_'] * len(new_list)
    if j != 'NA':
        for i in range(len(new_list)):
            if new_list[i] == max(new_list):
                m_list[i] = j / new_list[i]
    return m_list


def get(co_data):
    co_data = co_data.drop_duplicates()
    genelist = list(set(co_data['gene'].to_list()))
    line = []
    for gene in genelist:
        temp_df = co_data[co_data['gene'] == gene]
        sample = list(set(temp_df['tag'].to_list()))
        size = len(sample)
        var = dict(Counter(temp_df['protein'].to_list()))
        var_num = len(var)
        line.append([sample, size, gene, var, var_num])
    result = pd.DataFrame(line, columns=['sample', 'size', 'gene', 'variation', 'variation_number'])
    result = result[~result['variation'].isin([[]])]
    result = result.sort_values(by=['size'], ascending=False)
    result = result.reset_index(drop=True)
    return result


def gtfilter(a, lengthlimit=4, rid=[]):
    prepare.get_impact(a)
    a.candidate = []
    a.suppressor_group = []
    conc = pd.DataFrame()
    for i in range(len(a.taglist)):
        conc = pd.concat([conc, a.taglist[i]])
    filter = conc.groupby(['gene', 'protein']).filter(lambda x: len(x) <= lengthlimit)
    a.co_data = filter
    result = get(filter)
    result = result[~result['gene'].isin(rid)]
    n = result.shape[0]
    result.index = list(range(n))
    a.result = result
    return result


def search(df, col, kw):
    return df[col] == kw
