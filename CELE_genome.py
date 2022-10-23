import pandas as pd
import io
import os
import matplotlib
import matplotlib.pyplot as plt
import genetorch as gt
import numpy as np
import seaborn as sns
from collections import Counter

matplotlib.use('TkAgg')


# with open(r"C:\Users\guozh\Desktop\GCF_000002985.6_WBcel235_genomic.fna", 'r') as f:
#     m = f.readlines()
#
# genome = {}
# for i in m:
#     if i[0] == '>':
#         key = i
#         genome[key] = ''
#     else:
#         genome[key] += i.split('\n')[0]
# print('done')


def search(df, col, kw):
    return df[col] == kw


def read_vcf(path_a):
    with open(path_a, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def simp_file(raw_df):
    b = []
    for i, row in raw_df.iterrows():
        b.append(row['INFO'].split('|'))
    c = []
    for i in range(len(b)):
        c.append([b[i][3], b[i][4], b[i][1], b[i][9], b[i][10]])
    genelist = pd.DataFrame(c, columns=['gene', 'ID', 'type', 'base', 'protein'])
    genelist.insert(loc=0, column='chrom', value=raw_df['CHROM'].to_list())
    genelist.insert(loc=1, column='pos', value=raw_df['POS'].to_list())
    return genelist


class readfile:
    def __init__(self, path=None):
        self.path = path
        self.taglist = []
        self.names = []
        self.readfile()
        self.co_data = pd.DataFrame()
        self.result = pd.DataFrame()
        self.candidate = {}
        for i in self.names:
            self.candidate[i] = []

    def readfile(self):
        if self.path is not None:
            files = os.listdir(self.path)
            filelist = []
            filename = []

            # read vcf from file using pandas dataframe

            for file in files:
                if not os.path.isdir(file):
                    if os.path.basename(file).split('.')[-1] == 'vcf':
                        f = read_vcf(self.path + '/' + file)
                        filelist.append(f)
                        filename.append(os.path.basename(file).split('.')[0])
            self.names = filename

            simp_list = []
            for i in range(len(filelist)):
                simp_list.append(simp_file(filelist[i]))
            self.taglist = simp_list
            for i in range(len(self.taglist)):
                self.taglist[i].insert(loc=len(self.taglist[i].columns), column='tag', value=filename[i])


# read genome file
def read_genome(path):
    with open(path, 'r') as f:
        g = f.readlines()
        g = g[0].split('next')
        del g[7]
        genome = {}
        for i in g:
            k, v = i.split('\t')[0], i.split('\t')[1]
            genome[k] = v
    MtDNA = genome['X']
    x = genome['VI']
    genome['X'] = x
    genome['MtDNA'] = MtDNA
    del genome['VI']
    del g
    return genome


def find_codon(dict, chrom, pos):
    if pos == 0:
        raise ValueError('wrong position')
    else:
        chromosome = dict[chrom]
        pos -= 1
        if pos == 0:
            return ''.join([chromosome[pos], chromosome[pos + 1]])
        else:
            return ''.join([chromosome[pos - 1], chromosome[pos], chromosome[pos + 1]])


def get_codon_dict(co_data, genome):
    codon = {}
    pos = co_data['pos'].to_list()
    chrom = co_data['chrom'].to_list()
    gene = list(zip(pos, chrom))
    for i in gene:
        code = find_codon(genome, i[1], i[0])
        if code not in codon.keys():
            codon[code] = 1
        else:
            codon[code] += 1
    return codon


def plot_dict(dict):
    dict_sort = sorted(dict.items(), key=lambda item: item[1], reverse=True)
    data = list(zip(*dict_sort))
    label = data[0]
    val = list(data[1])
    sumer = sum(val)
    val = [n / sumer for n in val]
    plt.bar(range(len(val)), val, tick_label=label)
    plt.xticks(rotation=75)
    plt.show()


def run(path):
    a = readfile(path)
    gt.finder.filter(a, 8)
    codon = get_codon_dict(a.co_data, genome)
    plot_dict(codon)


comb = ['ATG', 'AGT', 'GTA', 'GAT', 'TAG', 'TGA',
        'ACG', 'AGC', 'CAG', 'CGA', 'GAC', 'GCA',
        'TGC', 'TCG', 'GCT', 'GTC', 'CGT', 'CTG',
        'TAC', 'TCA', 'CAT', 'CTA', 'ATC', 'ACT',
        'AAA', 'GGG', 'TTT', 'CCC',
        'ATA', 'ACA', 'AGA',
        'GAG', 'GTG', 'GCG',
        'CGC', 'CAC', 'CTC',
        'TAT', 'TGT', 'TCT',
        'AAT', 'AAC', 'AAG',
        'GGC', 'GGT', 'GGA',
        'TTC', 'TTG', 'TTA',
        'CCA', 'CCT', 'CCG',
        'TAA', 'TCC', 'TGG',
        'ATT', 'ACC', 'AGG',
        'CAA', 'CGG', 'CTT',
        'GCC', 'GAA', 'GTT']


def match(genome, comb):
    record = {}
    for i in list(genome.values()):
        for seq in comb:
            for m in range(len(i) - 3):
                if i[m] == seq[0]:
                    if i[m + 1] == seq[1] and i[m + 2] == seq[2]:
                        if seq not in record.keys():
                            record[seq] = 1
                        else:
                            record[seq] += 1
    return record


# I:15072434    II:15279421 III:13783801    IV:17493829 V:20924180  X:17718942  Mt:13794
def region_freq(co_data, n):
    record_dict = {}
    chrom = co_data['chrom'].to_list()
    pos = co_data['pos'].to_list()
    group = list(zip(chrom, pos))
    error = []
    for i in group:
        if i[0] != 'MtDNA':
            if i[0] not in record_dict.keys():
                record_dict[i[0]] = [0] * 420
            else:
                try:
                    ind = int(i[1]) // 50000
                    record_dict[i[0]][ind] += 1 / n
                except:
                    error.append(i)
    print(len(error))
    return record_dict


def heat_map_dict(dict):
    data = np.asarray(list(dict.values()))
    plt.imshow(data, cmap='Reds', interpolation='nearest')
    plt.show()


def gene_filter(a, lengthlimit=0.6):
    conc = pd.DataFrame()
    for i in range(len(a.taglist)):
        conc = pd.concat([conc, a.taglist[i]])
    filt = conc.groupby(['gene', 'protein', 'ID', 'base']).filter(lambda x: len(x) <= lengthlimit)
    return filt


def get_impact(a):
    newtaglist = []
    for i in range(len(a.taglist)):
        n = a.taglist[i]
        lista = []
        listb = []
        splice_donor = n[n['type'] == 'splice_donor_variant&intron_variant']['protein'].tolist()
        splice_acceptor = n[n['type'] == 'splice_acceptor_variant&intron_variant']['protein'].tolist()
        for i in range(len(splice_donor)):
            lista.append('X_donor')
        for i in range(len(splice_acceptor)):
            listb.append('X_acceptor')
        c = pd.DataFrame(n[n['type'] == 'splice_donor_variant&intron_variant'])
        d = pd.DataFrame(n[n['type'] == 'splice_acceptor_variant&intron_variant'])
        c['protein'] = lista
        d['protein'] = listb
        n = pd.concat([n, c, d])
        indp = search(n, 'protein', '')
        indq = search(n, 'protein', 'nan')
        ind = search(n, 'type', 'synonymous_variant')
        m = n.loc[ind, :]
        b = n.loc[indp, :]
        e = n.loc[indq, :]
        n = pd.concat([m, n, e]).drop_duplicates(keep=False)
        total = pd.concat([n, b]).drop_duplicates(keep=False)
        newtaglist.append(total)
    a.taglist = newtaglist
    return a.taglist


def get_heatmap(path):
    f = readfile(path)
    co = gene_filter(f, 8)
    n = co.shape[0]
    heatmap = region_freq(co, n)
    heat_map_dict(heatmap)
    return heatmap


def chrom_plot(heatmap, legend, key):
    d = []
    for i in range(len(heatmap)):
        data = prun(heatmap[i][key])
        plt.plot(list(range(len(data))), data, linestyle='--', label=legend[i], alpha=0.2)
        d.append(data)
    ave = np.mean(d, axis=0)
    plt.plot(list(range(len(ave))), ave, label="Average")
    plt.legend(loc='upper right')
    plt.title('Chromosome' + ' ' + key)
    plt.show()


def variation(heatmap):
    d = []
    chro = ['I', 'II', 'III', 'IV', 'V', 'X']
    for key in chro:
        for i in heatmap:
            d.append(i[key])
        std = np.std(d, axis=0)
        plt.plot(list(range(210)), std, linestyle='-', label=key)
        d = []
    plt.legend(loc='upper right')
    plt.show()


def his(heatmap, legend, key):
    data = []
    for i in heatmap:
        data.append(prun(i[key]))
    plt.hist(data, 15, label=legend)
    plt.legend(loc='upper right')
    plt.title('Chromosome' + ' ' + key)
    plt.show()


def bar(heatmap, key):
    data = []
    for i in heatmap:
        data.append(i[key])
    ave = np.mean(data, axis=0)
    plt.bar(range(len(ave)), ave)
    plt.show()


def kde(heatmap, legend, key):
    sns.set_style("whitegrid")
    for i in range(len(heatmap)):
        sns.kdeplot(prun(heatmap[i][key]), fill=True, label=legend[i])
    plt.legend(loc='upper right')
    plt.title('Chromosome' + ' ' + key)
    plt.show()


def kde_all(heatmap):
    chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
    for key in chrom:
        d = []
        for i in heatmap:
            d.append(i[key])
        mean = np.mean(d, axis=0)
        sns.kdeplot(prun(mean), fill=True, label=key)
    plt.legend(loc='upper right')
    plt.show()


# find low/high mutation region
def outlier(dict):
    low = {}
    high = {}
    for k, v in dict.items():
        v = prun(v)
        mean = sum(v) / len(v)
        sd = np.std(v)
        if k not in low.keys():
            low[k] = []
        if k not in high.keys():
            high[k] = []
        for i in enumerate(v):
            if i[1] <= mean - sd:
                low[k].append((i[0] * 50000, (i[0] + 1) * 50000))
            if i[1] >= mean + sd:
                high[k].append((i[0] * 50000, (i[0] + 1) * 50000))
    return low, high


def outlier_plot(map):
    chrom = {'I': 6, 'II': 5, 'III': 4, 'IV': 3, 'V': 2, 'X': 1}
    color = ['crimson', 'c', 'm', 'y']
    legend = ['uncoordinated', 'tbb-4', 'tba-5', 'dumpy']
    labeled = []
    fig, ax = plt.subplots()
    for k, y in chrom.items():
        for i in range(len(map)):
            if map[i][k] != []:
                for j in map[i][k]:
                    if legend[i] not in labeled:
                        ax.plot(list(j), [y, y], c=color[i], linewidth=5, label=legend[i])
                        labeled.append(legend[i])
                    else:
                        ax.plot(list(j), [y, y], c=color[i], linewidth=5)
                y -= 0.12
    ax.legend(loc='upper right')
    ax.set_yticks([1, 2, 3, 4, 5, 6])
    ax.set_yticklabels(('X', 'V', 'IV', 'III', 'II', 'I'))
    plt.show()


# find the gene in these regions:
# if a high/low region appear 3 times in different screening, select it!
def selection(out):
    res = {}
    chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
    for key in chrom:
        if key not in res.keys():
            res[key] = []
        for i in out:
            res[key].extend(i[key])
    for k, v in res.items():
        res[k] = [n for n, m in dict(Counter(v)).items() if m >= 3]
    return res


# get the sequence
def acquire_seq(res, genome):
    seq = {}
    for k, v in res.items():
        if k not in seq.keys():
            seq[k] = []
        for i in v:
            seq[k].append(genome[k][i[0]:i[1]])
    return seq


# count seq of all:
def count_seq(seq):
    res = {}
    for k, v in seq.items():
        for i in v:
            for j in range(len(i) - 2):
                code = ''.join([i[j], i[j + 1], i[j + 2]])
                if code not in res.keys():
                    res[code] = 1
                else:
                    res[code] += 1
    return res


# calculate percentage
def norm(dict):
    val = list(dict.values())
    s = sum(val)
    for k, v in dict.items():
        dict[k] = v / s


# A/T C/G percentage of seq
def atcg(seq):
    per = {'A/T': 0, 'C/G': 0}
    for i in seq:
        if i == 'A' or i == 'T':
            per['A/T'] += 1
        else:
            per['C/G'] += 1
    return per


def split_genome(genome, length):
    split = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    for k, v in genome.items():
        for i in range(int(21000000) // int(length)):
            start = i * length
            end = start + length
            n = len(v)
            if int(start) <= int(n):
                at = atcg(v[start:end])
                a = at['A/T'] / sum(list(at.values()))
                a -= 0.6
                split[k].append(a)
            else:
                break
    return split


# classify high and low freq
def classification(freq_dict):
    high = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': []}
    mid = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': []}
    low = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': []}
    value = []
    for k, v in freq_dict.items():
        line = [n[2] for n in v]
        value.extend(line)
    ave = np.mean(value)
    std = np.std(value)
    for k, v in freq_dict.items():
        for i in v:
            if i[2] < (ave - std):
                low[k].append(i)
            elif i[2] > (ave + std):
                high[k].append(i)
            else:
                mid[k].append(i)
    return high, mid, low


import numpy as np


def pca(X, k):  # k is the components you want
    # mean of each feature
    n_samples, n_features = X.shape
    mean = np.array([np.mean(X[:, i]) for i in range(n_features)])
    # normalization
    norm_X = X - mean
    # scatter matrix
    scatter_matrix = np.dot(np.transpose(norm_X), norm_X)
    # Calculate the eigenvectors and eigenvalues
    eig_val, eig_vec = np.linalg.eig(scatter_matrix)
    eig_pairs = [(np.abs(eig_val[i]), eig_vec[:, i]) for i in range(n_features)]
    # sort eig_vec based on eig_val from highest to lowest
    eig_pairs.sort(reverse=True)
    # select the top k eig_vec
    feature = np.array([ele[1] for ele in eig_pairs[:k]])
    # get new data
    data = np.dot(norm_X, np.transpose(feature))
    return data


# poly和count数序列中的连续片段
def poly(seq):
    dic = {'AAA': 0, 'CCC': 0}
    for k in dic.keys():
        dic[k] = count_continue(k, seq)
    return dic


def count_continue(pattern, seq):
    init = True
    count = 0
    res = 0
    for i in range(len(seq) - 2):
        match = ''.join([seq[i], seq[i + 1], seq[i + 2]])
        if match == pattern:
            count += 1
            init = False
        if match != pattern:
            if init == False:
                count += 2
                init = True
                if count >= 10:
                    res += count
                    count = 0
                else:
                    count = 0
    return res


# 把所有的AT变成A，GC变成G试一试
def transform_seq(seq):
    new_seq = []
    for i in seq:
        if i == 'A' or i == 'T':
            new_seq.append('A')
        else:
            new_seq.append('C')
    return ''.join(new_seq)


def fill(length, map):
    for i in map:
        while len(i) < length:
            i.append(0)


# 在一个dict中选取前28位d
def select(dic, rank):
    trans_dict = list(zip(list(dic.keys()), list(dic.values())))
    trans_dict.sort(key=lambda x: x[1], reverse=True)
    return (trans_dict[:rank])


# 比较两组的方差
def compare(s1, s2):
    k1 = [n[0] for n in s1]
    k2 = [n[0] for n in s2]
    k3 = []
    k3.extend(k1)
    k3.extend(k2)
    kf = list(set(k3))
    std = 0
    for i in kf:
        if i in k1:
            a = s1[k1.index(i)][1]
        else:
            a = 0
        if i in k2:
            b = s2[k2.index(i)][1]
        else:
            b = 0
        std += abs(a - b)
    return std


def impact_heatmap(lst, length, standard):
    genome = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': []}
    ind = np.arange(0, 71 * 300000, length)
    for i in genome.keys():
        for j in range(len(ind) - 1):
            block = {}
            start = ind[j]
            end = ind[j + 1]
            for k in lst:
                if k[0] == i and start <= k[1] < end:
                    key = k[2]
                    if key not in block.keys():
                        block[key] = 1
                    else:
                        block[key] += 1
            if block != {}:
                norm(block)
                s1 = select(block, 28)
                genome[i].append(compare(s1, standard))
            else:
                genome[i].append(0)
    return genome


def process(co):
    chrom = co['chrom'].to_list()
    pos = co['pos'].to_list()
    protein = co['protein'].to_list
    protein = co['protein'].to_list()
    protein = [transform(n) for n in protein]
    return list(zip(chrom, pos, protein))


def transform(description):
    if 'p.' in description:
        p = description.split('.')[1]
        if len(p) >= 7 and p[-1] != '*' and p[-2:] != 'fs':
            m = p[0:3] + '->' + p[-3:]
            return m
        elif p[-1] == '*':
            m = p[0:3] + '->*'
            return m
        elif p[-2:] == 'fs':
            m = p[0:3] + '->fs'
            return m
        else:
            return p
    else:
        return description


def freq_data(co_data):
    protein = co_data['protein'].to_list()
    protein = [transform(n) for n in protein]
    res = dict(Counter(protein))
    return res


def codon_heatmap(co, genome, length):
    lst = list(zip(co['chrom'].to_list(), co['pos'].to_list()))
    res = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': []}
    ind = np.arange(0, 71 * 300000, length)
    for i in res.keys():
        print(i)
        for j in range(len(ind) - 1):
            start = ind[j]
            end = ind[j + 1]
            block = {}
            for k in lst:
                chrom = k[0]
                pos = k[1]
                if chrom == i and start <= pos < end:
                    codon = ''.join([genome[chrom][pos - 2], genome[chrom][pos - 1], genome[chrom][pos]])
                    if codon not in block.keys():
                        block[codon] = 1
                    else:
                        block[codon] += 1
            norm(block)
            res[i].append(block)
    return res


# 把字典里比如TCC GGA求和作为一个数
# 1.三联密码子转换函数
def matching(codon):
    td = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = list(codon)
    seq = seq[::-1]
    return ''.join([td[seq[0]], td[seq[1]], td[seq[2]]])


# 2.输入字典转换为和
def combine(dic):
    new = {}
    for k, v in dic.items():
        m = matching(k)
        if m in dic.keys():
            a = v + dic[m]
        else:
            a = v + 0
        new[k] = a
    return new


# 2.5:填空
def fill_comb(comb, dic):
    for i in comb:
        if i not in dic.keys():
            dic[i] = 0


# 3.取得各key的排名
def sort_dict(dic):
    sort = sorted(dic.items(), key=lambda item: item[1], reverse=True)
    rank = 1
    rank_list = []
    for i in sort:
        i = list(i)
        i.append(rank)
        rank += 1
        rank_list.append(i)
    return rank_list


# 在list（codon,rate,rank）中输入codon返回rank
def search_codon(lst, codon):
    for i in lst:
        if i[0] == codon:
            return i[2]


# 比较前四位的指标
def ranking(lst1, lst2):
    key = ['TCT', 'AGA', 'GGA', 'TTC']
    gap = 0
    for i in key:
        gap += abs(search_codon(lst1, i) - search_codon(lst2, i))
    return gap


# 4.获取全局的排序
def rank_all(co, genome):
    lst = list(zip(co['chrom'].to_list(), co['pos'].to_list()))
    all_dic = {}
    for i in lst:
        chrom = genome[i[0]]
        codon = ''.join([chrom[i[1] - 2], chrom[i[1] - 1], chrom[i[1]]])
        if codon not in all_dic.keys():
            all_dic[codon] = 1
        else:
            all_dic[codon] += 1
    return all_dic


def prun(sta):
    i = len(sta) - 1
    while sta[i] == 0:
        del sta[i]
        i -= 1
    return sta


# 5.比较分区和全局
def ranking_compare(hm, std, comb):
    heatmap = {}
    for k, v in hm.items():
        if k not in heatmap.keys():
            heatmap[k] = []
        for i in v:
            if i != {}:
                fill_comb(comb, i)
                i_new = combine(i)
                i_new = sort_dict(i_new)
                val = ranking(i_new, std)
            else:
                val = 0
            heatmap[k].append(val)
    return heatmap


def slide(seq):
    record = []
    for i in range(len(seq) - 2):
        record.append(''.join([seq[i], seq[i + 1], seq[i + 2]]))
    return record


#计算按照突变频率每300000bp内会有多少个什么类型的突变
def standard_mut_block(genome,length,standard_rate,comb):
    split = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    for k, v in genome.items():
        for i in range(int(21000000) // int(length)):
            start = i * length
            end = start + length
            n = len(v)
            if int(start) <= int(n):
                at = v[start:end]
                a = {}
                for i in range(len(at)-2):
                    code = ''.join([at[i],at[i+1],at[i+2]])
                    if code not in a.keys():
                        a[code] = standard_rate[code]
                    else:
                        a[code] += standard_rate[code]
                fill_comb(comb,a)
                split[k].append(a)
            else:
                break
    return split

#计算一个co_data里按分区会有多少突变
def block_mut_co(co,genome,length,comb):
    res = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    chrom = co['chrom'].to_list()
    pos = co['pos'].to_list()
    for k in res.keys():
        c = genome[k]
        for i in range(int(21000000) // int(length)):
            start = i*length
            end = start + length
            block = {}
            for j in range(len(chrom)):
                if chrom[j] == k and start<=pos[j]<end:
                    code = ''.join([c[pos[j]-2],c[pos[j]-1],c[pos[j]]])
                    if code not in block.keys():
                        block[code] = 1
                    else:
                        block[code] += 1
            fill_comb(comb,block)
            res[k].append(block)
    return res

#一个简单函数，将两个字典中指定的值作方差，放入一个新字典:

compare_list = ['CCC','GGG','CGG','CCG','CGC','GCG','GGC','GCC','TGG','GGT','TCC','CCT','AGG','GGA','ACC','CCA','TCT','TGT','ACA'
                'AGA',]
def variation(dict1,dict2,compare_list):
    res = 0
    for k in compare_list:
        res += (dict1[k] - dict2[k])**2
    return (res/len(compare_list))**0.5