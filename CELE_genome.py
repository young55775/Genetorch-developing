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


def get_heatmap(path):
    f = readfile(path)
    co = gene_filter(f, 8)
    n = co.shape[0]
    heatmap = region_freq(co, n)
    heat_map_dict(heatmap)
    return heatmap


def prun(sta):
    new = []
    for i in range(len(sta) - 2):
        if not sta[i] == sta[i + 1] == sta[i + 2] == 0:
            new.append(sta[i])
    return new


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
            if start <= len(v):
                at = atcg(v[start:end])
                a = at['A/T'] / sum(list(at.values()))
                a -= 0.6
                split[k].append(a)
            else:
                break
    return split
