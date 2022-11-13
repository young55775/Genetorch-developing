import numpy as np
import scipy.stats as stats
import pandas as pd
import io
import os
import matplotlib
import matplotlib.pyplot as plt
import genetorch as gt
import numpy as np
import seaborn as sns
from collections import Counter
import random
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
#数据前处理，将突变平滑化

def norm_unit(mean): #根据每个染色体的间隔取正态分布seed
    x = np.linspace(0,1,2*int(mean)+1)
    g = stats.norm(0.5,0.15)
    seed = g.pdf(x)
    return list(seed)

def mut_spot(co): #返回一个list记录各个碱基的突变次数
    res = {}
    for k in genome.keys():
        if k not in res.keys():
            res[k] = [0] * len(genome[k])
    lst = list(zip(co['chrom'].to_list(), co['pos'].to_list()))
    for i in lst:
        res[i[0]][i[1] - 1] += 1
    return res

def smooth_lst(lst,seed):
    new_list = [0] * len(lst)
    for i in range(len(lst)):
        if lst[i] != 0:
            if i < (len(seed)-1)/2:  # 突变发生在一端
                gap = int((len(seed)-1)/2 - i)
                seed_gap = seed[gap:]
                k = 0
                for j in range(i, i + len(seed)- gap):
                    new_list[j] += seed_gap[k] * lst[i]
                    k += 1
            elif len(lst) - i - 1 < (len(seed)-1)/2: #另一端
                gap = (len(seed)-1)/2 + i + 1 - len(lst)
                k = 0
                for j in range(i - int((len(seed)-1)/2) - 1, len(lst)):
                    new_list[j] += seed[k] * lst[i]
                    k += 1
            elif (len(seed)-1)/2 < i < len(lst) - (len(seed)-1)/2:
                k = 0
                for j in range(i - int((len(seed)-1)/2), i + int((len(seed)-1)/2)-1):
                    new_list[j] += seed[k] * lst[i]
                    k += 1
    return new_list

def cal_dist(lst):
    dist = []
    start = 0
    for i in range(len(lst)):
        if lst[i] != 0:
            dist.append(i-start)
            start = i
    return np.mean(dist)

def ave_lst(lst,scale): #以100为单位作平均值
    new = []
    ori = len(lst)
    for i in range(len(lst)-(scale-1)):
        new.append(np.mean(lst[i:i+scale]))
    while len(new) < ori:
        new.append(0)
    return new
#1.计算间距均值
#2.使用均值，mean=0.5，sigma=0.15画正态曲线叠加

def preprocess(lst): #总流程，不改变参数的情况下调用这个函数就可以
    mean = cal_dist(lst)
    seed = norm_unit(mean)
    res = smooth_lst(lst,seed)
    ave_lst(lst,100)
    return res

#1.读取wig文件，在readwig函数上加以改动
def read_wig_sny(file_path,genome,chrom):
    with open('training_datas\\' + file_path, 'r') as f:
        m = f.readlines()
        data = [0] * len(genome[chrom])
        for k in m:
            if k[0] == chrom:
                i = k.split('\n')[0].split('\t')
                start = int(i[1])
                end = int(i[2])
                for n in range(start, end):
                    data[n] = int(i[3])
    return data

def mean_sny(lst,genome,chrom): #lst内放置成组的sny实验室的文件名，输出平均值
    data = []
    for i in lst:
        data.append(read_wig_sny(i,genome,chrom))
    data = np.asarray(data)
    return np.mean(data,axis = 0)

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

def read_wig_sny_float(file_path,genome,chrom):
    with open('training_datas\\' + file_path, 'r') as f:
        m = f.readlines()
        data = [0] * len(genome[chrom])
        for k in m:
            i = k.split('\n')[0].split('\t')
            if i[0] == chrom:
                start = int(i[1])
                end = int(i[2])
                for n in range(start, end):
                    try:
                        data[n] = float(i[3])
                    except:
                        print(i)
    return data

def mean_sny_float(lst,genome,chrom): #lst内放置成组的sny实验室的文件名，输出平均值
    data = []
    for i in lst:
        data.append(read_wig_sny_float(i,genome,chrom))
    data = np.asarray(data)
    data = list(np.mean(data,axis=0))
    data = ave_lst(data,100)
    return data

def read_wig_single(file_path,genome,chrom,splitkw): #读取只记录一个点信息的wig，并smooth化
    with open('training_datas\\' + file_path, 'r') as f:
        m = f.readlines()
        data = [0] * len(genome[chrom])
        start = False
        for k in m:
            if splitkw in k:
                ch = k.split(splitkw)[1].split(' ')[0]
                if ch == chrom:
                    start = True
                    print('start at chrom' + ch)
                if ch != chrom:
                    start = False
                    print('skip' + ch)
            if start:
                if k[0] != 'v' and k[0] != '#':
                    con = k.split('\n')[0].split('\t')
                    try:
                        data[int(con[0])-1] += float(con[1])
                    except:
                        print(con)
    return preprocess(data)


def mean_single_float(lst,genome,chrom,splitkw): #lst内放置成组的single_wig的文件名，输出平均值
    data = []
    for i in lst:
        data.append(read_wig_single(i,genome,chrom,splitkw))
    data = np.asarray(data)
    return list(np.mean(data,axis = 0))

def save(name,res):
    data = []
    for i in res.values():
        data.extend(i)
    if len(data) == 100272607:
        data = pd.DataFrame(data, columns=[name])
        data.to_csv('full dataset\\{}.csv'.format(name), index=False)
    else:
        print(name)


try:
    res = {}
    for i in chrom:
        res[i] = mean_single_float(['GSM1255289_SDQ0804_HCP3_N2_LTEMB_1_MA2Cscore.wig'],genome,i,'chrom=chr')
    save('HCP3_LTEMB',res)
except:
    print('oh my god')