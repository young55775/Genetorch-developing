from tqdm import tqdm
import pandas as pd
import io
import os
import matplotlib
import matplotlib.pyplot as plt
# import genetorch as gt
import numpy as np
import seaborn as sns
from collections import Counter
import random

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

genome = read_genome('genome')
dist_dict = {}
for chrom in genome.values():
    for j in tqdm(range(4,len(chrom)-4)):
        if chrom[j] not in dist_dict.keys():
            dist_dict[chrom[j]] = {-4:{'A':0,'T':0,'C':0,'G':0},-3:{'A':0,'T':0,'C':0,'G':0},-2:{'A':0,'T':0,'C':0,'G':0},-1:{'A':0,'T':0,'C':0,'G':0},0:{'A':0,'T':0,'C':0,'G':0},1:{'A':0,'T':0,'C':0,'G':0},2:{'A':0,'T':0,'C':0,'G':0},3:{'A':0,'T':0,'C':0,'G':0},4:{'A':0,'T':0,'C':0,'G':0}}
        for i in range(-4,5):
            base = chrom[j+i]
            dist_dict[chrom[j]][i][base] += 1

dist_dict = {}
for chrom in genome.values():
    for j in tqdm(range(4,len(chrom)-4)):
        if chrom[j] not in dist_dict.keys():
            dist_dict[chrom[j]] = {-4:{'A':0,'T':0,'C':0,'G':0},-3:{'A':0,'T':0,'C':0,'G':0},-2:{'A':0,'T':0,'C':0,'G':0},-1:{'A':0,'T':0,'C':0,'G':0},0:{'A':0,'T':0,'C':0,'G':0},1:{'A':0,'T':0,'C':0,'G':0},2:{'A':0,'T':0,'C':0,'G':0},3:{'A':0,'T':0,'C':0,'G':0},4:{'A':0,'T':0,'C':0,'G':0}}
        for i in range(-4,5):
            base = chrom[j+i]
            dist_dict[chrom[j]][i][base] += 1
def plot_dict(n_dict):
    fig,ax = plt.subplots()
    ax.bar(range(9),[n['A']+n['T']+n['C']+n['G'] for n in n_dict.values()], label = 'G',color='#FF7477')
    ax.bar(range(9), [n['A'] + n['T'] + n['C'] for n in n_dict.values()], label='C',color='#E69597')
    ax.bar(range(9), [n['A'] + n['T'] for n in n_dict.values()], label='T',color='#CEB5B7')
    ax.bar(range(9), [n['A']  for n in n_dict.values()], label='A',color='#B5D6D6')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks(range(0,9),'-4 -3 -2 -1 0 1 2 3 4'.split(' '))
    ax.set_xlabel('Position',fontsize=14)
    ax.set_ylabel('Count',fontsize=14)

ems_mutation = pd.read_csv('./EMS_mutation.csv')
mut_dict = {}
mut_dict = {}
for i in chrom:
    mut_dict[i] = ems_mutation[ems_mutation['chr'] == i]['pos'].to_list()

mutation_dict = {}
for chr in chrom:
    for j in mut_dict[chr]:
        if genome[chr][j-1] not in mutation_dict.keys():
            mutation_dict[genome[chr][j-1]] = {-4:{'A':0,'T':0,'C':0,'G':0},-3:{'A':0,'T':0,'C':0,'G':0},-2:{'A':0,'T':0,'C':0,'G':0},-1:{'A':0,'T':0,'C':0,'G':0},0:{'A':0,'T':0,'C':0,'G':0},1:{'A':0,'T':0,'C':0,'G':0},2:{'A':0,'T':0,'C':0,'G':0},3:{'A':0,'T':0,'C':0,'G':0},4:{'A':0,'T':0,'C':0,'G':0}}
        for i in range(-4,5):
            base = genome[chr][j-1+i]
            mutation_dict[genome[chr][j-1]][i][base] += 1

mutation_dict = {}
for chr in chrom:
    for j in mut_dict[chr]:
        if genome[chr][j-1] not in mutation_dict.keys():
            mutation_dict[genome[chr][j-1]] = {-4:{'A':0,'T':0,'C':0,'G':0},-3:{'A':0,'T':0,'C':0,'G':0},-2:{'A':0,'T':0,'C':0,'G':0},-1:{'A':0,'T':0,'C':0,'G':0},0:{'A':0,'T':0,'C':0,'G':0},1:{'A':0,'T':0,'C':0,'G':0},2:{'A':0,'T':0,'C':0,'G':0},3:{'A':0,'T':0,'C':0,'G':0},4:{'A':0,'T':0,'C':0,'G':0}}
        for i in range(-4,5):
            base = genome[chr][j-1+i]
            mutation_dict[genome[chr][j-1]][i][base] += 1

mut_count = {}
for k,v in genome.items():
    mut_count[k] = [0]*len(v)

for k,v in mut_dict.items():
    for i in v:
        mut_count[k][i-1] += 1

mut_dist_dic = {}
for k,v in mut_count.items():
    mut_dist_dic[k] = [sum(v[n:n+100000]) for n in range(0,len(v),100000)]