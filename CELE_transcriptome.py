import pandas as pd
import io, os

with open('WBcel235_annotation.gff', 'r') as f:
    lines = [l for l in f if not l.startswith('#')]
    anno = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

anno.columns = ['name', 'source', 'region', 'start', 'end', 'sig', 'chain', 'position', 'info']
info = anno['info'].to_list()
info = list(set(info))
k = []
v = []
for i in info:
    if 'gene=' in i and 'WormBase:WBGene' in i:
        line = i.split(';')
        for ch in line:
            if 'WormBase:WBGene' in ch:
                key = ch.split('WormBase:')[1].split(',')[0]
                k.append(key)
            elif 'gene=' in ch:
                val = ch.split('=')[1]
                v.append(val)

key = list(set(k))
val = list(set(v))
wormbase_id = {}
for i in range(len(k)):
    wormbase_id[k[i]] = v[i]

with open('WBcel235_cds.fa', 'r') as f:
    record = {}
    erro_list = []
    g = f.readlines()
    for i in g:
        if i[0] == '>':
            key = i.split('gene=')[1].split('\n')[0]
            if key in wormbase_id.keys():
                name = wormbase_id[key]
            else:
                name = key
                erro_list.append(key)
            transcript = i.split('>')[1].split(' ')[0]
            if name not in record.keys():
                record[name] = {}
            record[name][transcript] = ''
        else:
            record[name][transcript] += i.split('\n')[0].upper()

code = {'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': '*', 'UAG': '*',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': '*', 'UGG': 'Trp',
        'CGU': 'Arg', 'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}

record_dict = {}


def orf(seq):
    seq_list = []
    # start_ind = detect_start(seq)
    # end_ind = detect_stop(start_ind, seq)
    # if detect_start(seq) == 'none' or detect_stop(seq) == 'none' or (end_ind - start_ind + 1) % 3 != 0:
    #     return 'not good'
    # else:
    for i in range(0, len(seq), 3):
        seq_list.append(''.join([seq[i], seq[i + 1], seq[i + 2]]).replace('T', 'U'))
    return seq_list


def codon_usage(transcript, record_dict, codon_dict):
    for k, v in transcript.items():
        for key, val in v.items():
            codon = orf(val)
            if codon == 'not good':
                continue
            else:
                for i in codon:
                    protein = codon_dict[i]
                    if protein not in record_dict.keys():
                        record_dict[protein] = {}

                    if i not in record_dict[protein].keys():
                        record_dict[protein][i] = 1
                    else:
                        record_dict[protein][i] += 1
    return record_dict


record_dict = codon_usage(record, {}, code)
with open('codon_usage', 'w') as f:
    for k, v in record_dict.items():
        f.write(k + '\n')
        for key, val in v.items():
            f.write(key + '\t' + str(val) + '\n')


def chan(string, index, tar):
    strm = list(string)
    strm[index] = tar
    return ''.join(strm)


def mutation(codon):
    result = []
    c, g = count(codon)
    if c != []:
        for i in c:
            seq = codon[:]
            seq = chan(seq, i, 'U')
            result.append(seq)
    if g != []:
        for i in g:
            seq = codon[:]
            seq = chan(seq, i, 'A')
            result.append(seq)
    return result


def count(seq):
    g = []
    c = []
    for i in enumerate(seq):
        if i[1] == 'G':
            g.append(i[0])
        if i[1] == 'C':
            c.append(i[0])
    return c, g


with open('codon_mutation', 'w') as f:
    for k, v in record_dict.items():
        f.write('>' + k + '\n')
        for key, val in v.items():
            lst = mutation(key)
            mut = ['/'.join([n, code[n]]) for n in lst]
            f.write(key + '\t' + str(val) + '\t' + ';'.join(mut) + '\n')

mut_dict = {}
for k, v in record_dict.items():
    for key, val in v.items():
        lst = mutation(key)
        n = len(lst)
        for i in lst:
            mut = code[i]
            info = k + '->' + mut
            if info not in mut_dict.keys():
                mut_dict[info] = val/n
            else:
                mut_dict[info] += val/n

import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def plot_dict(dict):
    dict_sort = sorted(dict.items(), key=lambda item: item[1], reverse=True)
    data = list(zip(*dict_sort))
    label = data[0]
    val = list(data[1])
    sumer = sum(val)
    val = [n / sumer for n in val]
    x = plt.bar(range(len(val)), val, tick_label=label)
    plt.xticks(rotation=75)
    plt.show()

plot_dict(mut_dict)

def slide(seq):
    record = []
    for i in range(len(seq)-2):
        record.append(''.join([seq[i],seq[i+1],seq[i+2]]))
    return record

record = ''
with open(r'D:\annotation\Genetorch-developing\WBcel235_rna .fna','r') as f:
    file = f.readlines()
    for i in file:
        if i[0] == '>':
            record += '\n'
        else:
            record += i.split('\n')[0]

record = record.split('\n')
del record[0]
code = []
for i in record:
    code.extend(slide(i))
