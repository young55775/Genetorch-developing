import pandas as pd
import io
import os
import matplotlib
import matplotlib.pyplot as plt
import genetorch as gt

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
            for m in range(len(i)-3):
                if i[m] == seq[0]:
                    if i[m + 1] == seq[1] and i[m + 2] == seq[2]:
                        if seq not in record.keys():
                            record[seq] = 1
                        else:
                            record[seq] += 1
    return record
