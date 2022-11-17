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

comb_5 = ["GCCAC","TACCT","CCCCC","AACAG","TCTCC","CCCTG","ATCTT","AGTGC","TCCTT","TTCTT","CACCC","CCCCA","AGCCT","CGAAA","AATAA","CCACC","GATTT","TTCCA","AACTT","TTCCC","ATTCT","TTGTG","TTCGA","TGCTT","ATGCA","ACCAT","TCCCG","CTCAT","CTTTC","TTTAG","CTCTT","GCCCG","AGCAC","CTCTC","ACCAC","AGCTG","TTCTC","ATCCT","TGTGT","CCCCT","AACCC","CTCCA","TTCCG","TTCAA","TGCAA","ATCAG","ACTGC","TTCCT","ATCGC","TGCTG","GTTGG","TTGCG","GACAA","AAAAA","AAAAT","TCACC","ATTAA","TGCGG","GGTTT","TGGAT","GAGTA","GAGTG","ATGAG","CTGAT","AACAA","ACGCA","TTGTC","TTTCA","TCGAA","CGGAA","GGGGG","CATTG","GAGAA","GGGGT","TAGGT","AGGCA","CAGAG","GCGCC","AGGGG","GGGCG","AAGGT","CTGAA","GCGTG","AAACT","AAGAA","AAGGA","CAGCA","AAGCT","AAGCA","TAGTA","GAGAG","GTGCA","CCGTG","TAGAG","AAGAG","GAGGA","ACGGT","TGATA","TGATT","ATGAT","AAGGG","TGGGT","TAGGC","AGGGA","TCGCC","GTCGA","TCGAT","TTGAA","TCGCA","AGGAG","GCGGT","ACGCT","GTGCG","ATGTG","GGGAC","TGGAA","TTTCC","GCCCT","TTCTG","TGCGT","AACTG","CACAT","TGCAT","CACCG","GACTC","CACGA","CCCAT","TCCAA","CGCCG","GGCAC","CGTAC","GACGA","TTCAT","GCCGG","TCCTC","ATCGA","TGCGA","TACTC","ATCAC","ACCTA","AGCGA","CGCCC","ACCTT","CGCCA","AGCTT","AACAT","TTCTA","AACAC","ACCAA","CTCGC","GACCA","CCCCG","ACCTC","TGCCG","AGCGT","GCCGC","GTCTC","CCCGC","GACTT","TCCGG","TCCAT","ATGTT","GCCGA","CGCGA","CACCA","CTCGA","AGCGG","GCTTC","AACTC","CAAGG","TGGGG","AACCT","CAGAA","ATAGC","GTCGG","ACCTG","AGCAT","GTAGT","GTAGA","AAATG","AGTGG","AGGAA","GTGAA","ATGTA","AAGTT","CTGTT","GCGAA","GGGCA","GTGTA","CCGCA","ATACT","TGGTA","CGGAT","CGGCG","TTGCT","TACAT","GTGAT","TCGGT","TAGCC","CGGGA","CTCGT","CCGGG","GAGAC","TTTGC","TTGAT","AAGAT","GAGCG","CGGTT","CTGCC","AAGTC","ACGAA","ATATT","GGGTA","AGGTA","TCGCT","TTGGA","GAGCA","TAGAT","CGGTA","CCGGA","GCACC","TGGTT","CAGCT","TAGAA","ACGCG","GAGGT","GCGCT","TTGTT","TGGAG","ACGGG","TGGCA","GTGTC","CGGGG","AGGCT","TCGCG","CGGGT","TTGAC","AATCA","TTGAG","CCGAA","AGGCC","CCGCT","TTATG","GTGCT","CAGTA","CCAAG","TGGTG","ATGCT","GGGAG","TCGTT","TCGTA","TTGGG","GCGCA","GAGTT","GGGTG","GATAT","GGAGA","CAGAT","AAATA","AGGCG","ACGTG","ACGTA","ATAGT","TAGTT","TGCAC","AGTCT","ATCCG","AGCAG","GATAA","GCCAG","CCCGT","AGCTC","TGCTA","CACTT","TCCAG","CACAA","CCCTC","TACCG","GCAAC","TGCCT","GACCT","TCCGT","CACTG","GCCAT","GCCCC","CGCCT","ATCTG","ACCCC","ATAAT","TCCGA","CGCAC","ATCGT","CTCCG","TACAA","TGCCA","TAATA","TTCAG","AAACG","TCCCC","AACCA","AGCCA","TCCCT","GTCTG","TCCGC","CGCGT","AGGAT","CTGCG","CGCAT","TATTT","CTCTA","TGACA","AAAGG","CATAC","CATAT","TGCTC","AGTCC","CACGC","GGCAA","TCATC","GACAT","CGATT","GTCCT","TTTAA","TCTGG","CTGGA","CTGGT","GTGGT","GGGAT","TTGCA","AGTTT","TTACT","CATCC","GTTAT","CGGCT","AGGTC","TAGGA","CCCTT","ACGAC","AGGGT","TCGAG","GAGAT","CAGGG","AATTT","TTAGG","GGCCA","GGCGA","TTAAA","GGATG","TTAGT","GTCAG","ATGTC","CTGAG","CTGTG","CAGGT","CGAGA","CAACC","GTGGA","CACTC","TACTG","ATCAT","GACGT","GTAAA","TTCGC","ATCCA","TACGC","TACAG","ATTCC","CGCAA","GAGCT","AAGAC","GGATT","TTTTC","CTCCC","AATGA","AAATT","ATCTC","CCCAG","CGCTC","TTCGT","ACTAT","CGCTT","CGCAG","GCCTG","TTTGA","CCCTA","CACAC","GGCTG","TACCA","ATCAA","CACTA","CGTCA","ACCCG","GCCAA","ACCGG","TGAAA","AGATT","CACCT","GTCCC","GCGCG","TACCC","ACTGA","ACTTT","TGCCC","CTGTA","CCCGA","AATAC","TGAAT","ACGAG","ATCGG","TATCC","GCCTA","AGCGC","CTCAC","TACGG","TTGCC","GGAGG","GCGTA","GAGTC","TTGTA","GCATT","CAGTG","TAAAT","CTGCA","CTGTC","GTGAG","CTGAC","TGAAG","TGGGA","GGAGC","GTGTG","TCTAA","CTACA","CGATA","ACATT","GGGAA","GGGGA","CAGGA","ATGCC","AATCC","CCGGT","ATGAA","CGGCA","CAGTT","GGTCG","CTCTG","GGCAT","AGCCC","CTCCT","ATATA","AATTA","ATCTA","TACTT","GGCTT","ATTGT","ACAAT","CTATA","GGGTT","GTCCA","CTCAG","GTCTA","TTAAG","TAAAA","TTTTT","CTGGG","CTTAA","CTTCC","ATCCC","TCGTC","AATGT","AATGG","TCATA","ACCGT","CGCTA","ATAGG","CTTGA","GTACC","GTATT","AATCG","TCCTG","TGGCC","CAATT","TAGCT","ATGGC","GTGCC","CAGCG","ATTTT","ACACC","AAGGC","AGTGT","AAATC","GAGCC","GGGCT","ACGAT","AAGTA","CCGAT","CAGTC","AATAT","GCATG","TTGGT","ATGCG","TAGAC","CAGAC","AAGTG","CGGTG","TCAAT","TTTGT","GCGGC","CCGCC","ATGAC","CCATA","TAGTG","ATTTG","AGCAA","CTAGG","CATCA","GTCAC","TGTGC","TGTGA","TGACT","CGCGG","TTATC","GCCTT","GGTAT","CGCTG","ACCCT","TGTTG","AATAG","GAAAC","ACCCA","GACTA","GTCAT","GTGTT","GTCTT","GTAGG","TGCGC","CTCAA","ACCAG","ACCGC","CCTCC","ATTCA","TATGC","TATAG","AACGT","TTCAC","TGTAG","CCCAA","GCGAT","GCTGC","TCAAC","TCTTT","GGACC","CATTC","CCAGC","AAAGA","CCGTT","TCGGA","TAGCA","TGGCG","CCATT","AGGTT","TAGGG","GGTGA","GAAGG","TCAAA","AAACA","CGGAG","AACCG","GCATA","AGATC","ATTTA","TTACC","CGTGG","GGCCC","CCAGA","GGCTC","TGCAG","AGCCG","ACTGG","ATTAT","GTTTG","CCCAC","AACGG","GCTTT","TAAGG","TCCCA","ATTAG","CTACC","CGTTA","ATACC","GAAAA","CTATT","TTTAC","TAAAC","CCGAG","GCTCC","GTATC","GAGGC","TCAGA","TTTTG","TTAGC","GGTGC","CTGGC","TCTAG","TATTC","TCTTA","GCAGC","GAACT","AAAAC","CAAAT","TAAGC","AGATA","GTTTA","AACTA","TGTAA","CACAG","AGAAA","TCCAC","CGCGC","TTTGG","TCATT","TTAGA","TGTTA","CCGTA","CGGCC","ATTGG","TTTTA","TCCTA","CACGT","TTTAT","GGCAG","CCTGT","GCGGA","ACGTT","TTGGC","GACAC","CGGAC","GACCC","TATAT","GTCGC","GCCCA","CAGCC","AAAGT","AAGCC","CGAGG","AAAGC","AGGTG","ATGGA","ACGCC","AGGGC","TAATG","TCTAC","AGGAC","CTACT","GCGAG","GTATG","TAACA","AGAAT","TTTCT","AGTAC","GAATA","TCGGG","TTATT","ACAGA","GAATT","GTCGT","GTCCG","CGTGA","TGTCT","TCTCA","AGTCA","AATTC","GAACC","GACTG","AGACA","GGCGT","GGCGC","CATGT","TGAGC","AGACG","GTACT","AGTAA","CCAAT","ATAGA","TCACG","GACGG","GTTCC","TTATA","GTCAA","ACGGA","CAAAA","TCAAG","GATTA","TTACG","TGGTC","AAGCG","GCAGT","CAAAC","GTGAC","CCGAC","GACGC","TGAGT","GGCTA","TACTA","TTCGG","CCCGG","GGAGT","ACAAC","TCTCT","AACGC","CATTT","GCTAC","CAGGC","CGTTT","ACTCC","TGTCC","CTGCT","GAGGG","TAGCG","CGATG","TAATT","AGACT","GTGGC","CCACT","TCAGT","GATCC","TGAGA","GTTTT","CGAAT","GCGGG","GGGGC","ATACA","GATGG","TACAC","GACCG","TCTGT","AATCT","TGTAT","GCTGA","ACAGG","GCAAT","AGAGG","GGGTC","GCGTT","TATTA","TGTTT","ATGGG","TCACA","GGTAA","ATTTC","CAAGT","TATAA","GATAC","CATTA","ACTAA","AGTAT","CTTCG","TGGAC","TCAGG","CGACT","TCTCG","AGACC","ATTGA","TCGTG","AGTTC","CTAGT","AAACC","AATTG","TGTGG","AGCTA","TGAGG","ATTGC","TACGA","TGACC","TACGT","GGTGG","CATAA","GCACT","GTGGG","GAAGA","AAAAG","TCGAC","CGACC","GGACA","CCGTC","GGCCT","CGACG","TCTTG","GAAAT","CTAGA","GGAAT","GACAG","GTTGC","TATCT","AACGA","ACAGT","CTATG","TAACC","GTTGA","CCGCG","CTCGG","TTAAT","TAGTC","CCGGC","GATGA","CGTCG","ACGTC","CCTGG","ATGGT","ATTAC","GTATA","CAATG","ACATA","TATGA","GGCGG","ACACT","ATAAA","CCACA","CGTAT","GGTCA","GGAAG","GCTGT","TGGCT","AGATG","TATCA","CAATA","AGTTG","ACATG","CACGG","GCCTC","TGATG","TGGGC","ACAGC","CCTTC","GCACA","TCGGC","CGGGC","GCCGT","TGTCA","TTACA","ACAAA","CAAAG","CGGTC","TCTAT","GGGCC","CTTCA","GAACG","GAAAG","GATTG","GCAAG","TGTTC","CCTTG","GCAGA","TCTGA","TGACG","GCAAA","ACCGA","GTACA","GAATC","GATCT","AGTGA","TATGT","GCTGG","ACTTA","ATATG","TGTCG","TGTAC","AGTCG","AATGC","TATTG","ATATC","CTTGG","CCTCA","GCGTC","ATAAC","CGTGT","CCTTT","CCTCG","TCACT","AGTAG","CCAGT","TAATC","GGAAC","CAATC","CCTAG","CTAAA","TATGG","GAATG","GTTGT","GGACG","CATAG","AGAAG","CCAAA","CCTAA","TATAC","CAAGA","CTTTG","CGACA","CGATC","GCATC","GTTAA","CCTAC","ACTAG","CCAAC","AGTTA","ACTCA","GGAAA","ACACG","ACGGC","GTTCT","GGCCG","GTAAT","CTTTT","GCTCA","ACTAC","GCGAC","CAAGC","TTAAC","CCATC","GTTTC","TCTGC","ACAAG","CTTAT","CATCG","GAACA","CAACA","CTTTA","TGAAC","TCATG","GCAGG","ATAAG","CGTTG","GGTGT","GGTCC","CCAGG","TATCG","GCTTA","AGAGC","GAAGC","ACTTG","GGTTA","ATTCG","CATGA","GATTC","CGTGC","CCATG","CTTGT","CCTCT","CGAAG","CTATC","GTAAG","GCTCG","GTTCA","CCTAT","GTACG","GGATA","ACACA","GATGC","TCTTC","GGTCT","CGAGC","CATGG","CTAAT","GATCG","ACTCG","GTAAC","CTAAC","GGTTG","CAACG","CTTCT","TAAGA","TAACG","CATCT","GATAG","GCTAA","GTTAG","CGTAG","TAAAG","CTTAC","GGATC","ACTGT","CGTAA","GAAGT","CTAGC","GGTAG","CCACG","GATGT","ACTTC","GGACT","CATGC","GGTTC","TTTCG","CCTGA","GTTCG","GCTCT","CGTCC","GCTAT","GCACG","AGAAC","CTAAG","GTAGC","CAACT","CTACG","CGTCT","TCAGC","GGTAC","GCTTG","CTTGC","CTTAG","ACATC","CCTTA","CGTTC","GATCA","CGAGT","TAACT","TAAGT","GTTAC","CCTGC","TGATC","ACTCT","AGAGT","AGAGA","ATACG","CGAAC"
]
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
    plt.bar(range(len(val)), val, tick_label=label)
    plt.xticks(rotation=75)
    plt.show()




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
def region_freq(co_data, n,block):
    record_dict = {}
    chrom = co_data['chrom'].to_list()
    pos = co_data['pos'].to_list()
    group = list(zip(chrom, pos))
    error = []
    for i in group:
        if i[0] != 'MtDNA':
            if i[0] not in record_dict.keys():
                record_dict[i[0]] = [0] * block
            else:
                try:
                    ind = int(i[1]) // n
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
    filt = conc.groupby(['chrom', 'pos']).filter(lambda x: len(x) <= lengthlimit)
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
    return False


# 比较前四位的指标
def ranking(lst1, lst2, comb):
    gap = 0
    for i in comb:
        if i[2] == 'C' or i[2] == 'G' and search_codon(lst1, i) and search_codon(lst2, i):
            k = abs(search_codon(lst1, i) - search_codon(lst2, i))
            if k > 4:
                gap += k
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


# 计算一个co内每个triple base的突变率
def standard_rate(co, genome,scale):
    chrom = co['chrom'].to_list()  # 计数有多少个triple被突变了
    pos = co['pos'].to_list()
    res = []
    for i in range(len(chrom)):
        code = ''.join([genome[chrom[i]][pos[i] - 2], genome[chrom[i]][pos[i] - 1], genome[chrom[i]][pos[i]]])
        res.append(code)
    count = dict(Counter(res))
    all = []
    for i in genome.keys():  # 计算全基因的triple数
        all.extend(slide(genome[i]))
    count_all = dict(Counter(all))
    re = {}
    for k in count.keys():
        val = count[k] / count_all[k]
        re[k] = val/scale
    return re,count,count_all


# 计算按照突变频率每300000bp内会有多少个什么类型的突变
def standard_mut_block(genome, length, standard_rate, comb):
    split = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    for k, v in genome.items():
        for i in range(int(21000000) // int(length)):
            start = i * length
            end = start + length
            n = len(v)
            if int(start) <= int(n):
                at = v[start:end]
                a = {}
                for i in range(len(at) - 2):
                    code = ''.join([at[i], at[i + 1], at[i + 2]])
                    if code not in a.keys():
                        a[code] = standard_rate[code]
                    else:
                        a[code] += standard_rate[code]
                fill_comb(comb, a)
                split[k].append(a)
            else:
                break
    return split


# 计算一个co_data里按分区会有多少突变
def block_mut_co(co, genome, length, comb):
    res = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    chrom = co['chrom'].to_list()
    pos = co['pos'].to_list()
    for k in res.keys():
        c = genome[k]
        for i in range(int(21000000) // int(length)):
            start = i * length
            end = start + length
            block = {}
            for j in range(len(chrom)):
                if chrom[j] == k and start <= pos[j] < end:
                    code = ''.join([c[pos[j] - 2], c[pos[j] - 1], c[pos[j]]])
                    if code not in block.keys():
                        block[code] = 1
                    else:
                        block[code] += 1
            fill_comb(comb, block)
            res[k].append(block)
    return res


# 一个简单函数，将两个字典中指定的值作方差，放入一个新字典:
def variation(dict1, dict2):
    res = 0
    n = 0
    for k in dict1.keys():
        if k not in dict2.keys():
            dict2[k] = 0
        if k[1] == 'C' or k[1] == 'G':
            res += (dict1[k] - dict2[k]) ** 2
            n += 1
    return (res / n) ** 0.5


# 把预测的和实际获得的作标准差
def predict_var(std, mut):
    res = {}
    for k, v in std.items():
        if k not in res.keys():
            res[k] = []
        for i in range(len(v)):
            res[k].append(variation(v[i], mut[k][i]))
    return res


# 计算五连密码子
def standard_rate_5(co, genome,scale):
    chrom = co['chrom'].to_list()  # 计数有多少个triple被突变了
    pos = co['pos'].to_list()
    res = []
    for i in range(len(chrom)):
        c = genome[chrom[i]]
        p = pos[i]
        code = ''.join([c[p - 3], c[p - 2], c[p - 1], c[p], c[p + 1]])
        res.append(code)
    count = dict(Counter(res))
    all = []
    for i in genome.keys():  # 计算全基因的triple数
        all.extend(slide_5(genome[i]))
    count_all = dict(Counter(all))
    re = {}
    for k in count.keys():
        val = count[k] / count_all[k]
        re[k] = val/scale
    return re


def slide_5(seq):
    record = []
    for i in range(len(seq) - 4):
        record.append(''.join([seq[i], seq[i + 1], seq[i + 2], seq[i + 3], seq[i + 4]]))
    return record


def block_mut_co_5(co, genome, length, comb):
    res = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    chrom = co['chrom'].to_list()
    pos = co['pos'].to_list()
    for k in res.keys():
        c = genome[k]
        for i in range(int(21000000) // int(length)):
            start = i * length
            if start > len(genome[k]):
                break
            end = start + length
            block = {}
            for j in range(len(chrom)):
                if chrom[j] == k and start <= pos[j] < end:
                    code = ''.join([c[pos[j] - 3], c[pos[j] - 2], c[pos[j] - 1], c[pos[j]], c[pos[j] + 1]])
                    if code not in block.keys():
                        block[code] = 1
                    else:
                        block[code] += 1
            fill_comb(comb, block)
            res[k].append(block)
    return res


def standard_mut_block_5(genome, length, standard_rate, comb,scale):
    split = {'I': [], 'II': [], 'III': [], 'IV': [], 'V': [], 'X': [], 'MtDNA': []}
    for k, v in genome.items():
        for i in range(int(21000000) // int(length)):
            start = i * length
            end = start + length
            n = len(v)
            if int(start) <= int(n):
                at = v[start:end]
                a = {}
                for i in range(len(at) - 4):
                    code = ''.join([at[i], at[i + 1], at[i + 2],at[i+3],at[i+4]])
                    if code not in a.keys():
                        a[code] = standard_rate[code]*scale
                    else:
                        a[code] += standard_rate[code]*scale
                fill_comb(comb, a)
                split[k].append(a)
            else:
                break
    return split

def analyze(seq,model): #输入5碱基的突变频率，预测整个基因的突变曲线
    pro = []
    for i in range(2,len(seq)-2):
        code = ''.join([seq[i-2],seq[i-1],seq[i],seq[i+1],seq[i+2]])
        if code in model.keys():
            pro.append(model[code])
        else:
            pro.append(0)
    pro.insert(0,0)
    pro.insert(0,0)
    pro.insert(-1,0)
    pro.insert(-1,0)
    return pro


def get_seq(gene_name,range_path,genome): #获得seq的DNA序列，由于评分是从第三个碱基开始到倒数第三个碱基终止，各往前后多取两个碱基
    with open(range_path,'r') as f: #例如start = 50， end = 60的情况，取genome[47:62],也就是从genome的第48位取到第62位，经过analyze函数后依旧是10的长度
        m = f.readlines()
        for i in m:
            if i.split('\t')[0] == gene_name:
                chrom = i.split('\t')[1]
                start = int(i.split('\t')[2])-3
                end = int(i.split('\t')[3])+1
                seq = genome[chrom][start:end+1]
                return seq
    return False


def get_gene(gene_name,co,range_path): #从co中获取基因的突变信息
    gene = co[co['gene'] == gene_name]
    pos = gene['pos'].to_list()
    with open(range_path, 'r') as f: #获得基因起点
        m = f.readlines()
        for i in m:
            if i.split('\t')[0] == gene_name:
                start = int(i.split('\t')[2])
                end = int(i.split('\t')[3])
                length = end - start + 1
    pos = [int(n) - start for n in pos] #第一位变成0
    plot = [0] * length
    for i in pos:
        if 0< i - 1 < length - 1:
            plot[i] += 1
    return plot

def re(a): #阶乘
    m = 1
    for i in range(1,a+1):
        m *= i
    return m

def combination(all,num): #求组合数
    return re(all)/(re(num)*re(all-num))

def evaluate(probability,count,scale): #评估一个位点的概率
    if len(probability) != len(count):
        raise (ValueError)
    else:
        color = []
        for i in range(len(probability)):
            p = probability[i]
            c = count[i]
            if combination(scale,c)*(p**c)*((1-p)**(scale-c)) > 0.5:
                color.append(0.02)
            else:
                color.append(combination(scale,c)*(p**c)*((1-p)**(scale-c)))
    return color


def norm_color(lst): # 根据cmap需要生成0-1的数据
    maxi = max(lst)
    mini = min(lst)
    multi = 1/(maxi-mini)
    return [(1-n*multi) for n in lst]


def accumulate(pro_list):
    ac = 0
    for i in range(len(pro_list)):
        line = pro_list[i]
        for j in range(len(pro_list)):
            if j != i:
                line *= 1 - pro_list[j]
        ac += line
    return ac

# plt.scatter(x = list(range(48457)),y = dig_1_cm, c = plt.get_cmap('Reds')(norm_color(dig_1_cm)))


# if __name__ == '__main__':
#     tbb = readfile('tbb-4')
#     tbb_co = gene_filter(tbb, 8)
#     genome = read_genome('genome')
#     standard = standard_rate(tbb_co, genome)
#     tbb_mut = block_mut_co(tbb_co, genome, 300000, comb)
#     tbb_std = standard_mut_block(genome, 300000, standard, comb)
#     res = predict_var(tbb_std, tbb_mut)
#     for k in res.keys():
#         while len(res[k]) < 70:
#             res[k].append(0)
#     del res['MtDNA']
#     heat_map_dict(res)

# tbb = readfile('tbb-4')
# tbb_co = gene_filter(tbb, 8)
# genome = read_genome('genome')
# standard = standard_rate_5(tbb_co, genome,437)
# standard['GCTAG'] = 0
# tbb_mut = block_mut_co_5(tbb_co, genome, 300000, comb)
# tbb_std = standard_mut_block_5(genome, 300000, standard, comb)
# res = predict_var(tbb_std, tbb_mut)

#获得预测的图片
# cla_1_seq = get_seq('cla-1','gene_range',genome)
# cla_1_pro = analyze(cla_1_seq,tbb_standard)
# cla_1_count = get_gene('cla-1',tbb_co,'gene_range')
# cla_1_cm = evaluate(cla_1_pro,cla_1_count,437)

def centri(count9):
    res = {}
    for j in ['A', 'T', 'C', 'G']:
        c = {'-4': {}, '-3': {}, '-2': {}, '-1': {}, '0': {}, '1': {}, '2': {}, '3': {}, '4': {}}
        for i in count9.keys():
            if i[4] == j:
                if i[0]:
                    if i[0] not in c['-4'].keys():
                        c['-4'][i[0]] = count9[i]
                    else:
                        c['-4'][i[0]] += count9[i]
                if i[1]:
                    if i[1] not in c['-3'].keys():
                        c['-3'][i[1]] = count9[i]
                    else:
                        c['-3'][i[1]] += count9[i]
                if i[2]:
                    if i[2] not in c['-2'].keys():
                        c['-2'][i[2]] = count9[i]
                    else:
                        c['-2'][i[2]] += count9[i]
                if i[3]:
                    if i[3] not in c['-1'].keys():
                        c['-1'][i[3]] = count9[i]
                    else:
                        c['-1'][i[3]] += count9[i]
                if i[4]:
                    if i[4] not in c['0'].keys():
                        c['0'][i[4]] = count9[i]
                    else:
                        c['0'][i[4]] += count9[i]
                if i[5]:
                    if i[5] not in c['1'].keys():
                        c['1'][i[5]] = count9[i]
                    else:
                        c['1'][i[5]] += count9[i]
                if i[6]:
                    if i[6] not in c['2'].keys():
                        c['2'][i[6]] = count9[i]
                    else:
                        c['2'][i[6]] += count9[i]
                if i[7]:
                    if i[7] not in c['3'].keys():
                        c['3'][i[7]] = count9[i]
                    else:
                        c['3'][i[7]] += count9[i]
                if i[8]:
                    if i[8] not in c['4'].keys():
                        c['4'][i[8]] = count9[i]
                    else:
                        c['4'][i[8]] += count9[i]
        res[j] = c
    return res


def norm_count9(count9):
    for k,v in count9.items():
        total = v['0'][k]
        for a,b in v.items():
            for m,n in b.items():
                count9[k][a][m] = n/total
    return count9

def barh_dict(dict):
    key = ['A','T','C','G']  #前处理，让中间那一列添加别的碱基为0
    for i in key:
        if i not in dict['0'].keys():
            dict['0'][i] = 0
    a = []
    t = []
    c = []
    g = []

    for v in dict.values():
        a.append(v['A'])
        t.append(v['T'])
        c.append(v['C'])
        g.append(v['G'])
    labels = list(dict.keys())
    plt.bar(labels,a,label = 'A')
    plt.bar(labels,t,bottom=a,label = 'T')
    plt.bar(labels,c,bottom=np.asarray(a)+np.asarray(t),label = 'C')
    plt.bar(labels,g,bottom=np.asarray(a)+np.asarray(t)+np.asarray(c),label = 'G')
    plt.legend(loc = 'upper right')
    plt.show()




def standard_rate_9(co, genome,scale):
    chrom = co['chrom'].to_list()  # 计数有多少个triple被突变了
    pos = co['pos'].to_list()
    res = []
    for i in range(len(chrom)):
        c = genome[chrom[i]]
        p = pos[i]
        code = ''.join([c[p-5],c[p-4],c[p - 3], c[p - 2], c[p - 1], c[p], c[p + 1],c[p+2],c[p+3]])
        res.append(code)
    count = dict(Counter(res))
    all = []
    for i in genome.keys():  # 计算全基因的triple数
        all.extend(slide_9(genome[i]))
    count_all = dict(Counter(all))
    re = {}
    for k in count.keys():
        val = count[k] / count_all[k]
        re[k] = val/scale
    return re,count,count_all
def slide_9(seq):
    record = []
    for i in range(len(seq) - 8):
        record.append(''.join([seq[i], seq[i + 1], seq[i + 2],seq[i+3],seq[i+4],seq[i+5],seq[i+6],seq[i+7],seq[i+8]]))
    return record


def predict_gap(range_path,co,genome,model,scale): #看每个基因预测值与实际值的差距
    gene = list(set(co['gene'].to_list()))
    res = {}
    for i in gene:
        count = co[co['gene'] == i].shape[0]
        seq = get_seq(i, range_path, genome)
        if seq:
            pro = analyze(seq, model)
            count = get_gene(i,co,range_path)
            predict = sum(pro)*scale
            c = sum(count)
            if c > 4:
                print((predict-c)/c)
                res[i] = (c-predict)/predict
                print(i + ' done')
    return res


def get_polII(polII_path,genome):  #读取polII均值文件
    with open(polII_path, 'r') as f:
        m = f.readlines()
        res = {'I': [0] * len(genome['I']), 'II': [0] * len(genome['II']), 'III': [0] * len(genome['III']), 'IV': [0] * len(genome['IV']),
               'V': [0] * len(genome['V']), 'X': [0] * len(genome['X'])}
        for i in m:
            if i[0] != 't' and i[0] != '#' and i[0] != 'M':
                n = i.split('\n')[0].split('\t')
                for k in range((int(n[1]) - 1), (int(n[2]) - 1), 1):
                    res[n[0]][k] = n[3]
    for k in res:
        res[k] = [int(n) for n in res[k]]
    return res



def get_wig(wig_path,genome,step):
    error = []
    with open(wig_path, 'r') as f:
        m = f.readlines()
        res = {'I': [0] * len(genome['I']), 'II': [0] * len(genome['II']), 'III': [0] * len(genome['III']), 'IV': [0] * len(genome['IV']),
               'V': [0] * len(genome['V']), 'X': [0] * len(genome['X']),'MtDNA':[0]*len(genome['MtDNA'])}
        for i in m:
            if i[0] != 't' and i[0] != '#' and i[0] != 'M':
                if i[0] == 'v':
                    try:
                        key = i.split('chrom=chr')[1].split(' ')[0]
                    except:
                        key = 'MtDNA'
                else:
                    n = i.split('\n')[0].split('\t')
                    start = int(n[0])-1
                    try:
                        for k in range(start,start+step):
                            res[key][k] = n[1]
                    except:
                        error.append(key)
    del res['MtDNA']
    for k in res:
        res[k] = [float(n) for n in res[k]]
    return res



def heatmap_pol(polII_path,wig_path,genome,length):
    if polII_path:
        res = get_polII(polII_path,genome)
    if wig_path:
        res = get_wig(wig_path,genome,step=10)
    map = {}
    for k in res.keys():
        if k not in map.keys():
            map[k] = []
        while len(res[k]) < 21000000:
            res[k].append(0)
        for i in range(0,len(res[k]),length):
            start = i
            end = i+length
            map[k].append(sum(res[k][start:end])/length)
    heat_map_dict(map)
    return map

