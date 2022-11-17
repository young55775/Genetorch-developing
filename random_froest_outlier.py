import pandas as pd

# [15072433, 30351854, 44135655, 61629484, 82553664, 100272606, 100286400]
# 输入在全基因组模型中的碱基位置，返回基因名或者不在基因上。
# 若不在基因上则返回染色体的位置以及上下游基因名
# 计算dumpy的mutation factor并作散点图

chrom = {}


class OutLier:
    def __init__(self, pos):
        self.pos = pos
        self.chrom = str()
        self.property = str()
        self.relative_pos = 0
        self.get_chrom()

    def get_chrom(self):  # 取得该点的染色体位置
        if 0 <= self.pos < 15072433:
            self.chrom = 'I'
            self.relative_pos = self.pos
        if 15072433 <= self.pos < 30351854:
            self.chrom = 'II'
            self.relative_pos = self.pos - 15072433
        if 30351854 <= self.pos < 44135655:
            self.chrom = 'III'
            self.relative_pos = self.pos - 30351854
        if 44135655 <= self.pos < 61629484:
            self.chrom = 'IV'
            self.relative_pos = self.pos - 44135655
        if 61629484 <= self.pos < 82553664:
            self.chrom = 'V'
            self.relative_pos = self.pos - 61629484
        if 82553664 <= self.pos < 100272606:
            self.chrom = 'X'
            self.relative_pos = self.pos - 82553664


def find_gene(OutLier, range_path):
    info = pd.read_csv(range_path, sep='\t', names=['name', 'chrom', 'start', 'end'])
    info = info[info['chrom'] == OutLier.chrom]
    start = info['start'].to_list()
    end = info['end'].to_list()
    name = info['name'].to_list()
    i = 0
    j = 0
    for i in range(len(start)):
        if start[i] > OutLier.relative_pos:
            break
    for j in range(len(end)):
        if end[j] >= OutLier.relative_pos:
            break
    i -= 1
    if i == j:
        return name[i]
    else:
        return name[i], name[j]


def find_outlier(file_path, range_path):
    out = pd.read_csv(file_path)
    pos = out['out'].to_list()
    res = []
    gene = []
    for i in pos:
        gene.append(find_gene(OutLier(i), range_path))
    return gene

