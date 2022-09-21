import pandas as pd
import io
import os
import shutil


def read_vcf(path_a):
    with open(path_a, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )


# read info
def simp_file(raw_df):
    b = []
    for i, row in raw_df.iterrows():
        b.append(row['INFO'].split('|'))
    c = []
    for i in range(len(b)):
        c.append([b[i][9], b[i][10], b[i][1], b[i][2], b[i][3], b[i][4], b[i][6]])
    genelist = pd.DataFrame(c, columns=['nucleotide', 'amino acid', 'type', 'impact', 'gene', 'WBGeneID',
                                        'transcript'])
    return genelist


class stockfile:
    def __init__(self, filepath, outpath):
        self.filepath = filepath
        self.filenames = os.listdir(filepath)
        self.temp = filepath + '\\temp'
        self.outpath = outpath
        self.simplist = []
        self.names = []
        self.stock()

    def stock(self):
        os.mkdir(self.outpath)
        file = []
        for i in self.filenames:
            file_s = read_vcf(self.filepath + '\\' + i)
            file_si = simp_file(file_s)
            ind = file_s.drop(['INFO', 'FILTER', 'ID', 'REF', 'ALT', 'QUAL', 'FORMAT', 'unknown'], axis=1)
            co = pd.concat([ind, file_si], axis=1, join='outer')
            co.to_csv(self.outpath + '\\' + i.split('.')[0] + '.csv')


if __name__ == "__main__":
    stockfile(r"E:\g\tbb-4", r"C:\Users\YOUNG\Desktop\tbb-4_stock")
