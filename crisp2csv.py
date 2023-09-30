import pandas as pd
import io
import os
import shutil
import argparse



def read_vcf(path_a):
    with open(path_a, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def simp_file(vcf):
    mix_lst = ['mixA', 'mixB', 'mixC', 'mixD', 'mixE']
    c = []
    for ind, row in vcf.iterrows():
        b = row['INFO'].split('|')
        chrom, pos = row['CHROM'], row['POS']
        gene, id, type, base, protein, transcript = b[3], b[4], b[1], b[9], b[10], b[7]
        mix_row = []
        mix = [row['mixA'], row['mixB'], row['mixC'], row['mixD'], row['mixE']]
        print('row:{}'.format(ind))
        for i in range(len(mix)):
            # print(mix[i].split(':')[0])
            if mix[i].split(':')[0] != '.' and int(mix[i].split(':')[0].split(',')[0]) > 0:
                mix_row.append(mix_lst[i])
        mix_res = ':'.join(mix_row)
        c.append([chrom, pos, gene, id, type, base, protein, transcript, mix_res])
    return pd.DataFrame(c,columns = ['chrom','pos','gene', 'ID', 'type', 'base', 'protein','transcript','mix'])
# def simp_file(raw_df):
#     b = []
#     for i, row in raw_df.iterrows():
#         b.append(row['INFO'].split('|'))
#     c = []
#     for i in range(len(b)):
#         c.append([b[i][3], b[i][4], b[i][1], b[i][9], b[i][10],b[i][7]])
#     genelist = pd.DataFrame(c, columns=['gene', 'ID', 'type', 'base', 'protein','transcript'])
#     return genelist



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file',help='vcf file')
    args = parser.parse_args()
    vcf = read_vcf(args.file)
    simp = simp_file(vcf)
    simp.to_csv(args.file.split('.vcf')[0]+'.csv',index=False)