import pandas as pd
import io
import os
import shutil
import argparse

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

def simp_file(vcf):
    c = []
    for ind,row in vcf.iterrows():
        b = row['INFO'].split('|')
        gene,id,type,base,protein,transcript = b[3],b[4],b[1],b[9],b[10],b[7]
        chrom,pos = row['#CHROM'],row['POS']
        c.append([chrom,pos,gene,id,type,base,protein,transcript])
    return pd.DataFrame(c,columns=['chrom','pos','gene', 'ID', 'type', 'base', 'protein','transcript'])

# def simp_file(raw_df):
#     b = []
#     for i, row in raw_df.iterrows():
#         b.append(row['INFO'].split('|'))
#     c = []
#     for i in range(len(b)):
#         c.append([b[i][3], b[i][4], b[i][1], b[i][9], b[i][10],b[i][7]])
#     genelist = pd.DataFrame(c, columns=['gene', 'ID', 'type', 'base', 'protein','transcript'])
#     return genelist

def main(folder):
    files = os.listdir(folder)
    result_path = folder+'/{}'.format('result')
    if not os.path.isdir(result_path):
        os.mkdir(result_path)
    for i in files:
        vcf = read_vcf(folder+'/{}'.format(i))
        csv = simp_file(vcf)
        csv.to_csv(result_path+'/{}.csv'.format(i.split('.')[0]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder',help='folder containing vcf files')
    args = parser.parse_args()
    main(args.folder)