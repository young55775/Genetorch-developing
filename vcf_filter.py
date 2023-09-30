import pandas as pd
import io
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file',help = 'raw.vcf file')
parser.add_argument('--af',help = 'allele frequency')
parser.add_argument('--dp',help = 'sequence depth')
args = parser.parse_args()
def read_vcf(path_a):
    with open(path_a, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
anno = read_vcf(args.file)
anno_filter = anno[anno['ALT']!='.']

anno_fin = pd.DataFrame()
for i,v in tqdm(anno_filter.iterrows()):
    val = v['unknown'].split(':')
    dp = val[1]
    ao = val[4]
    if ',' not in ao:
        af = int(ao)/int(dp)
        if af > float(args.af) and int(dp) > int(args.dp):
            anno_fin = pd.concat([anno_fin,pd.DataFrame(v)],axis = 1)
    else:
        continue
anno_fin = anno_fin.T
anno_fin.to_csv('./filtered.vcf',sep='\t',index=False)