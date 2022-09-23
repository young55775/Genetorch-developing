import os
import sys


def model(gene):
    lst = os.listdir("apf")
    if gene+'.apf' in lst:
        return "apf\\{}.apf".format(gene)
    else:
        print("file not exist, please refer to https://alphafold.com/",file=sys.stderr)
