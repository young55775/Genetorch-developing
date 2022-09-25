import os
import sys


def model(gene):
    current_dir = os.path.dirname(__file__)
    lst = os.listdir(os.path.join(current_dir,"apf"))
    if gene+'.apf' in lst:
        return os.path.join(current_dir, "apf/{}.apf".format(gene))
    else:
        print("file not exist, please refer to https://alphafold.com/",file=sys.stderr)