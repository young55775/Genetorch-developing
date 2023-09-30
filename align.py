import os
import argparse


def split_name(name):
    elements = name.split('_')
    ind = elements[1]
    num = elements[2].split('.')
    return ind, num[0], num[1]


def main(folder, ref, out):
    file = os.listdir(folder)
    history = []
    for i in file:
        element1 = split_name(i)
        for j in file:
            if i != j:
                element2 = split_name(j)
                if element1[0] == element2[0] and element1[1] == element2[1] and i not in history and  j not in history:
                    history.extend([i, j])
                    os.system("bwa-mem2 mem -t 25 {} {} {} > {}".format(ref,i,j,out))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder',help = 'folder')
    parser.add_argument('--ref',help='ref')
    parser.add_argument('--out',help='out')
    args = parser.parse_args()
    main(args.folder,args.ref,args.out)
