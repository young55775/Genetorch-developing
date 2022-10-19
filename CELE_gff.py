def gff_key(path,keyword):
    with open(path) as f:
        m = f.readlines()
        anno = {}
        keys = list(anno.keys())
        num = 0
        for i in m:
            if i[0] != '#':
                line = i.split('\t')
                if line[0] not in anno.keys():
                    anno[line[0]] = []
                if line[2] == keyword:
                    anno[line[0]].append((line[3], line[4]))
    return anno



def exon_freq(anno,length):
    heatmap = []
    res = {}
    for k, v in anno.items():
        line = []
        i = length
        count = 0
        for j in v:
            if int(j[1]) <= i and j != v[-1]:
                count += (int(j[1]) - int(j[0]))
            else:
                line.append(count)
                count = 0
                i += length
        heatmap.append(line)
    key = ['I','II','III','IV','V','X','MtDNA']
    for i in range(len(heatmap)):
        res[key[i]] = heatmap[i]
    return res

