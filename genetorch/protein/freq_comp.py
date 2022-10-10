from collections import Counter
from matplotlib import pyplot as plt
import matplotlib


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


def plot_freq(co_data):
    matplotlib.use('TkAgg')
    res = freq_data(co_data)
    res_sort = sorted(res.items(), key=lambda items: items[1], reverse=True)
    res_sum = [n[1] for n in res_sort]
    res_sum = sum(res_sum)
    lst = res_sort[:50]
    val = [n[1] / res_sum for n in lst]
    label = [n[0] for n in lst]
    plt.bar(range(len(val)), val, tick_label=label)
    plt.xticks(rotation=75, fontsize=9)
    plt.show()


def mark(readfile):
    sample = len(readfile.taglist)
    snp = readfile.co_data.shape[0]
    return snp / sample
