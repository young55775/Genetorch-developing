import matplotlib
from collections import Counter
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

record = ''
with open(r'D:\annotation\Genetorch-developing\WBcel235_rna .fna', 'r') as f:
    file = f.readlines()
    for i in file:
        if i[0] == '>':
            record += '\n'
        else:
            record += i.split('\n')[0]

record = record.split('\n')
del record[0]

def slide(seq):
    record = []
    for i in range(len(seq) - 2):
        record.append(''.join([seq[i], seq[i + 1], seq[i + 2]]))
    return record


code = []
for i in record:
    code.extend(slide(i))



graph = dict(Counter(code))


def plot_dict(dict):
    dict_sort = sorted(dict.items(), key=lambda item: item[1], reverse=True)
    data = list(zip(*dict_sort))
    label = data[0]
    val = list(data[1])
    sumer = sum(val)
    val = [n / sumer for n in val]
    plt.bar(range(len(val)), val, tick_label=label)
    plt.xticks(rotation=75)
    plt.show()


# plot_dict(graph)


def exchange(seq):
    s = list(seq[:])
    a = s[0]
    s[0] = s[2]
    s[2] = a
    return ''.join(s)

graph2 = {}
for k, v in graph.items():
    done_list = []
    peer = exchange(k)
    if peer not in done_list and k not in done_list:
        if peer == k:
            graph2[k] = v
            done_list.append(k)
        else:
            graph2[k] = v + graph[peer]
            done_list.extend([k,peer])
plot_dict(graph2)
