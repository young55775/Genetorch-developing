#几个想法。1.避免连续造成的误差，将数据以1500bp为单位取均值，x，y作图后分析
        #2.进行四维统计后用异常值检测算法如IF取10%
        #3.只计算C->T G->A的突变率试一试，在鉴定得的突变中将别的突变都过滤掉
import pandas as pd
def data_cleaner(co): #3
    data = list(zip(co['chrom'].to_list(),co['pos'].to_list(),co['base'].to_list()))
    res = {}
    count = 0
    for i in data:
        if i[2][-3:] == 'C>T' or i[2][-3:] == 'G>A' or i[2][-3:] == 'C>A':
            if i[0] not in res.keys():
                res[i[0]] = []
            res[i[0]].append(i[1])
            count += 1
    print(count)
    return res