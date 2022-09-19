import pandas as pd


def add_candidate(readfile, gene):
    ge = readfile.candidate
    co = readfile.co_data
    data = readfile.result
    name = data[data['gene'] == gene]['sample'].to_list()[0]
    for j in name:
        df1 = co[co['tag'] == j]
        df2 = df1[df1['gene'] == gene]['protein'].to_list()[0]
        ge[j].append([gene,df2])
    readfile.candidate = ge


def del_candidate(readfile, gene):
    ge = readfile.candidate
    for k, v in ge.items():
        if gene in v:
            m = list(filter((1).__ne__, v))
            ge[k] = m
    readfile.candidate = ge

def genelist_out(readfile,path):
    ge = readfile.candidate
    ge = pd.DataFrame.from_dict(ge,orient='index')
    ge.to_csv(path)

def genelist_cls(readfile):
    readfile.candidate = {}
    for i in readfile.names:
        readfile.candidate[i] = []
