protein_dic = {}
with open(r"C:\Users\guozh\Desktop\GCF_000002985.6_WBcel235_protein.faa", "r") as f:
    m = f.readlines()
    for j in m:
        if j[0] != '>':
            for i in j:
                if i != '\n':
                    if i not in protein_dic.keys():
                        protein_dic[i] = 1
                    else:
                        protein_dic[i] += 1

# tran = {'G': 'Gly', 'P': 'Pro', 'Q': 'Gln', 'L': 'Leu', 'V': 'Val', 'E': 'Glu',
#         'M': 'Met', 'D': 'Asp', 'N': 'Asn', 'K': 'Lys', 'R': 'Arg', 'H': 'His',
#         'C': 'Cys', 'T': 'Thr', 'S': 'Ser', 'W': 'Trp', 'F': 'Phe', 'I': 'Ile',
#         'A': 'Ala', 'Y': 'Tyr'}

code = {'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': '*', 'UAG': '*',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': '*', 'UGG': 'Trp',
        'CGU': 'Arg', 'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}


# 改变string在index的位置的字母变成tar
def chan(string, index, tar):
    strm = list(string)
    strm[index] = tar
    return ''.join(strm)


# G>A,C>U
def mutation(seq, code_mat, mut_dic):
    c, g = count(seq)
    orig = code_mat[seq]
    mut_list = []
    cp = seq[:]
    if c != [] or g != []:
        if not c:
            for i in g:
                cp = chan(cp, i, 'A')
                print(cp)
                change = cp[:]
                mut_list.append((change, code_mat[change]))
                cp = chan(cp, i, 'G')
        if not g:
            for i in c:
                cp = chan(cp, i, 'U')
                change = cp[:]
                mut_list.append((change, code_mat[change]))
                cp = chan(cp, i, 'C')
        if g != [] and c != []:
            for i in c:
                cp = chan(cp, i, 'U')
                change = cp[:]
                mut_list.append((change, code_mat[change]))
                cp = chan(cp, i, 'C')
            for j in g:
                cp = chan(cp, j, 'A')
                change = cp[:]
                mut_list.append((change, code_mat[change]))
                cp = chan(cp, j, 'G')
    mut_list.insert(0, (seq, 'origin'))
    if orig not in mut_dic.keys():
        mut_dic[orig] = mut_list
    else:
        mut_dic[orig].extend(mut_list)
    return mut_dic


# 返回一个序列的C和G的index位置
def count(seq):
    g = []
    c = []
    for i in enumerate(seq):
        if i[1] == 'G':
            g.append(i[0])
        if i[1] == 'C':
            c.append(i[0])
    return c, g


# 数每个变化经历了几次
def type_count(mut_dic):
    record = {}
    for k, v in mut_dic.items():
        for i in v:
            if i[1] != 'origin':
                key = k + '->' + i[1]
                if key not in record.keys():
                    record[key] = 1
                else:
                    record[key] += 1
    return record


# 去除intron
def remove_low(string):
    line = filter(lambda ch: ch not in 'atcg', string)
    return line


# 计数codon
def orf(string):
    protein = {}
    string = string.split('\n')
    for k in string:
        l = len(k)
        n = l - (l % 3) - 2
        for i in range(0, n, 3):
            codon = ''.join([k[i], k[i + 1], k[i + 2]])
            if codon not in protein.keys():
                protein[codon] = 1
            else:
                protein[codon] += 1
    return protein

