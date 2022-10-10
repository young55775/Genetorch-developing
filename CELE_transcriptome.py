# 读取转录组信息
with open('WBcel235_rna.fna', 'r') as f:
    transcript = {}
    g = f.readlines()
    for i in g:
        if i[0] == '>':
            name = i.split('(')[1].split(')')[0]
            transcript[name] = ''
        else:
            transcript[name] += i.split('\n')[0]

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


def codon_usage(transcript, record_dict, codon_dict):
    for k, v in transcript.items():
        codon = orf(v)
        if codon == 'not good':
            continue
        else:
            for i in codon:
                protein = codon_dict[i]
                if protein not in record_dict.keys():
                    record_dict[protein] = {}

                if i not in record_dict[protein].keys():
                    record_dict[protein][i] = 1
                else:
                    record_dict[protein][i] += 1
    return record_dict


def orf(seq):
    seq_list = []
    start_ind = detect_start(seq)
    end_ind = detect_stop(start_ind, seq)
    if detect_start(seq) == 'none' or detect_stop(seq) == 'none' or (end_ind - start_ind + 1) % 3 != 0:
        return 'not good'
    else:
        for i in range(start_ind, end_ind - 1, 3):
            seq_list.append(''.join([seq[i], seq[i + 1], seq[i + 2]]).replace('T', 'U'))
    return seq_list


def detect_start(seq):
    for i in range(len(seq) - 2):
        if seq[i] == 'A':
            if seq[i + 1] == 'T':
                if seq[i + 2] == 'G':
                    return i
    return 'none'


def detect_stop(start_ind, seq):
    for i in range(start_ind, len(seq) - 2):
        if seq[i] == 'T':
            if seq[i + 1] == 'A':
                if seq[i + 2] == 'A' or seq[i + 2] == 'G':
                    return i
            elif seq[i + 1] == 'G' and seq[i + 2] == 'A':
                return i
    return 'none'
