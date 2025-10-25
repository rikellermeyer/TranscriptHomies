#!/usr/bin/env python3

import pandas as pd

sample1_dict = {'gene1': 0.0, 'gene2': 10.0, 'gene3': 400.0, 'gene4': 8, 'gene5': 200}



def normalizelist(sample_list):
    Z_scores = []
    for item in sample_list:
        score = (item - min(sample_list))/(max(sample_list) - min(sample_list))
        Z_scores.append("{:.3f}".format(score))
    return Z_scores

def normalizedict(sample_dict):
    Z_scores = []
    gene_names = list(sample_dict.keys())
    raw_reads = list(sample_dict.values())
    for value in raw_reads:
        score = (value - min(raw_reads))/(max(raw_reads) - min(raw_reads))
        Z_scores.append("{:.3f}".format(score))
    normal_dict = dict(zip(gene_names, Z_scores))
    return normal_dict

dict1_zscore = normalizedict(sample1_dict)
df = pd.DataFrame(dict1_zscore)
with open('write_dummyz.txt', 'w') as fw:
    fw.write(df.to_string(index=False))
