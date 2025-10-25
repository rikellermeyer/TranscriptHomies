#!/usr/bin/env python3

#import modules needed for generating heatmap
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import sys

## adding stars
def add_stars(fdr):
    annotations = []
    for i in range(len(fdr)): #length of rows
        row = []
        for j in range(len(fdr.columns)): #length of columns
            p_val = fdr.iloc[i, j] #gets value at row i, column j

            if p_val < 0.01:
                text = '**'
            elif p_val < 0.05:
                text = '*'
            else:
                text = ''
            
            row.append(text)
        annotations.append(row)
    return annotations

#load data set 
r2_cancer = pd.read_csv("gene_correlation_R2_linear_DS.csv", index_col=0)
r2_norm = pd.read_csv("gene_correlation_R2_linear_NM.csv", index_col=0)
fdr_cancer = pd.read_csv("gene_correlation_FDR_linear_DS.csv", index_col=0)
fdr_norm = pd.read_csv("gene_correlation_FDR_linear_NM.csv", index_col=0)

#create star only annotations
ann_c = add_stars(fdr_cancer)
ann_n = add_stars(fdr_norm)

#plot correlation heatmap
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
sns.heatmap(r2_cancer, vmin=0, vmax=1, ax=axs[0], annot=ann_c, fmt='', cmap='RdBu_r')
axs[0].set_title('cancer')
sns.heatmap(r2_norm, vmin=0, vmax=1, ax=axs[1], annot=ann_n, fmt='', cmap='RdBu_r')
axs[1].set_title('normal')

#Display heatmap
plt.tight_layout()
plt.show()