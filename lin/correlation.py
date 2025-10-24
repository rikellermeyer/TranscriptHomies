#! /usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.stats import pearsonr, spearmanr

# input
# expression_df = pd.read_csv("gene_expression.csv", index_col=0)

# testing
data = {
    "Sample1": [120, 90, 300, 10],
    "Sample2": [150, 100, 250, 15],
    "Sample3": [130, 95, 280, 12],
    "Sample4": [160, 110, 260, 14],
    "Sample5": [125, 98, 275, 13],
    "Sample6": [140, 105, 290, 16],
    "Sample7": [135, 102, 265, 12],
    "Sample8": [145, 108, 280, 15],
    "Sample9": [150, 95, 270, 14],
    "Sample10": [128, 100, 285, 11],
}
genes = ["GeneA", "GeneB", "GeneC", "GeneD"]
df = pd.DataFrame(data, index=genes)
df = df.astype(float)

# Z-score normalization across genes for each sample (axis=0)
expression_df = df.apply(lambda x: zscore(x), axis=1, result_type='broadcast')
#print(expression_df.dtypes)
#expression_df = df
#print(expression_df)

# gene-gene correlation matrix
genes = expression_df.index
print(genes)
n_genes = len(genes)
print(n_genes)

pearson_corr = pd.DataFrame(np.zeros((n_genes, n_genes)), index=genes, columns=genes)
pearson_pval = pd.DataFrame(np.zeros((n_genes, n_genes)), index=genes, columns=genes)

spearman_corr = pd.DataFrame(np.zeros((n_genes, n_genes)), index=genes, columns=genes)
spearman_pval = pd.DataFrame(np.zeros((n_genes, n_genes)), index=genes, columns=genes)

# pairwise correlations and p-values
for i, g1 in enumerate(genes):
    for j, g2 in enumerate(genes):
        if j >= i:
            r, p = pearsonr(expression_df.loc[g1], expression_df.loc[g2])
            pearson_corr.loc[g1, g2] = pearson_corr.loc[g2, g1] = r
            pearson_pval.loc[g1, g2] = pearson_pval.loc[g2, g1] = p

            r_s, p_s = spearmanr(expression_df.loc[g1], expression_df.loc[g2])
            spearman_corr.loc[g1, g2] = spearman_corr.loc[g2, g1] = r_s
            spearman_pval.loc[g1, g2] = spearman_pval.loc[g2, g1] = p_s


# FDR (Benjamini-Hochberg)
def fdr_bh(pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    fdr = sorted_pvals * n / (np.arange(1, n+1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]  # enforce monotonicity
    fdr[fdr > 1] = 1
    unsorted_fdr = np.empty_like(fdr)
    unsorted_fdr[sorted_idx] = fdr
    return unsorted_fdr

pearson_fdr = pd.DataFrame(
    fdr_bh(pearson_pval.values.flatten()).reshape(n_genes, n_genes),
    index=genes, columns=genes
)

spearman_fdr = pd.DataFrame(
    fdr_bh(spearman_pval.values.flatten()).reshape(n_genes, n_genes),
    index=genes, columns=genes
)


print(pearson_corr)
print(pearson_fdr)
print(spearman_corr)
print(spearman_fdr) 
pearson_corr.to_csv("gene_correlation_matrix_linear.csv")
pearson_fdr.to_csv("gene_correlation_fdr_linear.csv")

spearman_corr.to_csv("gene_correlation_matrix_non_linear.csv")
spearman_fdr.to_csv("gene_correlation_fdr_non_linear.csv")


# heatmap
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
sns.heatmap(pearson_corr, annot=True, cmap="coolwarm", center=0, fmt=".2f")
plt.title("Gene–Gene Correlation Heatmap (Linear)")
plt.show()

plt.figure(figsize=(8, 6))
sns.heatmap(spearman_corr, annot=True, cmap="coolwarm", center=0, fmt=".2f")
plt.title("Gene–Gene Correlation Heatmap (Non-Linear)")
plt.show()



