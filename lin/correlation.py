#! /usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.stats import pearsonr, spearmanr
from scipy.stats import false_discovery_control

# input 1: gene expression matrix (genes x samples)
df = pd.read_csv("dummy_zscore2.txt", sep=r"\s+", index_col=0)
print(df.head())
# print(df.shape)        # should be (5, 6)
# print(df.index)        # should be ['gene1', 'gene2', ...]
# print(df.columns)      # should be ['s1NM', 's2NM', ...]
# print(df.dtypes)       # should all be float

# input 2: pathway
pathway = pd.read_csv("reactome_pathways_genes.csv")
print(pathway.head()) 


# Split columns by substring
nm_cols = [c for c in df.columns if "NM" in c]
ds_cols = [c for c in df.columns if "DS" in c]

# Create two separate DataFrames
df_nm = df[nm_cols]
df_ds = df[ds_cols]

print("NM subset:\n", df_nm.head(), "\n")
print("DS subset:\n", df_ds.head(), "\n")


# df_nm.to_csv("dummy_NM_only.csv")
# df_ds.to_csv("dummy_DS_only.csv")


# Z-score normalization across genes for each sample (axis=0)
#df = df.astype(float)
#expression_df = df.apply(lambda x: zscore(x), axis=1, result_type='broadcast')
#print(expression_df.dtypes)

expression_df = df_ds # using only NM samples for correlation
name = "DS"
expression_df = expression_df.loc[expression_df.var(axis=1) > 1e-6] 
print(expression_df)

# gene-gene correlation matrix
genes = expression_df.index
#print(genes)
n_genes = len(genes)

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
# explain why FDR correction is needed
pearson_fdr = pd.DataFrame(
    false_discovery_control(pearson_pval.values.flatten(), method='bh').reshape(n_genes, n_genes),
    index=genes, columns=genes
)

spearman_fdr = pd.DataFrame(
    false_discovery_control(spearman_pval.values.flatten(), method='bh').reshape(n_genes, n_genes),
    index=genes, columns=genes
)

print("Pearson Correlation Matrix:")
print(pearson_corr)
print("Pearson FDR Matrix:")
print(pearson_fdr)
print("Spearman Correlation Matrix:")
print(spearman_corr)
print("Spearman FDR Matrix:")
print(spearman_fdr) 

# Define FDR significance threshold
fdr_threshold = 0.05

# Count significant gene–gene pairs (excluding self-pairs)
pearson_sig = np.sum((pearson_fdr < fdr_threshold).values) - len(genes)
spearman_sig = np.sum((spearman_fdr < fdr_threshold).values) - len(genes)

print(f"Significant Pearson pairs (FDR<{fdr_threshold}): {pearson_sig}")
print(f"Significant Spearman pairs (FDR<{fdr_threshold}): {spearman_sig}")

# Compare counts
# highlight pearson vs spearman choosen method
if pearson_sig > spearman_sig:
    chosen_method = "pearson"
elif spearman_sig > pearson_sig:
    chosen_method = "spearman"
else:
    # Tie → compare average |r²|
    pearson_r2 = np.nanmean(np.square(pearson_corr.values))
    spearman_r2 = np.nanmean(np.square(spearman_corr.values))
    chosen_method = "pearson" if pearson_r2 > spearman_r2 else "spearman"

print(f"Chosen method: {chosen_method.upper()}")

if chosen_method == "pearson":
    pearson_corr.to_csv(f"pearson_gene_correlation_r2_{name}.csv")
    pearson_fdr.to_csv(f"pearson_gene_correlation_fdr_{name}.csv")
else:
    spearman_corr.to_csv(f"spearman_gene_correlation_r2_{name}.csv")
    spearman_fdr.to_csv(f"spearman_gene_correlation_fdr_{name}.csv")


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

#####################
# PCA plot based on expression_df
######################
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

data = df.loc[df.var(axis=1) > 1e-6] 
print("Data for PCA:")
print(data)

# expression_df: rows = genes, columns = samples
data_T = data.T  # now rows = samples, columns = genes
print("Transposed Data for PCA:")
print(data_T)

# Define a function to assign group based on column names
def assign_group(sample_name):
    if 'DS' in sample_name:
        return 'Disease'
    elif 'NM' in sample_name:
        return 'Normal'
    else:
        return 'Unknown'
# Apply to all sample names (index after transpose)
groups = [assign_group(s) for s in data_T.index]
print(groups)

# Initialize PCA (2 components for a 2D plot)
pca = PCA(n_components=2)

# Fit PCA to the data
pca_result = pca.fit_transform(data_T)

# Create a DataFrame for plotting
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'], index=data_T.index)

# annotation
pca_df['Group'] = groups
print("PCA Result DataFrame:")
print(pca_df)

plt.figure(figsize=(8,6))
sns.scatterplot(x='PC1', y='PC2', hue='Group', data=pca_df, s=100, alpha=0.7)  # added alpha
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA of Gene Expression")
plt.legend(title='Group')
plt.show()



plt.figure(figsize=(6,4))
plt.bar(range(1, len(pca.explained_variance_ratio_)+1), pca.explained_variance_ratio_*100)
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained (%)')
plt.title('Variance Explained by Each PC')
plt.show()
