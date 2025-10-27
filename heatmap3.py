#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 

def add_stars(fdr):
    """Add significance stars, excluding diagonal"""
    annotations = []
    for i in range(len(fdr)):
        row = []
        for j in range(len(fdr.columns)):
            if i == j:
                text = ''  # No annotation for diagonal
            else:
                p_val = fdr.iloc[i, j]
                if p_val < 0.01:
                    text = '**'
                elif p_val < 0.05:
                    text = '*'
                else:
                    text = ''
            row.append(text)
        annotations.append(row)
    return annotations

def find_genes_with_stars(fdr):
    """Find genes that have significant correlations"""
    genes_with_stars = set()
    
    for i in range(len(fdr)):
        gene = fdr.index[i]
        for j in range(len(fdr.columns)):
            if i != j:  # Exclude diagonal
                p_val = fdr.iloc[i, j]
                if p_val < 0.05:  # Has * or **
                    genes_with_stars.add(gene)
                    break
    return genes_with_stars

# Load datasets
r2_cancer = pd.read_csv("lin/spearman_gene_correlation_r2_DS.csv", index_col=0)
r2_norm = pd.read_csv("lin/spearman_gene_correlation_r2_NM.csv", index_col=0)
fdr_cancer = pd.read_csv("lin/spearman_gene_correlation_fdr_DS.csv", index_col=0)
fdr_norm = pd.read_csv("lin/spearman_gene_correlation_fdr_NM.csv", index_col=0)

# Find genes with significant correlations
cancer_genes_with_stars = find_genes_with_stars(fdr_cancer)
normal_genes_with_stars = find_genes_with_stars(fdr_norm)

print(f"Cancer genes with significant correlations: {len(cancer_genes_with_stars)}")
print(f"Normal genes with significant correlations: {len(normal_genes_with_stars)}")

# Create annotations
ann_c = add_stars(fdr_cancer)
ann_n = add_stars(fdr_norm)

# Cancer clustered heatmap
cancer_cluster = sns.clustermap(r2_cancer, 
                               vmin=-1, vmax=1,
                               annot=ann_c, 
                               fmt='', 
                               cmap='RdBu_r',
                               figsize=(12, 10),
                               method='ward',
                               metric='euclidean',
                               cbar_pos=(0.01, 0.3, 0.01, 0.05),  # [left, bottom, width, height]
                               cbar_kws={'label': 'R²'},
                               annot_kws={'size': 4},
                               xticklabels=True,
                               yticklabels=True)

# Adjust title positioning - move it higher and add padding
cancer_cluster.fig.suptitle('Cancer - Clustered Gene Correlations', 
                           y=0.95, fontsize=14, fontweight='bold')

# Add extra space at the top for the title
cancer_cluster.fig.subplots_adjust(top=0.92)

# Get the reordered gene names after clustering
cancer_reordered_genes = cancer_cluster.data2d.index.tolist()
cancer_reordered_cols = cancer_cluster.data2d.columns.tolist()

# Set ticks for all genes first, then modify labels
cancer_cluster.ax_heatmap.set_yticks(range(len(cancer_reordered_genes)))
cancer_cluster.ax_heatmap.set_xticks(range(len(cancer_reordered_cols)))

# Create labels showing only significant genes
cancer_ylabels = [gene if gene in cancer_genes_with_stars else '' for gene in cancer_reordered_genes]
cancer_xlabels = [gene if gene in cancer_genes_with_stars else '' for gene in cancer_reordered_cols]

cancer_cluster.ax_heatmap.set_yticklabels(cancer_ylabels, fontsize=4)
cancer_cluster.ax_heatmap.set_xticklabels(cancer_xlabels, rotation=90, fontsize=4)

# Normal clustered heatmap  
normal_cluster = sns.clustermap(r2_norm, 
                               vmin=-1, vmax=1,
                               annot=ann_n, 
                               fmt='', 
                               cmap='RdBu_r',
                               figsize=(12, 10),
                               method='ward',
                               metric='euclidean',
                               cbar_pos=(0.01, 0.3, 0.01, 0.05),  # [left, bottom, width, height]
                               cbar_kws={'label': 'R²'},
                               annot_kws={'size': 4},
                               xticklabels=True,
                               yticklabels=True)

# Adjust title positioning - move it higher and add padding
normal_cluster.fig.suptitle('Normal - Clustered Gene Correlations', 
                           y=0.95, fontsize=14, fontweight='bold')

# Add extra space at the top for the title
normal_cluster.fig.subplots_adjust(top=0.92)

# Get the reordered gene names after clustering
normal_reordered_genes = normal_cluster.data2d.index.tolist()
normal_reordered_cols = normal_cluster.data2d.columns.tolist()

# Set ticks for all genes first, then modify labels
normal_cluster.ax_heatmap.set_yticks(range(len(normal_reordered_genes)))
normal_cluster.ax_heatmap.set_xticks(range(len(normal_reordered_cols)))

# Create labels showing only significant genes
normal_ylabels = [gene if gene in normal_genes_with_stars else '' for gene in normal_reordered_genes]
normal_xlabels = [gene if gene in normal_genes_with_stars else '' for gene in normal_reordered_cols]

normal_cluster.ax_heatmap.set_yticklabels(normal_ylabels, fontsize=4)
normal_cluster.ax_heatmap.set_xticklabels(normal_xlabels, rotation=90, fontsize=4)

plt.show()

# Optional: Print which genes are being shown
print(f"\nGenes displayed in cancer heatmap: {[g for g in cancer_reordered_genes if g in cancer_genes_with_stars]}")
print(f"Genes displayed in normal heatmap: {[g for g in normal_reordered_genes if g in normal_genes_with_stars]}")