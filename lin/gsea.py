#! /usr/bin/env python3
import pandas as pd

# Hallmark gene sets
gmt_file = "h.all.v2025.1.Hs.symbols.gmt"

# Parse the GMT file
pathway_genes = {}
with open(gmt_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        pathway_name = parts[0]
        genes = parts[2:]  # Skip the description/URL
        pathway_genes[pathway_name] = genes

# Convert to a tidy DataFrame
rows = [(p, g) for p, gene_list in pathway_genes.items() for g in gene_list]
pathway_df = pd.DataFrame(rows, columns=['Pathway', 'Gene'])


print(pathway_df.head())

pathway_df.to_csv('reactome_pathways_genes.csv', index=False)




