# TranscriptHomies 🧬🤝

**Members: Grace Beggs, Caroline Harrer, HeaJin Hong, Tess Kelly, Zilin Xianyu**

**TAs: Riley Kellermeyer, Bekah Kim**


## Objective: 
Build a tool that identifies and visualizes gene–gene expression correlations between two biological groups (e.g., diseased vs. normal samples). Currently, there are analysis workflows for analyzing RNA-seq data to determine differential gene expression between 2 groups.
However, what if you want to look at how differential gene expression is correlated between multiple genes from a dataset?
Our tool addresses this problem.

### General Pipeline Outline: 

1. Input RNA-seq Dataset of raw counts (genes as rows, samples as columns)
2. Determine which genes are differentially expressed between two groups
3. Perform correlation statistics
4. Represent correlation matrices as heatmaps
5. Users can view the correlation between two genes as scatter plots generated from raw reads

## Data Organization (Grace)

Pipeline workflow
![alt text](TH_flowchart.png)

## Input data search and formatting (Caroline)
Public databases (including NCBI Gene Expression Omnibus (GEO)) were screened to identify bulk RNA-seq datasets comparing control and experimental conditions, with an emphasis on cancer-related research. A breast cancer dataset comprising paired samples of adjacent normal (control) and tumor (experimental) tissues from n = 6 patients was selected for subsequent analysis (GEO accession: GSE280284, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280284 , accessed 10/25/2025).

>**Generation of a dictionary of lists for visualisation of raw counts**
For the visualization of gene expression levels across samples from raw data:
- the `pandas` package was imported
the `raw_data.txt` file was formatted into two separate (“\t”-delimited) tables — listing the gene IDs and raw gene counts for samples from the control (normal adjacent tissue, Table 1) versus the experimental group (cancer tissue, Table 2)
- a dictionary of lists was created for each table, using the gene identifier (here, "gene_id") as the key and a list of expression values for each gene across samples as the corresponding value
- the function returns a tuple of both dictionaries containing `gene_dict_control` and `gene_dict_experimental`
- to access each dictionary separately, the desired dictionary can be called from the tuple

```javascript
#!/usr/bin/env python3

import pandas as pd

gene_identifier = "gene_id"
ends_with_experimental = "CP"
ends_with_control = "C"

def make_dictionary_from_raw_data_for_visualisation (filename):

    df = pd.read_csv(filename, sep="\t", engine="python")

    #  print(df.columns.tolist())
    #  df.columns = df.columns.str.strip()

    ds_columns_control = [col for col in df.columns if col.endswith(ends_with_control)]
    df_ds = df[[gene_identifier] + ds_columns_control]
    gene_dict_control = df_ds.set_index(gene_identifier)[ds_columns_control].apply(list, axis=1).to_dict()

    ds_columns_experimental = [col for col in df.columns if col.endswith(ends_with_experimental)]
    df_ds = df[[gene_identifier] + ds_columns_experimental]
    gene_dict_experimental = df_ds.set_index(gene_identifier)[ds_columns_experimental].apply(list, axis=1).to_dict()

    return gene_dict_control, gene_dict_experimental 

# filename = "GSE280284_Processed_data_files.txt"
# dictionary_complete_tuple = make_dictionary_from_raw_data_for_visualisation(filename)
# control_dict = dictionary_complete_tuple[0] 
# experimental_dict = dictionary_complete_tuple[1]

# print(f'control_dict: {control_dict}')
# print(f'experimental: {experimental_dict}')

```



## Correlation Analysis (Zilin)

Inputs: 
1.  gene expression matrix (genes x samples) from bulk RNA-Seq
2.  Pathway of interest
Output:
Gene-gene correlation matrix

```
```javascript
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

```


```

## Visualization and Output (HeaJin & Lin)

>**Net Graph**

```javascript
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

G = nx.DiGraph()

# Add edges with weights (probabilities)
n = 0
maxN = 10000
for from_gene in transition_df.index:
    for to_gene in transition_df.columns:
        fdr = transition_fdr.loc[from_gene, to_gene]
        weight = transition_df.loc[from_gene, to_gene]
        if abs(weight) > 0 and from_gene != to_gene and fdr < 0.01:
            G.add_edge(from_gene, to_gene, weight=weight)
            n += 1
            if n >= maxN:
                break
        if n >= maxN:
            break

# adjusted spacing for nodes
pos = nx.spring_layout(G, seed=42, k=0.5)  # Adjust 'k' for more spacing, larger k = more spacing
# or circular_layout for even spacing
# pos = nx.circular_layout(G)

# Get edge weights
weights = [G[u][v]['weight'] * 5 for u, v in G.edges()]  # scale up for visibility

# Set up colormap for edge colors
cmap = plt.cm.coolwarm

# Normalize weights between 0 and 1 for colormap scaling
norm = plt.Normalize(vmin=min(weights), vmax=max(weights))
edge_colors = [cmap(norm(w)) for w in weights]

# Set up the plot
fig, ax = plt.subplots(figsize=(12, 10), constrained_layout=True)

# Draw the network
nx.draw(
    G,
    pos,
    with_labels=True,
    node_size=1000,
    node_color="lightgrey",
    edge_color=edge_colors,
    width=2,
    font_size=7,  # Adjust font size for labels
    arrows=True,
    ax=ax
)

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
cbar.set_label("Correlation Strength")

# Edge labels with smaller font size
edge_labels = {(u, v): f"{G[u][v]['weight']:.2f}" for u, v in G.edges()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="black", font_size=6, ax=ax)

# Title
plt.title("Gene Interaction Network (Edge Color = Correlation Strength)", fontsize=14)

# Show the plot
plt.show()

```
<img width="798" height="668" alt="Screenshot 2025-10-27 at 12 06 08 PM" src="https://github.com/user-attachments/assets/ba27385e-2013-4a81-aa11-a54697937038" />
<img width="816" height="677" alt="Screenshot 2025-10-27 at 12 01 35 PM" src="https://github.com/user-attachments/assets/05ad2a6b-4767-4528-9f7b-75a7752d2d45" />


>**HEATMAP**
* The heat map will be a 2D graphical representation of how genes are correlated to each other (based on the R² values).
* The heat map will be generated on cancer and adjacent normal cells. 
* Since the heatmap will not distinguish + vs. - correlation, we will plot them separately
* The heatmap script was adapted from the following sources: 
[Link1](https://www.geeksforgeeks.org/python/how-to-create-a-seaborn-correlation-heatmap-in-python/)
[Link2](https://medium.com/@szabo.bibor/how-to-create-a-seaborn-correlation-heatmap-in-python-834c0686b88e)
[Link3](https://python-graph-gallery.com/92-control-color-in-seaborn-heatmaps/)

* The heat map will be generated using Seaborn module
* To install Seaborn, type the following command in the terminal: 
`pip install seaborn`

```javascript
#!/usr/bin/env python3

#import modules needed for generating heatmap
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 

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

```
>**SCATTERPLOT**
* This will generate a close-up view of how the raw read counts (expression level) of the two genes are correlated.
* The scatterplot will be fitted using the linear regression model

* To install scipy, type the following command in the terminal: 
`pip install scipy`

```javascript
#!/usr/bin/env python

#!/usr/bin/env python
#pip install scipy ## this is a unix command
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys
from dictionary_for_visualistion_final import make_dictionary_from_raw_data_for_visualisation

#Raw data will be prepared in a dictionary format. 
raw_orig = make_dictionary_from_raw_data_for_visualisation('GSE280284_Processed_data_files.txt')
ctrl_dict = raw_orig[0] #this is the parsed dictionary from the raw reads
exp_dict = raw_orig[1]
# print(ctrl_dict)
# print(exp_dict)

geo1 = exp_dict[sys.argv[0]] #geneID of interest
geo2 = exp_dict[sys.argv[1]] #geneID of interest
# print(geo1)
# print(geo2)
data = {'X': geo1, 'Y': geo2}

df = DataFrame(data, columns = ['X', 'Y'])
m, b, r_value, p_value, std_err = scipy.stats.linregress(df['X'], df['Y'])

fig, ax = plt.subplots() #this creates a blank figure and axis for plotting 
                         #(fig= entire plot window) (ax = the actual plot area,labels,styling)
#plot the data points
ax.scatter(df['X'], df['Y'], color = 'steelblue', label = 'title of the scatterplot')

#plot the regression line
ax.plot(df['X'], m*df['X'] + b, color = 'red', linewidth = 2, label = 'regression')

#add statistics 
ax.annotate(f'R² = {r_value**2:.3f}', xy=(0.05, 0.95), xycoords='axes fraction')
ax.annotate(f'y = {m:.2f}x + {b:.2f}', xy=(0.05, 0.88), xycoords='axes fraction')
ax.annotate(f'p = {p_value:.2e}', xy=(0.05, 0.81), xycoords='axes fraction')

#labels and styling
ax.set_xlabel('Gene 1 Expression (Raw Counts)')
ax.set_ylabel('Gene 2 Expression (Raw Counts)')
ax.set_title('Gene Expression Correlation')

plt.show()

print(f"R² = {r_value**2:.4f}")
print(f"P-value = {p_value:.2e}")

```
>**MASTER SCRIPT USING JUPITER NOTEBOOK**

**Acknowledgement**
![alt text](group.JPG)






[def]: /Users/pfb/transcripthomies/TranscriptHomies/TH_flowchart.png
