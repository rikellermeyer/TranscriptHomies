# TranscriptHomies 🧬🤝

**Members: Grace Beggs, Caroline Harrer, HeaJin Hong, Tess Kelly, Zilin Xianyu**

**TAs: Riley Kellermeyer, Bekah Kim**


## Objective: 
Build a tool that identifies and visualizes gene–gene expression correlations between two biological groups (e.g., diseased vs. normal samples).

### General Project Outline: 

1. 
2.
3.
4.

## Input data search and formatting (Caroline)

(Include text here)
(if you want to insert an image, put images in folder, insert using [!filename.png])

```
insert code here




```


## Data Organization (Grace)

(include text here)

```
insert code here




```


## Correlation Analysis (Zilin and Tess)

(include text here)

```
insert code here




```

## Visualization and Output (HeaJin)
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

```
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

```
#!/usr/bin/env python

from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

#Raw data in a dictionary format. 
raw_data = {'gene_1': [30, 20, 50, 60, 100, 200], 'gene_2': [40, 40, 50, 60, 100, 500]} #this is the parsed dictionary from the raw reads
data = {'X': raw_data['gene_1'], 'Y': raw_data['gene_2']}

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

##Acknowledgement

(group photo goes here)





