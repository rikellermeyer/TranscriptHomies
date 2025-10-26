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
