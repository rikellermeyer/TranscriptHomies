#!/usr/bin/env python3

import pandas as pd

filename = "GSE280284_Processed_data_files.txt"
gene_identifier = "gene_id"
ends_with_experimental = "CP"
ends_with_control = "C"

#create variable datafile to read the text document with Panda
# sep="\s+" is a regular expression (regex) that tells pandas how to split the columns when reading the file
#s= means one or more white spaces (dummy data set may have more than one, had issue running it)
#engine = python is needed when regular expression is needed for reading the document (here s+)
df = pd.read_csv(filename, sep="\t", engine="python")

#check headers to understand how to organise experimental versus control columns
print(df.columns.tolist())

#remove any white space in column names that may interfere with naming
df.columns = df.columns.str.strip()

#divide it into two different data tables to run for loop: control (normal) & experimental (cancer)
#using panda package you can merge the gene_ID as first column with any columns that contain sample identifier

columns_to_keep_in_table_control = [gene_identifier] + [col for col in df.columns if col.endswith(ends_with_control)]
columns_to_keep_in_table_experimental = [gene_identifier] + [col for col in df.columns if col.endswith(ends_with_experimental)]

#df slices the DataFrame to keep only specific columns, that are here listed as columns_to_keep_in_table_experimental and columns_to_keep_in_table_control
dummy_list_experimental = df[columns_to_keep_in_table_experimental]
dummy_list_control = df[columns_to_keep_in_table_control]

# Display or use the trimmed DataFrame
print(dummy_list_experimental)
print(dummy_list_control)

#make a list in a dictionary: gene_id as key, list of counts of each sample as corresponding value
#it iterates through all columns that are not "g_IDs" and adds to list
# for this prepare one table with columns for the list and one with where gene_id is determined as key

ds_columns_control = [col for col in df.columns if col.endswith(ends_with_control)]
df_ds = df[[gene_identifier] + ds_columns_control]
gene_dict_control = df_ds.set_index(gene_identifier)[ds_columns_control].apply(list, axis=1).to_dict()


ds_columns_experimental = [col for col in df.columns if col.endswith(ends_with_experimental)]
df_ds = df[[gene_identifier] + ds_columns_experimental]
gene_dict_experimental = df_ds.set_index(gene_identifier)[ds_columns_experimental].apply(list, axis=1).to_dict()

# Print the result
print(gene_dict_control)
print(gene_dict_experimental)