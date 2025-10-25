#!/usr/bin/env python3

import pandas as pd

# GEO_file = "dummy_raw.txt"

# bulk_file = pd.read_table(GEO_file)
# print(bulk_file)

# dict_control_datasets = {}
# dict_experimental_datasets = {}
# list_of_expression_values = ()

# for line in bulk_file:
# append ()

# dummy_list = 'dummy_raw.txt'

# print(dummy_list_control)


# print(dummy_list_experimental)


# df = pd.read_csv("dummy_raw.txt", sep="\t")
df = pd.read_csv("dummy_raw.txt", sep="\s+", engine="python")


print(df.columns.tolist())


df.columns = df.columns.str.strip()

#divide it into two different data tables. One control (normal), one experimental (cancer)

columns_to_keep_control = ['g_IDs'] + [col for col in df.columns if col.endswith('NM')]
columns_to_keep_experimental = ['g_IDs'] + [col for col in df.columns if col.endswith('DS')]


dummy_list_experimental = df[columns_to_keep_experimental]
dummy_list_control = df[columns_to_keep_control]

# Display or use the trimmed DataFrame


print(dummy_list_experimental)
print(dummy_list_control)

#make a dictionary with gene_id as key and then safe the rest as list
if not 'g_IDs' in df.columns:
    for col in df.columns:
        if 'g_IDs' in col:
            df.rename(columns={col: 'g_IDs'}, inplace=True)
            break

ds_columns = ['s7DS', 's8DS', 's9DS', 's10DS', 's11DS', 's12DS']
df_ds = df[['g_IDs'] + ds_columns]

gene_dict = df_ds.set_index('g_IDs')[ds_columns].apply(list, axis=1).to_dict()

# Print the result
print(gene_dict)
# build dict: key = column 1, value = list of columns 2-7

# for index in dummy_list.iterrows():
#     {row['g_IDs']: [row['s1NM'], row['s6NM']]
#     for _, row in dummy_list.iterrows()}

# dict_of_table = dummy_list.to_dict(orient='dict', *, into=<class 'dict'>, index=True)
# print(dict_of_table)


# import pandas as pd

# # Load the file and split columns manually
# df = pd.read_csv("dummy_raw.txt", sep="\s+", engine="python")

# # Strip any extra whitespace from column names
# df.columns = df.columns.str.strip()

# # Create column lists
# columns_to_keep_control = ['g_IDs'] + [col for col in df.columns if col.endswith('NM')]
# columns_to_keep_experimental = ['g_IDs'] + [col for col in df.columns if col.endswith('DS')]

# # Slice the DataFrame
# dummy_list_control = df[columns_to_keep_control]
# dummy_list_experimental = df[columns_to_keep_experimental]

# # Display the results
# print("Control (Normal) Data:")
# print(dummy_list_control)

# print("\nExperimental (Cancer) Data:")
# print(dummy_list_experimental)
