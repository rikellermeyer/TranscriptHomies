#!/usr/bin/env python3

import pandas as pd

def make_dictionary_from_raw_data_for_visualisation(filename):

    # Read the file
    df = pd.read_csv(filename, sep="\t", engine="python")

    # Clean column names
    df.columns = df.columns.str.strip()

    # Extract gene identifier column
    gene_identifier = df['symbol']

    # Extract experimental and control columns
    df_experimental = df.iloc[:, 1:7]  # columns 1 to 6
    df_control = df.iloc[:, 7:12]      # columns 7 to 11

    # Combine gene identifier with experimental/control data
    df_ds_experimental = pd.concat([gene_identifier, df_experimental], axis=1)
    df_ds_control = pd.concat([gene_identifier, df_control], axis=1)

    # Create dictionaries
    gene_dict_experimental = df_ds_experimental.set_index('symbol').apply(list, axis=1).to_dict()
    gene_dict_control = df_ds_control.set_index('symbol').apply(list, axis=1).to_dict()

    return gene_dict_control, gene_dict_experimental

# test run
filename = "GSE280284_Processed_data_files.txt"
dictionary_complete_tuple = make_dictionary_from_raw_data_for_visualisation(filename)
# call dictionaries from tuple that is returned by function
control_dict = dictionary_complete_tuple[0] 
experimental_dict = dictionary_complete_tuple[1]

print(f'control_dict: {control_dict}')
print(f'experimental: {experimental_dict}')
