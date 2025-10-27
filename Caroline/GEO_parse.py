# #!/usr/bin/env python3


# # pip install biopython geoparse

# #ftp
import GEOparse

gse = GEOparse.get_GEO("GSE280284", destdir="./GEOdata", how="full")
print(gse.metadata)

# #GEOparse

# import GEOparse
# import pandas as pd

# # Load the GEO Series
# gse = GEOparse.get_GEO("GSE280284", destdir="./Caroline/GEOdata")

# expression_data = {}

# for gsm_id, gsm in gse.gsms.items():
#     df = gsm.table
#     columns = [col.upper() for col in df.columns]

#     id_col = next((col for col in df.columns if col.upper() in ["ID_REF", "ID", "PROBENAME"]), None)
#     val_col = next((col for col in df.columns if col.upper() in ["VALUE", "SIGNAL", "INTENSITY"]), None)


#     if id_col and val_col:
#         expression_data[gsm_id] = df.set_index(id_col)[val_col]
#     else:
#         print(f"⚠️ Skipping {gsm_id}: couldn't find expected columns")

# # Combine into one DataFrame
# combined_df = pd.DataFrame(expression_data)

# # Save to tab-separated TXT
# combined_df.to_csv("combined_expression_matrix.txt", sep="\t")

#rename columns and let people decide on rearrangement
#or possibility to check metadata to select based on control versus experimental
#add control or experimental at the end of the headers name and sort by this
#indexing columns and choose/merge columns of interest


#rearrange column names by index

import pandas as pd

filename = "GSE280284_Processed_data_files.txt"

df = pd.read_csv(filename, sep="\t")
# print(df)
df_index = df.iloc[:, [0]]
df_experimental = df.iloc[:, [1, 2, 3, 4, 5, 6]]
df_control = df.iloc[:, list(range(7, 12))]

df_ds = pd.concat([df_index, df_experimental], axis=1)




def make_dictionary_from_raw_data_for_visualisation (filename):

#create variable datafile to read the text document with Panda
# sep="\s+" is a regular expression (regex) that tells pandas how to split the columns when reading the file
#s= means one or more white spaces (dummy data set may have more than one, had issue running it)
#engine = python is needed when regular expression is needed for reading the document (here s+)
    df = pd.read_csv(filename, sep="\t", engine="python")

#check headers to understand how to organise experimental versus control columns
    #  print(df.columns.tolist())

#remove any white space in column names that may interfere with naming
    df.columns = df.columns.str.strip()

#make a list in a dictionary: gene_id as key, list of counts of each sample as corresponding value
#it iterates through all columns that are not "g_IDs" and adds to list
# for this prepare one table with columns for the list and one with where gene_id is determined as key
    gene_identifier = df.iloc[:, [12]]
    df_experimental = df.iloc[:, [1, 2, 3, 4, 5, 6]]
    df_control = df.iloc[:, list(range(7, 12))]

    df_ds_experimental = pd.concat([gene_identifier, df_experimental], axis=1)
    df_ds_control = pd.concat([gene_identifier, df_control], axis=1)

    gene_dict_control = df_ds.set_index(gene_identifier)[ df_ds_control].apply(list, axis=1).to_dict()
   
    gene_dict_experimental = df_ds.set_index(gene_identifier)[df_ds_experimental ].apply(list, axis=1).to_dict()

    return gene_dict_control, gene_dict_experimental #provides result as tuple, to have each dictionary we have to extract each dictionary from tuple (see below)


# Print the result as control
# print(gene_dict_control)
# print(gene_dict_experimental)


filename = "GSE280284_Processed_data_files.txt"
dictionary_complete_tuple = make_dictionary_from_raw_data_for_visualisation(filename)
control_dict = dictionary_complete_tuple[0] 
experimental_dict = dictionary_complete_tuple[1]

print(f'control_dict: {control_dict}')
print(f'experimental: {experimental_dict}')

# print(df_ds)

# import yaml

# # Load YAML config
# with open("config.yaml", "r") as file:
#     config = yaml.safe_load(file)

# filename = config["filename"]

# # Read file using config
# df = pd.read_csv(filename, sep="delimiter")

# # Slice columns based on config
# df_index = df.iloc[:, config["columns"]["index"]]
# df_experimental = df.iloc[:, config["columns"]["experimental"]]
# df_control = df.iloc[:, config["columns"]["control"]]

# # Combine index and experimental columns
# df_ds = pd.concat([df_index, df_experimental], axis=1)
# print(df_ds)