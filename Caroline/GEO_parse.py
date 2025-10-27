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

print(df_ds)

