#!/usr/bin/env python3

import pandas as pd

GEO_file = "GSE280284_Processed_data_files.txt"

bulk_file = pd.read_table(GEO_file)
print(bulk_file)

head_list = bulk_file.head(n = 50)

print(bulk_file.columns)
#just use the upper 50 gene IDs to have a smaller sample size to work with
#for this use head 

# for line in head_list
#     print (line)




#calculate Z score





#safe file output   