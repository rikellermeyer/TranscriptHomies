#!/usr/bin/env python3

import pandas as pd

# GEO_file = "GSE156728_bulk_RNA.counts.txt"

print(pd.read_table("GSE156728_bulk_RNA.counts.txt"))




# with open ("/Users/pfb/group_project/TranscriptHomies/GSE156728_bulk_RNA.counts.txt", "r") as bulk_seq_text_file:

#     text_file_for_priorToZScore = pd.read_table(bulk_seq_text_file, delimiter=',', index_col=0)

#     print(text_file_for_priorToZScore)