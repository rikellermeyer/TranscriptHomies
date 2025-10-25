#!/usr/bin/env python3
import sys
import pandas as pd
input_file = sys.argv[1]

dfraw = pd.read_csv(input_file, sep=r",\s*", index_col=0)

dfz = (dfraw - dfraw.min())/(dfraw.max()-dfraw.min())

dfz.to_csv("dfz_out.txt", sep="\t", float_format="%.3f")
print(f'Converted reads to Zscore and wrote to file dfz_out.txt')