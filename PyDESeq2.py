#!/usr/bin/env python3
import os, sys, pandas as pd
os.environ["PYDESEQ2_NO_PLOTS"] = "1"

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

if len(sys.argv) < 3:
    sys.exit("Usage: python PyDESeq2.py <counts.csv> <factor_name>")

counts_file, factor = sys.argv[1], sys.argv[2]

print("🔹 Loading counts...")
df = pd.read_csv(counts_file, index_col=0)

# Transpose so samples = rows, genes = columns
df_t = df.T

print("🔹 Building metadata...")
meta = pd.DataFrame({
    "cell_line": [c.split("_")[0] for c in df_t.index],
    "condition": [c.split("_")[1] for c in df_t.index]
}, index=df_t.index)

meta["cancer_type"] = meta["cell_line"].map({
    "PC3": "prostate", "LNCAP": "prostate"
}).fillna(meta["cell_line"])

print("🔹 Filtering low-count genes...")
df_t = df_t.loc[:, df_t.sum(axis=0) >= 10]

print(f"🔹 Running DESeq2 for factor: {factor}")
dds = DeseqDataSet(counts=df_t, metadata=meta, design_factors=factor)
dds.deseq2()

# Choose a valid contrast (e.g., "condition", "Hypoxia", "Normoxia")
contrast = (factor, meta[factor].unique()[0], meta[factor].unique()[1])
print(f"🔹 Computing stats for contrast: {contrast}")

stats = DeseqStats(dds, contrast=contrast)
stats.summary()

res = stats.results_df
out_file = "deseq_results.csv"
res.to_csv(out_file)
print(f"✅ Done! Results saved to {out_file}")




