#!/usr/bin/env python3
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# === Load counts ===
df = pd.read_csv("testfile_PyDeseq2.csv", index_col=0)

# === Build metadata ===
meta = pd.DataFrame({
    "sample": df.columns,
})
meta["cell_line"] = meta["sample"].str.split("_").str[0]
meta["condition"] = meta["sample"].str.split("_").str[1]
meta = meta.set_index("sample")

# === Subset to Normoxia only ===
norm_meta = meta[meta["condition"] == "Normoxia"]
df_t = df[norm_meta.index].T  # DESeq2 expects samples as rows

# === Run DESeq2: PC3 vs LNCAP under Normoxia ===
dds = DeseqDataSet(
    counts=df_t,
    metadata=norm_meta,
    design_factors="cell_line"
)
dds.deseq2()

stats = DeseqStats(dds, contrast=("cell_line", "PC3", "LNCAP"))
stats.summary()
res = stats.results_df

# === Save results ===
res.to_csv("DE_PC3_vs_LNCAP_Normoxia.csv")
print("Saved results to DE_PC3_vs_LNCAP_Normoxia.csv")

