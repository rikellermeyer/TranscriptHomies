#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# -----------------------------
# Command-line argument check
# -----------------------------
if len(sys.argv) != 4:
    print("Usage: python breast_cancer_DESeq2_tabdelim_remzeros_v2.py <counts.txt> <tumor_label> <normal_label>")
    sys.exit(1)

file, tumor_label, normal_label = sys.argv[1], sys.argv[2], sys.argv[3]

# -----------------------------
# Load counts file
# -----------------------------
print(f"📂 Reading file: {file}")
df = pd.read_csv(file, sep="\t", index_col=0)

# Identify sample columns (containing tumor or normal labels)
samples = [c for c in df.columns if tumor_label in c or normal_label in c]
counts = df[samples]

# -----------------------------
# Remove genes with all-zero counts
# -----------------------------
zero_genes = counts[(counts == 0).all(axis=1)]
filtered_counts = counts.loc[~(counts == 0).all(axis=1)]

print(f"🧹 Removed {len(zero_genes)} genes with zero counts across all samples.")
removed_file = "removed_zero_count_genes_v2.txt"
zero_genes.to_csv(removed_file, sep="\t", index=True)
print(f"📁 Saved list of removed genes to: {removed_file}")

# -----------------------------
# Detect and pair samples
# -----------------------------
patients = []
conditions = []

for sample in filtered_counts.columns:
    if sample.endswith(tumor_label):
        patient_id = sample[: -len(tumor_label)]
        condition = "tumor"
    elif sample.endswith(normal_label):
        patient_id = sample[: -len(normal_label)]
        condition = "normal"
    else:
        print(f"⚠️ Skipping unrecognized sample name: {sample}")
        continue

    patients.append(patient_id)
    conditions.append(condition)

meta = pd.DataFrame({
    "sample": filtered_counts.columns,
    "patient": patients,
    "condition": conditions
}, index=filtered_counts.columns)

# Ensure pairs exist for both conditions
paired_patients = (
    set(meta.loc[meta["condition"] == "tumor", "patient"]) &
    set(meta.loc[meta["condition"] == "normal", "patient"])
)

if len(paired_patients) == 0:
    print("❌ No properly matched pairs found. Check sample names.")
    print(f"🔍 Tip: ensure tumor samples end with '{tumor_label}' and normals end with '{normal_label}'")
    sys.exit(1)

meta = meta[meta["patient"].isin(paired_patients)]
filtered_counts = filtered_counts[meta.index]

print(f"✅ {len(paired_patients)} matched patient pairs found.")

# -----------------------------
# Run DESeq2 (paired design)
# -----------------------------
print("⚙️ Running DESeq2 paired analysis...")
dds = DeseqDataSet(counts=filtered_counts.T, metadata=meta, design_factors=["patient", "condition"])
dds.deseq2()

# -----------------------------
# Compute DESeq2 stats
# -----------------------------
print("📊 Calculating statistics...")
stat_res = DeseqStats(dds, contrast=["condition", "tumor", "normal"])
stat_res.summary()

results = stat_res.results_df.copy()
results.sort_values("padj", inplace=True)
results.dropna(subset=["padj"], inplace=True)

# If gene names/symbols exist in original file, append them
if "symbol" in df.columns:
    results["symbol"] = df.loc[results.index, "symbol"]
if "name" in df.columns:
    results["name"] = df.loc[results.index, "name"]

# -----------------------------
# Save outputs (tab-delimited)
# -----------------------------
out_all = f"DE_{tumor_label}_vs_{normal_label}_paired_all_v2.txt"
out_sig = f"DE_{tumor_label}_vs_{normal_label}_paired_significant_v2.txt"

results.to_csv(out_all, sep="\t", index=True, float_format="%.5f")
results[results["padj"] < 0.05].to_csv(out_sig, sep="\t", index=True, float_format="%.5f")

print(f"\n🧬 Total genes analyzed: {len(results)}")
print(f"✨ Significant genes (padj < 0.05): {(results['padj'] < 0.05).sum()}")
print(f"📁 Saved all results to: {out_all}")
print(f"📁 Saved significant results to: {out_sig}")
print(f"📁 Saved removed genes to: {removed_file}")

# -----------------------------
# Show top 10 most significant genes
# -----------------------------
print("\n🔝 Top 10 significant genes:")
top10 = (
    results.head(20)[["log2FoldChange", "padj", "symbol"]]
    if "symbol" in results.columns
    else results.head(20)[["log2FoldChange", "padj"]]
)
print(top20)
