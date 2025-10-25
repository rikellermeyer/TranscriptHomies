#!/usr/bin/env python3
import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Usage example:
# python breast_cancer_DESeq2_paired.py GSE280284_Processed_data_files.txt C CP

if len(sys.argv) != 4:
    print("Usage: python breast_cancer_DESeq2_paired.py <counts.tsv> <condition1> <condition2>")
    sys.exit(1)

file, cond1, cond2 = sys.argv[1], sys.argv[2], sys.argv[3]

# Load tab-delimited file
df = pd.read_csv(file, sep="\t", index_col=0)

# Keep annotation columns if present
gene_symbols = df["symbol"] if "symbol" in df.columns else None
gene_names = df["name"] if "name" in df.columns else None

# Drop annotation columns before DESeq2
drop_cols = [c for c in ["symbol", "name"] if c in df.columns]
df = df.drop(columns=drop_cols)

# Select only tumor (cond1) and normal (cond2) columns
samples = [c for c in df.columns if c.endswith(cond1) or c.endswith(cond2)]
counts = df[samples]

# Build patient IDs (strip the trailing condition label)
patients = [s.replace(cond1, "").replace(cond2, "") for s in samples]

# Create metadata table
meta = pd.DataFrame({
    "sample": samples,
    "patient": patients,
    "condition": [cond1 if s.endswith(cond1) else cond2 for s in samples]
}, index=samples)

# --- Sanity check: verify pairing structure ---
print("\n🔍 Checking detected tumor-normal pairs:")
tumor_samples = [s for s in samples if s.endswith(cond1)]
normal_samples = [s for s in samples if s.endswith(cond2)]

paired_patients = []
missing_pairs = []

for t in tumor_samples:
    pid = t.replace(cond1, "")
    n = pid + cond2
    if n in normal_samples:
        print(f"  ✅ {t}  ↔  {n}")
        paired_patients.append(pid)
    else:
        print(f"  ⚠️  Missing normal pair for tumor sample: {t}")
        missing_pairs.append(pid)

for n in normal_samples:
    pid = n.replace(cond2, "")
    t = pid + cond1
    if t not in tumor_samples:
        print(f"  ⚠️  Missing tumor pair for normal sample: {n}")
        missing_pairs.append(pid)

if len(missing_pairs) > 0:
    print(f"\n⚠️ Warning: {len(missing_pairs)} unmatched samples detected. "
          f"Only complete pairs will be used.")
    # Keep only samples with both tumor and normal present
    paired_patients = list(set(paired_patients))
    meta = meta[meta["patient"].isin(paired_patients)]
    counts = counts[meta.index]
else:
    print("\n✅ All pairs detected successfully.")

# --- Run DESeq2 with paired design (~ patient + condition) ---
dds = DeseqDataSet(
    counts=counts.T,
    metadata=meta[["patient", "condition"]],
    design_factors=["patient", "condition"]
)
dds.deseq2()

# Compute differential expression stats
stat_res = DeseqStats(dds, contrast=["condition", cond1, cond2])
stat_res.summary()

# Combine DESeq2 results with gene info
results = stat_res.results_df.copy()
results.index.name = "ensembl_id"

if gene_symbols is not None:
    results = results.merge(gene_symbols, left_index=True, right_index=True, how="left")
if gene_names is not None:
    results = results.merge(gene_names, left_index=True, right_index=True, how="left")

# Save all results
out_file = f"DE_{cond1}_vs_{cond2}_paired.csv"
results.to_csv(out_file)

# Save significant genes (padj < 0.05)
sig_results = results[results["padj"] < 0.05].sort_values("padj")
sig_file = f"DE_{cond1}_vs_{cond2}_paired_significant.csv"
sig_results.to_csv(sig_file)

# Summary
print(f"\n✅ Saved paired DESeq2 results to {out_file}")
print(f"✅ Saved significant genes (padj < 0.05) to {sig_file}")
print(f"🧬 Total genes analyzed: {len(results)}")
print(f"✨ Significant genes found: {len(sig_results)}")
