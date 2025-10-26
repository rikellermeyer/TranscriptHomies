#!/usr/bin/env python3
import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# -------------------------------
# Usage:
# python breast_cancer_DESeq2_paired_tab_remzeros_clean.py <counts.tsv> <condition1> <condition2>
# Example:
# python breast_cancer_DESeq2_paired_tab_remzeros_clean.py GSE280284_Processed_data_files.txt C CP
# -------------------------------

if len(sys.argv) != 4:
    print("Usage: python breast_cancer_DESeq2_paired_tab_remzeros_clean.py <counts.tsv> <condition1> <condition2>")
    sys.exit(1)

file, cond1, cond2 = sys.argv[1], sys.argv[2], sys.argv[3]

# --- Load expression data ---
print(f"📂 Reading file: {file}")
df = pd.read_csv(file, sep="\t", index_col=0)

# Save annotation columns if they exist
gene_symbols = df["symbol"] if "symbol" in df.columns else None
gene_names = df["name"] if "name" in df.columns else None

# Drop annotation columns before DESeq2
drop_cols = [c for c in ["symbol", "name"] if c in df.columns]
df = df.drop(columns=drop_cols)

# Select relevant samples (ending with cond1 or cond2)
samples = [c for c in df.columns if c.endswith(cond1) or c.endswith(cond2)]
counts = df[samples]

# --- Remove genes with all-zero counts ---
zero_genes = counts[(counts == 0).all(axis=1)]
filtered_counts = counts.loc[~(counts == 0).all(axis=1)]

removed_file = f"removed_zero_count_genes_{cond1}_vs_{cond2}_clean.txt"
zero_genes.to_csv(removed_file, sep="\t")
print(f"🧹 Removed {len(zero_genes)} genes with zero counts across all samples.")
print(f"📁 Saved list of removed genes to: {removed_file}")

# --- Helper: extract patient ID robustly ---
def extract_patient_id(s, cond1, cond2):
    if s.endswith(cond2):
        base = s[: -len(cond2)]
    elif s.endswith(cond1):
        base = s[: -len(cond1)]
    else:
        base = s
    if base.endswith("0") and (base[:-1] + cond2) in df.columns:
        base = base[:-1]
    return base

# --- Build metadata with patient IDs ---
patients = [extract_patient_id(s, cond1, cond2) for s in samples]
meta = pd.DataFrame({
    "sample": samples,
    "patient": patients,
    "condition": [cond1 if s.endswith(cond1) else cond2 for s in samples]
}, index=samples)

# --- Check and print tumor-normal pairing ---
print("\n🔍 Checking detected tumor-normal pairs:")
tumor_samples = [s for s in samples if s.endswith(cond1)]
normal_samples = [s for s in samples if s.endswith(cond2)]
paired_patients, missing_pairs = [], []

for t in tumor_samples:
    pid = extract_patient_id(t, cond1, cond2)
    n_candidates = [n for n in normal_samples if extract_patient_id(n, cond1, cond2) == pid]
    if n_candidates:
        print(f"  ✅ {t}  ↔  {n_candidates[0]}")
        paired_patients.append(pid)
    else:
        print(f"  ⚠️  Missing normal pair for tumor sample: {t}")
        missing_pairs.append(pid)

for n in normal_samples:
    pid = extract_patient_id(n, cond1, cond2)
    t_candidates = [t for t in tumor_samples if extract_patient_id(t, cond1, cond2) == pid]
    if not t_candidates:
        print(f"  ⚠️  Missing tumor pair for normal sample: {n}")
        missing_pairs.append(pid)

# Keep only complete pairs
if missing_pairs:
    print(f"\n⚠️ Warning: {len(missing_pairs)} unmatched samples detected. Only complete pairs will be used.")
    paired_patients = list(set(paired_patients))
    meta = meta[meta["patient"].isin(paired_patients)]
    filtered_counts = filtered_counts[meta.index]
else:
    print("\n✅ All pairs detected successfully.")

# --- Run DESeq2 with paired design (~ patient + condition) ---
dds = DeseqDataSet(
    counts=filtered_counts.T,
    metadata=meta[["patient", "condition"]],
    design_factors=["patient", "condition"]
)
dds.deseq2()

# --- Compute differential expression statistics ---
stat_res = DeseqStats(dds, contrast=["condition", cond1, cond2])
stat_res.summary()

# --- Combine results with gene info ---
results = stat_res.results_df.copy()
results.index.name = "ensembl_id"

if gene_symbols is not None:
    results = results.merge(gene_symbols, left_index=True, right_index=True, how="left")
if gene_names is not None:
    results = results.merge(gene_names, left_index=True, right_index=True, how="left")

# --- Save results (scientific notation, 3 decimal places) ---
out_all = f"DE_{cond1}_vs_{cond2}_paired_all_clean.txt"
results.to_csv(out_all, sep="\t", float_format="%.3e")

# --- Significant genes (padj < 0.05) ---
sig_results = results[results["padj"] < 0.05].sort_values("padj")
sig_file = f"DE_{cond1}_vs_{cond2}_paired_significant_clean.txt"
sig_results.to_csv(sig_file, sep="\t", float_format="%.3e")

# --- Top 20 genes ---
top20_file = f"DE_{cond1}_vs_{cond2}_paired_top20_clean.txt"
top20 = sig_results.head(20)
top20.to_csv(top20_file, sep="\t", float_format="%.3e")

# --- Summary ---
print(f"\n✅ Saved all DESeq2 results to {out_all}")
print(f"✅ Saved significant genes (padj < 0.05) to {sig_file}")
print(f"✅ Saved top 20 genes to {top20_file}")
print(f"🧬 Total genes analyzed: {len(results)}")
print(f"✨ Significant genes found: {len(sig_results)}")


#this script takes raw counts from ~100,000 genes, removes the genes with raw counts of zero, does the DESeq2 analysis, saves files with all genes, significant genes, and removed genes. It also prints out a list of the top 10 statistically significant genes. All numbers are rounded two decimal places. 