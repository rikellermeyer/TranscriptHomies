#!/usr/bin/env python3
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# -------------------------------
# Hardcoded parameters
# -------------------------------
file = "GSE280284_Processed_data_files.txt"
cond1 = "C"
cond2 = "CP"

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

# --- Filter low-count genes ---
# Filtering criteria: keep genes with >= 10 counts in >= 3 samples
MIN_COUNTS = 10
MIN_SAMPLES = 3

# Identify genes that pass the filter
genes_pass_filter = (counts >= MIN_COUNTS).sum(axis=1) >= MIN_SAMPLES
filtered_counts = counts[genes_pass_filter]
removed_genes = counts[~genes_pass_filter]

# Save removed genes
removed_file = f"removed_low_count_genes_{cond1}_vs_{cond2}.txt"
removed_genes.to_csv(removed_file, sep="\t")

print(f"\n🧹 Gene Filtering Summary:")
print(f"   Total genes before filtering: {len(counts)}")
print(f"   Genes passing filter (>={MIN_COUNTS} counts in >={MIN_SAMPLES} samples): {len(filtered_counts)}")
print(f"   Genes removed: {len(removed_genes)}")
print(f"   Percentage kept: {len(filtered_counts)/len(counts)*100:.1f}%")
print(f"📁 Saved list of removed genes to: {removed_file}")

# --- Helper: extract patient ID robustly ---
def extract_patient_id(s, cond1, cond2):
    """Strip condition suffix and normalize IDs with or without trailing zeros."""
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
patients = [extract_patient_id(s, cond1, cond2) for s in filtered_counts.columns]
meta = pd.DataFrame({
    "sample": filtered_counts.columns,
    "patient": patients,
    "condition": [cond1 if s.endswith(cond1) else cond2 for s in filtered_counts.columns]
}, index=filtered_counts.columns)


# --- Check and print tumor-normal pairing ---
print("\n🔍 Checking detected tumor-normal pairs:")
tumor_samples = [s for s in filtered_counts.columns if s.endswith(cond1)]
normal_samples = [s for s in filtered_counts.columns if s.endswith(cond2)]
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
print("\n⚙️ Running DESeq2 paired analysis...")
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

# --- Save all results (tab-delimited) ---
out_file = f"DE_{cond1}_vs_{cond2}_paired_filtered.txt"
results.to_csv(out_file, sep="\t")

# --- Save only significant results (padj < 0.05) ---
sig_results = results[results["padj"] < 0.05].sort_values("padj")
sig_file = f"DE_{cond1}_vs_{cond2}_paired_filtered_significant.txt"
sig_results.to_csv(sig_file, sep="\t")

# --- Top 20 genes (same format as significant) ---
top20_file = f"DE_{cond1}_vs_{cond2}_paired_filtered_top20.txt"
top20 = sig_results.head(20)
top20.to_csv(top20_file, sep="\t")

# Save filtered counts with gene annotations
filtered_counts_with_annotations = filtered_counts.copy()
if gene_symbols is not None:
    filtered_counts_with_annotations = filtered_counts_with_annotations.merge(
        gene_symbols, left_index=True, right_index=True, how="left"
    )
if gene_names is not None:
    filtered_counts_with_annotations = filtered_counts_with_annotations.merge(
        gene_names, left_index=True, right_index=True, how="left"
    )
filtered_counts_with_annotations.to_csv('final_input_filtered.csv', index=True)

# --- Print summary ---
print(f"\n✅ DESeq2 Analysis Complete!")
print(f"📊 Results saved to: {out_file}")
print(f"✨ Significant genes (padj < 0.05) saved to: {sig_file}")
print(f"🧬 Total genes analyzed: {len(results)}")
print(f"🔬 Significant genes found: {len(sig_results)}")
print(f"🏆 Top 20 genes saved to: {top20_file}")
print(f"📁 Filtered count matrix saved to: final_input_filtered.csv")