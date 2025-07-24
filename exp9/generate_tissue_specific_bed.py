import pandas as pd
import numpy as np
import os

from alphagenome.data import gene_annotation
from alphagenome.data import transcript as transcript_utils

# -----------------------------
# CONFIG
# -----------------------------
tissue = "Brain_Substantia_nigra"
output_bed = "Brain_Substantia_nigra_promoters.bed"
ratio_cutoff = 5

# -----------------------------
# Load GTF
# -----------------------------
gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
gtf_path = 'gencode.v46.annotation.gtf.feather'
if not os.path.exists(gtf_path):
    print("Downloading GTF...")
    gtf = pd.read_feather(gtf_url)
    gtf.to_feather(gtf_path)
else:
    print("Loading cached GTF...")
    gtf = pd.read_feather(gtf_path)

gtf = gene_annotation.filter_protein_coding(gtf)
gtf = gene_annotation.filter_to_longest_transcript(gtf)

# -----------------------------
# Extract promoters
# -----------------------------
promoters = []
for _, row in gtf.iterrows():
    chrom = row["Chromosome"]
    strand = row["Strand"]
    gene_id = row["gene_id"]
    transcript_id = row["transcript_id"]

    tss = row["Start"] if strand == "+" else row["End"]
    start = max(0, tss - 1024)
    end = tss + 1024

    promoters.append((chrom, start, end, gene_id, transcript_id, strand))

promoters_df = pd.DataFrame(
    promoters,
    columns=["chr", "start", "end", "gene_id", "transcript_id", "strand"]
)
promoters_df["gene_id"] = promoters_df["gene_id"].str.replace(r"\.\d+", "", regex=True)
promoters_df = promoters_df.drop_duplicates(subset="gene_id", keep="first")

# -----------------------------
# Load GTEx & compute other tissue median
# -----------------------------
gtex_expr_path = "~/LargeFiles/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
gtex_df = pd.read_csv(os.path.expanduser(gtex_expr_path), sep='\t', skiprows=2)

# Clean gene_id
gtex_df["gene_id"] = gtex_df["Name"].str.replace(r'\.\d+', '', regex=True)

# Convert all tissue columns to numeric (avoids accidental string types)
tissue_cols = [c for c in gtex_df.columns if c not in ["Name", "Description", "gene_id"]]
gtex_df[tissue_cols] = gtex_df[tissue_cols].apply(pd.to_numeric, errors="coerce")

# Median across other tissues
other_tissues = [c for c in tissue_cols if c != tissue]
gtex_df["Other_Median_TPM"] = gtex_df[other_tissues].median(axis=1, skipna=True)

# -----------------------------
# Merge with promoters
# -----------------------------
merged = promoters_df.merge(
    gtex_df[["gene_id", tissue, "Other_Median_TPM"]],
    on="gene_id",
    how="inner"
)

# -----------------------------
# Classify tissue-specific
# -----------------------------

merged["Specificity_Ratio"] = merged[tissue] / (merged["Other_Median_TPM"] + 1e-3)

merged["Tissue_Specific"] = (
    (merged["Specificity_Ratio"] > ratio_cutoff) &
    (merged[tissue] >= 1)
).astype(int)

# -----------------------------
# BED formatting
# -----------------------------
merged["name"] = (
    merged["gene_id"]
    + "|TPM=" + merged[tissue].round(3).astype(str)
    + "|OtherMed=" + merged["Other_Median_TPM"].round(3).astype(str)
    + "|TissueSpecific=" + merged["Tissue_Specific"].astype(str)
)

bed_df_final = merged[["chr", "start", "end", "name", tissue, "Other_Median_TPM", "strand"]]
bed_df_final.to_csv(output_bed, sep="\t", header=False, index=False)

print(f"BED file saved: {output_bed}")
print(f"Total promoters: {len(merged)}")
print(f"Tissue-specific promoters: {merged['Tissue_Specific'].sum()}")

