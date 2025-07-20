import pandas as pd
import numpy as np
import os

from alphagenome.data import gene_annotation

# ===============================
# 1. Load GENCODE annotation
# ===============================
gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
gtf_path = 'gencode.v46.annotation.gtf.feather'

if not os.path.exists(gtf_path):
    print("Downloading GTF...")
    gtf = pd.read_feather(gtf_url)
    gtf.to_feather(gtf_path)
else:
    print("Loading cached GTF...")
    gtf = pd.read_feather(gtf_path)

# Filter to protein-coding longest transcripts
gtf = gene_annotation.filter_protein_coding(gtf)
gtf = gene_annotation.filter_to_longest_transcript(gtf)

# ===============================
# 2. Extract promoter regions
# ===============================
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
    columns=["seqname", "start_promoter", "end_promoter", "gene_id", "transcript_id", "strand"]
)
promoters_df["gene_id_clean"] = promoters_df["gene_id"].str.replace(r"\.\d+", "", regex=True)
promoters_df = promoters_df.drop_duplicates(subset="gene_id_clean", keep="first")

# ===============================
# 3. Load GTEx median TPMs
# ===============================
gtex_expr_path = "~/LargeFiles/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
gtex_df = pd.read_csv(os.path.expanduser(gtex_expr_path), sep='\t', skiprows=2)
gtex_df["gene_id_clean"] = gtex_df["Name"].str.replace(r'\.\d+', '', regex=True)

# Keep only protein-coding genes from GTF
gtex_df = gtex_df[gtex_df["gene_id_clean"].isin(promoters_df["gene_id_clean"])]

# ===============================
# 4. Compute R_cutoff and assign labels
# ===============================
tissue_cols = [c for c in gtex_df.columns if c not in ["Name", "Description", "gene_id_clean"]]
sn_col = "Brain_Substantia_nigra"

# Define expression thresholds
TPM_CUTOFF = 10
R_CUTOFF = 5

gtex_df["median_other_tissues"] = gtex_df[[c for c in tissue_cols if c != sn_col]].median(axis=1)
gtex_df["R_fold"] = gtex_df[sn_col] / (gtex_df["median_other_tissues"] + 1)
gtex_df["label"] = ((gtex_df[sn_col] >= TPM_CUTOFF) & (gtex_df["R_fold"] >= R_CUTOFF)).astype(int)

print(f"Positive (1) genes: {gtex_df['label'].sum()} / {len(gtex_df)}")

# ===============================
# 5. Merge with promoters and save BED
# ===============================
merged = promoters_df.merge(
    gtex_df[["gene_id_clean", sn_col, "R_fold", "label"]],
    on="gene_id_clean",
    how="inner"
)

merged["name"] = (
    merged["gene_id_clean"]
    + f"|TPM_SN=" + merged[sn_col].round(3).astype(str)
    + "|R_fold=" + merged["R_fold"].round(2).astype(str)
    + "|label=" + merged["label"].astype(str)
)

# ===============================
# 5. Filter negatives (Strategy 1) and balance
# ===============================

# Separate positives and low-SN negatives
pos_df = merged[merged["label"] == 1].copy()
neg_candidates = merged[(merged["label"] == 0) & (merged["Brain_Substantia_nigra"] < 2)].copy()

print(f"Positives available: {len(pos_df)}")
print(f"Low-SN negatives available: {len(neg_candidates)}")

# Downsample negatives to match positives
neg_sample = neg_candidates.sample(len(pos_df), random_state=42)

# Combine and shuffle
balanced_df = pd.concat([pos_df, neg_sample]).sample(frac=1, random_state=42).reset_index(drop=True)

# ===============================
# 6. Final BED file
# ===============================
balanced_df["name"] = (
    balanced_df["gene_id_clean"]
    + "|TPM_SN=" + balanced_df["Brain_Substantia_nigra"].round(3).astype(str)
    + "|R_fold=" + balanced_df["R_fold"].round(2).astype(str)
    + "|label=" + balanced_df["label"].astype(str)
)

bed_df_final = balanced_df[["seqname", "start_promoter", "end_promoter",
                            "name", "Brain_Substantia_nigra", "R_fold", "strand", "label"]]
bed_df_final.columns = ["chr", "start", "end", "name", "tpm_sn", "R_fold", "strand", "label"]

output_bed = "promoters_brain_substantia_nigra_balanced_lowSN.bed"
bed_df_final.to_csv(output_bed, sep='\t', header=False, index=False)

print(f"Balanced BED file saved: {output_bed}")
print(bed_df_final['label'].value_counts())


