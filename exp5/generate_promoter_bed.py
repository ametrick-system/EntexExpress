import pandas as pd
import numpy as np
import os

from alphagenome.data import gene_annotation
from alphagenome.data import transcript as transcript_utils

# Download GENCODE v46 annotation (AlphaGenome-compatible)
gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
gtf_path = 'gencode.v46.annotation.gtf.feather'
if not os.path.exists(gtf_path):
    print("Downloading GTF...")
    gtf = pd.read_feather(gtf_url)
    gtf.to_feather(gtf_path)
else:
    print("Loading cached GTF...")
    gtf = pd.read_feather(gtf_path)

# GTF columns: ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'gene_id', 'gene_type', 'gene_name', 'level', 'tag', 'transcript_id', 'transcript_type', 'transcript_name', 'transcript_support_level', 'havana_transcript', 'exon_number', 'exon_id', 'hgnc_id', 'havana_gene', 'ont', 'protein_id', 'ccdsid', 'artif_dupl']

# Filter to protein-coding genes and longest transcripts
gtf = gene_annotation.filter_protein_coding(gtf)
gtf = gene_annotation.filter_to_longest_transcript(gtf)

# Extract promoter region from TSS (±1kb)
promoters = []
for _, row in gtf.iterrows():
    chrom = row["Chromosome"]
    strand = row["Strand"]
    gene_id = row["gene_id"]
    transcript_id = row["transcript_id"]

    tss = row["Start"] if strand == "+" else row["End"]
    start = max(0, tss - 1000)
    end = tss + 1000

    promoters.append((chrom, start, end, gene_id, transcript_id, strand))

# Convert to DataFrame
promoters_df = pd.DataFrame(
    promoters,
    columns=["seqname", "start_promoter", "end_promoter", "gene_id", "transcript_id", "strand"]
)
promoters_df["gene_id_clean"] = promoters_df["gene_id"].str.replace(r"\.\d+", "", regex=True)

# Keep only the first promoter per gene
promoters_df = promoters_df.drop_duplicates(subset="gene_id_clean", keep="first")

# Load the full GTEx GCT file (skip first 2 header rows)
gtex_expr_path = "~/LargeFiles/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
gtex_df = pd.read_csv(os.path.expanduser(gtex_expr_path), sep='\t', skiprows=2)

# Filter to expressed genes in Brain_Substantia_nigra
gtex_df = gtex_df[gtex_df["Brain_Substantia_nigra"] > 0].copy()
gtex_df["log_tpm"] = np.log(gtex_df["Brain_Substantia_nigra"] + 1)
gtex_df["gene_id"] = gtex_df["Name"].str.replace(r'\.\d+', '', regex=True)

# Merge
merged = promoters_df.merge(
    gtex_df[['gene_id', 'Brain_Substantia_nigra', 'log_tpm']],
    left_on='gene_id_clean',
    right_on='gene_id'
)

# Write BED file
bed_df = merged[['seqname', 'start_promoter', 'end_promoter', 'gene_id_clean', 'Brain_Substantia_nigra', 'log_tpm', 'strand']]
bed_df.columns = ['chr', 'start', 'end', 'gene_id', 'tpm', 'log_tpm', 'strand']

bed_df["name"] = (
    bed_df["gene_id"]
    + "|chr=" + bed_df["chr"]
    + ":" + bed_df["start"].astype(str) + "-" + bed_df["end"].astype(str)
    + "|TPM=" + bed_df["tpm"].round(3).astype(str)
    + "|logTPM=" + bed_df["log_tpm"].round(3).astype(str)
)

# Select and rename columns for BED
bed_df_final = bed_df[['chr', 'start', 'end', 'name', 'tpm', 'log_tpm', 'strand']]

# Write BED file
bed_df_final.to_csv("promoters_brain_substantia_nigra.bed", sep='\t', header=False, index=False)

print("BED file saved: promoters_brain_substantia_nigra.bed")
