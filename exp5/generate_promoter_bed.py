import pandas as pd
import numpy as np
import os

from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components

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

# Filter to protein-coding genes and longest transcripts
gtf = gene_annotation.filter_protein_coding(gtf)
gtf = gene_annotation.filter_to_longest_transcript(gtf)

# Extract promoter region from TSS (±1kb)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf)
promoters = transcript_extractor.extract_promoters(flank_up=1000, flank_down=1000)

# Load the full GTEx GCT file (skip first 2 header rows)
gtex_expr_path = "~/LargeFiles/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
gtex_df = pd.read_csv(os.path.expanduser(gtex_expr_path), sep='\t', skiprows=2)

# Filter to expressed genes in Brain_Substantia_nigra
gtex_df = gtex_df[gtex_df["Brain_Substantia_nigra"] > 0].copy()

# Add log-transformed TPM
gtex_df["log_tpm"] = np.log(gtex_df["Brain_Substantia_nigra"] + 1)

# Strip version suffix from Ensembl IDs
gtex_df["gene_id"] = gtex_df["Name"].str.replace(r'\.\d+', '', regex=True)

merged = promoters.merge(
    gtex_df[['gene_id', 'Brain_Substantia_nigra', 'log_tpm']],
    left_on='gene_id_clean', right_on='gene_id'
)

# Write BED file with TPM values in name or score field
bed_df = merged[['seqname', 'start_promoter', 'end_promoter', 'gene_id', 'Brain_Substantia_nigra', 'log_tpm', 'strand']]
bed_df.columns = ['chr', 'start', 'end', 'gene_id', 'tpm', 'log_tpm', 'strand']
bed_df.to_csv("promoters_brain_substantia_nigra.bed", sep='\t', header=False, index=False)

print("BED file saved: promoters_brain_substantia_nigra.bed")