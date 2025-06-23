'''
Script to extract the TSS for each gene in pc.txt

- pc.txt: Contains all the protein coding genes and how many tissues each is expressed in
- TSS: Transcription Start Site; the position in the genome where RNA polymerase starts transcribing a gene
    - Promoters are centered around the TSS
    - For + strand genes, TSS = start of chromosome
    - For - strand genes, TSS = end of chromosome
'''

import pandas as pd

# Load pc.txt
pc = pd.read_csv("./entex_data/expressed_gene.tissue_specificity/pc.txt", sep="\t", header=None, names=["gene", "num_tissues"])
pc["gene_id"] = pc["gene"].str.extract(r"(ENSG\d+)")

# Load GTF
gtf = pd.read_csv("~/LargeFiles/gencode.v48.annotation.gtf", sep="\t", comment="#", header=None)
gtf.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

# Filter for protein-coding gene entries
genes = gtf[gtf["feature"] == "gene"]
genes = genes[genes["attribute"].str.contains("gene_type \"protein_coding\"")]

# Extract Ensembl gene ID from attribute field
genes["gene_id"] = genes["attribute"].str.extract(r'gene_id "([^"]+)"')

# Merge to get TSS
merged = pc.merge(genes, on="gene_id")

# Compute TSS and promoter window
def get_promoter(row):
    if row["strand"] == "+":
        tss = row["start"]
        return pd.Series([row["chr"], max(tss - 1000, 0), tss + 200, row["gene_id"], row["num_tissues"], row["strand"]])
    else:
        tss = row["end"]
        return pd.Series([row["chr"], max(tss - 200, 0), tss + 1000, row["gene_id"], row["num_tissues"], row["strand"]])

promoters = merged.apply(get_promoter, axis=1)
promoters.columns = ["chr", "start", "end", "gene_id", "num_tissues", "strand"]

# Save as BED file
promoters.to_csv("promoters.bed", sep="\t", header=False, index=False)