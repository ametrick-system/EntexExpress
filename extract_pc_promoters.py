'''
Script to extract the TSS and promoter region for each protein-coding gene in pc.txt

- pc.txt: Contains protein-coding genes and the number of tissues each is expressed in
- TSS (Transcription Start Site): 
    - The genomic coordinate where transcription begins
    - For + strand genes: TSS = start coordinate of gene
    - For - strand genes: TSS = end coordinate of gene
- Promoter region: 
    - For + strand: TSS - 1000 to TSS + 200
    - For - strand: TSS - 200 to TSS + 1000
'''

import pandas as pd

# Load pc.txt
pc = pd.read_csv("./entex_data/expressed_gene.tissue_specificity/pc.txt", sep="\t", header=None, names=["gene", "num_tissues"])
pc["gene_id"] = pc["gene"].str.extract(r"(ENSG\d+)")

# Load GTF
gtf = pd.read_csv("~/LargeFiles/gencode.v48.annotation.gtf", sep="\t", comment="#", header=None)
gtf.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

# Filter to 'gene' lines only
genes = gtf[gtf["feature"] == "gene"].copy()

# Extract gene_id, gene_type, gene_name from attributes
genes["gene_id"] = genes["attribute"].str.extract(r'gene_id "([^"]+)"')
genes["gene_id"] = genes["gene_id"].str.replace(r"\.\d+$", "", regex=True) # remove Ensembl version numbers -> stable ID
genes["gene_type"] = genes["attribute"].str.extract(r'gene_type "([^"]+)"')
genes["gene_name"] = genes["attribute"].str.extract(r'gene_name "([^"]+)"')

# Drop rows with missing values
genes.dropna(subset=["gene_id", "gene_type", "gene_name"], inplace=True)

# Filter to protein-coding genes only that are in pc.txt
print(f"# genes in pc.txt: {len(pc)}")
print(f"# protein-coding genes in GTF: {len(genes[genes['gene_type'] == 'protein_coding'])}")
pc_genes = pd.merge(pc, genes[genes["gene_type"] == "protein_coding"], on="gene_id") # merge
print(f"# genes matched: {len(pc_genes)}")

# Compute promoter region (returns new columns only)
def compute_promoter(row):
    if row["strand"] == "+":
        tss = row["start"]
        start = max(tss - 1000, 0)
        end = tss + 200
    else:
        tss = row["end"]
        start = max(tss - 200, 0)
        end = tss + 1000
    return pd.Series({
        "chr": row["chr"],
        "start": start,
        "end": end,
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "strand": row["strand"]
    })

promoters = pc_genes.apply(compute_promoter, axis=1)

# Save BED file
promoters.to_csv("promoters.bed", sep="\t", header=False, index=False)
print(f"Wrote {len(promoters)} promoter regions to promoters.bed")