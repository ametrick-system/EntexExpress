'''
Script to extract the TSS and promoter region for each protein-coding gene in pc.txt

- pc.txt: Contains protein-coding genes and the number of tissues each is expressed in
- TSS (Transcription Start Site): 
    - The genomic coordinate where transcription begins
    - For + strand genes: TSS = start coordinate of gene
    - For - strand genes: TSS = end coordinate of gene
'''

import pandas as pd
import numpy as np

# Choosing 100 bp upstream of TSS, 200 bp downstream into gene
upbp = 1000
downbp = 200

# Load pc.txt
pc = pd.read_csv("../entex_data/expressed_gene.tissue_specificity/pc.txt", sep="\t", header=None, names=["gene", "num_tissues"])
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
    tss = row["start"] if row["strand"] == "+" else row["end"]
    start = max(0, tss - downbp)
    end = tss + upbp
    return pd.Series({
        "chr": row["chr"],
        "start": start,
        "end": end,
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "strand": row["strand"],
        "num_tissues": row["num_tissues"]
    })

promoters = pc_genes.apply(compute_promoter, axis=1)

promoters["label"] = np.select(
    [
        promoters["num_tissues"] <= 15,
        (promoters["num_tissues"] > 15) & (promoters["num_tissues"] < 29),
        promoters["num_tissues"] == 29
    ],
    [0, 1, 2]
) # 3 classes, higher class => more ubiquitous

promoters["name"] = (
        promoters["gene_id"]
        + "|chr=" + promoters["chr"]
        + "|gene_id=" + promoters["gene_id"]
        + "|gene_name=" + promoters["gene_name"]
        + "|strand=" + promoters["strand"]
        + "|num_tissues=" + promoters["num_tissues"].astype(str)
        + "|label=" + promoters["label"].astype(str)
    )


# Order columns correctly for bedtools and save BED file
bed_cols = ["chr", "start", "end", "name", "gene_name", "strand", "num_tissues", "label"]
promoters[bed_cols].to_csv("promoters.bed", sep="\t", header=False, index=False)
print(f"Wrote {len(promoters)} promoter regions with expression metadata to promoters.bed")

