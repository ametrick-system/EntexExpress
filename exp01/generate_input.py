from Bio import SeqIO
import pandas as pd

from entex_express.utils import config_dnabert2_input

# Read BED with labels
bed = pd.read_csv("promoters.bed", sep="\t", header=None, names=["chr", "start", "end", "gene_id", "gene_name", "strand", "num_tissues", "label"])
bed["key"] = bed["gene_id"]

# Read original FASTA
records = list(SeqIO.parse("promoters.fa", "fasta"))
new_records = []

for record in records:
    gene_id = record.id.split("::")[0]
    label_row = bed[bed["key"] == gene_id]
    if len(label_row) == 1:
        label = label_row.iloc[0]["label"]
        record.description = f"{gene_id}|label={label}"
        new_records.append(record)
    else:
        print(f"Could not find label for {gene_id}, skipping...")

# Save updated FASTA
SeqIO.write(new_records, "promoters_labeled.fa", "fasta")

config_dnabert2_input(
    fasta='promoters_labeled.fa',
    label_key='label',
    save_prefix='input',
    task='classification',
    cutoffs=[0, 1]
)
