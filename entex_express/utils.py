from Bio import SeqIO
import pandas as pd
import random
import os

def convert_to_csv_input(fasta_file, bed_file, output_prefix, task, threshold, split_ratio=(0.8, 0.1, 0.1)):
    # Load gene → num_tissues from BED
    gene_to_count = {}
    with open(bed_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                gene_id = parts[3].split(".")[0] # strip Ensembl version if present
                try:
                    num_tissues = int(parts[6])
                    gene_to_count[gene_id] = num_tissues
                except ValueError:
                    continue  # skip lines with bad counts

    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        if "::" in header:
            gene_id, _ = header.split("::")
        else:
            continue

        # Get tissue count from BED
        if gene_id not in gene_to_count:
            continue  # skip genes not in the BED file

        num_tissues = gene_to_count[gene_id]
        sequence = str(record.seq).upper()
        if "N" in sequence:
            continue

        if task == "binary":
            label = 1 if num_tissues > threshold else 0
        elif task == "regression":
            label = num_tissues
        else:
            raise ValueError("task must be 'binary' or 'regression'")

        records.append((sequence, label))

    random.shuffle(records)
    total = len(records)
    n_train = int(split_ratio[0] * total)
    n_dev = int(split_ratio[1] * total)

    datasets = {
        "train.csv": records[:n_train],
        "dev.csv": records[n_train:n_train + n_dev],
        "test.csv": records[n_train + n_dev:]
    }

    os.makedirs(output_prefix, exist_ok=True)
    for name, data in datasets.items():
        df = pd.DataFrame(data, columns=["sequence", "label"])
        df.to_csv(f"{output_prefix}/{name}", index=False)

    print(f"Saved {len(records)} total examples to {output_prefix}/")