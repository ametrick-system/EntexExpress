from Bio import SeqIO
import pandas as pd
import random
import os

def convert_to_csv_input_tpms(fasta_file, output_prefix, task, cutoffs=None, split_ratio=(0.8, 0.1, 0.1)):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        gene_id, _, _, log_tpm_str = header.split("|")
        log_tpm = float(log_tpm_str.split("=")[1].split("::")[0])
        
        sequence = str(record.seq).upper()
        if "N" in sequence:
            continue

        if task == "bins": # cutoffs[i] stores the lowest value of the i-th bin
            found = False
            for i in range(len(cutoffs)-1):
                if log_tpm >= cutoffs[i] and log_tpm < cutoffs[i+1]:
                    label = i
                    found = True
            if found == False:
                label = len(cutoffs)-1
        elif task == "regression":
            label = log_tpm
        else:
            raise ValueError("task must be 'bins' or 'regression'")

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

def convert_to_csv_input(fasta_file, bed_file, output_prefix, task, cutoffs, split_ratio=(0.8, 0.1, 0.1)):
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

        if task == "bins": # cutoffs[i] stores the lowest value of the i-th bin
            found = False
            for i in range(len(cutoffs)-1):
                if num_tissues >= cutoffs[i] and num_tissues < cutoffs[i+1]:
                    label = i
                    found = True
            if found == False:
                label = len(cutoffs)-1
        elif task == "regression":
            label = num_tissues
        else:
            raise ValueError("task must be 'bins' or 'regression'")

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

def convert_to_csv_input_random_labels(fasta_file, bed_file, output_prefix, labels, split_ratio=(0.8, 0.1, 0.1)):
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

        ########### RANDOM LABELLING #########
        label = random.choice(labels)
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