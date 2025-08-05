# Must be run while in alphagenome-env conda virtual environment

from Bio import SeqIO
import pandas as pd
import numpy as np
import random
import os
import csv

from alphagenome.data import gene_annotation
from alphagenome.data import transcript as transcript_utils

''' function to convert a tsv row to a BED entry'''
def make_bed_entry(row, window=2048):
    chrom = row["chr"]
    snv_pos = int(row["ref_start"])
    half = window // 2
    start = snv_pos - half
    end = snv_pos + half + 1  # BED is half-open [start, end)

    # Parse score, or return None to skip
    try:
        ratio = float(row['ref_allele_ratio'])
        name = (
            f"|chr={chrom}"
            f"|SNV_pos={str(snv_pos)}"
            f"|ref_allele_ratio={str(ratio)}"
        )
    except (ValueError, KeyError):
        return None

    return (chrom, start, end, name, ratio)

def generate_het_snvs_bed(tissue, assay, output_bed, window):
    tsv_file = os.path.expanduser("~/LargeFiles/hetSNVs.tsv")
    rows = []

    # open file and search
    with open(tsv_file, newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            # check if this row is a positive or negative training example, add to lists accordingly
            if row["assay"] == assay and row["tissue"] == tissue:
                rows.append(row)

    # make BED file with these rows
    written = 0
    with open(output_bed, "w", newline='') as bed:
        writer = csv.writer(bed, delimiter='\t')
        for row in rows:
            writer.writerow(make_bed_entry(row))
            written += 1

    print(f"Successfully saved {output_bed} with {written} rows")

# Includes column of binary tissue specific/not, with equal 1s and 0s unless all_genes=True
# If all_genes=True, does not downsize the negatives (0s)
def generate_tissue_specific_bed(tissue, output_bed, downbp, upbp, ratio_cutoff, tpm_floor, all_genes=False):
    # -----------------------------
    # (1) Load GTF
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
    # (2) Extract promoters
    # -----------------------------
    promoters = []
    for _, row in gtf.iterrows():
        chrom = row["Chromosome"]
        strand = row["Strand"]
        gene_id = row["gene_id"]
        transcript_id = row["transcript_id"]

        tss = row["Start"] if strand == "+" else row["End"]
        start = max(0, tss - downbp)
        end = tss + upbp

        promoters.append((chrom, start, end, gene_id, transcript_id, strand))

    promoters_df = pd.DataFrame(
        promoters,
        columns=["chr", "start", "end", "gene_id", "transcript_id", "strand"]
    )
    promoters_df["gene_id"] = promoters_df["gene_id"].str.replace(r"\.\d+", "", regex=True)
    promoters_df = promoters_df.drop_duplicates(subset="gene_id", keep="first")
    
    # -----------------------------------------------
    # (3) Load GTEx & compute other tissue median TPM
    # -----------------------------------------------
    gtex_expr_path = "~/LargeFiles/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
    gtex_df = pd.read_csv(os.path.expanduser(gtex_expr_path), sep='\t', skiprows=2)

    # Clean gene_id
    gtex_df["gene_id"] = gtex_df["Name"].str.replace(r'\.\d+', '', regex=True)

    # Convert all tissue columns to numeric (avoids accidental string types)
    tissue_cols = [c for c in gtex_df.columns if c not in ["Name", "Description", "gene_id"]]
    gtex_df[tissue_cols] = gtex_df[tissue_cols].apply(pd.to_numeric, errors="coerce")
    
    # Compute median across other tissues
    other_tissues = [c for c in tissue_cols if c != tissue]
    gtex_df["Other_Median_TPM"] = gtex_df[other_tissues].median(axis=1, skipna=True)

    # -----------------------------
    # (4) Merge with promoters
    # -----------------------------
    merged = promoters_df.merge(
        gtex_df[["gene_id", tissue, "Other_Median_TPM"]],
        on="gene_id",
        how="inner"
    )

    # -----------------------------
    # (5) Classify tissue-specific
    # -----------------------------
    merged["Specificity_Ratio"] = merged[tissue] / (merged["Other_Median_TPM"] + 1e-3)
    
    merged["Tissue_Specific"] = (
        (merged["Specificity_Ratio"] > ratio_cutoff) &
        (merged[tissue] >= tpm_floor)
    ).astype(int)

    # ------------------------------------------------------------------
    # (6) Pick negatives (most ubiquitous) if want balanced binary input
    # ------------------------------------------------------------------
    if not all_genes:
        specificity_file = "~/EntexExpress/entex_data/expressed_gene.tissue_specificity/pc.txt"

        spec_df = pd.read_csv(
            os.path.expanduser(specificity_file),
            sep="\t",
            header=None,
            names=["gene_raw", "Num_Tissues"]
        )

        # Clean gene_id (strip version, match promoters)
        spec_df["gene_id"] = spec_df["gene_raw"].str.split("|").str[0].str.replace(r"\.\d+", "", regex=True)

        # Merge with our dataset
        merged = merged.merge(spec_df[["gene_id", "Num_Tissues"]], on="gene_id", how="left")

        # Separate positives and negatives
        positives = merged[merged["Tissue_Specific"] == 1]
        negatives = merged[merged["Tissue_Specific"] == 0]

        # Pick the most ubiquitous negatives (highest Num_Tissues)
        negatives_sorted = negatives.sort_values("Num_Tissues", ascending=False)
        negatives_balanced = negatives_sorted.head(len(positives))

        # Combine back and shuffle
        merged = pd.concat([positives, negatives_balanced], ignore_index=True)
        merged = merged.sample(frac=1)
 

    # Add log(TPM)
    merged["log_tpm"] = np.log(merged[tissue] + 1)

    # -----------------------------
    # (7) Format output BED
    # -----------------------------
    merged["name"] = (
        merged["gene_id"]
        + "|TPM=" + merged[tissue].round(3).astype(str)
        + "|log(TPM)=" + merged["log_tpm"].round(3).astype(str)
        + "|OtherMedian=" + merged["Other_Median_TPM"].round(3).astype(str)
        + "|SpecificityRatio=" + merged["Specificity_Ratio"].round(3).astype(str)
        + "|TissueSpecific=" + merged["Tissue_Specific"].astype(str)
    )

    bed_df_final = merged[["chr", "start", "end", "name", tissue, "log_tpm","Other_Median_TPM", "Specificity_Ratio", "strand"]]

    bed_df_final.to_csv(output_bed, sep="\t", header=False, index=False)

    print(f"BED file saved: {output_bed}")
    print(f"Total promoters: {len(merged)}")
    print(f"Tissue-specific promoters: {merged['Tissue_Specific'].sum()}")

def config_dnabert2_input(fasta, label_key, save_prefix, task, int_regression=False, cutoffs=None, split_ratio=(0.8, 0.1, 0.1)):
    '''
    Function to convert FASTA (with headers containing '|{label_key}={val}') 
    into train/dev/test CSVs for DNABERT2 finetuning
    '''
    records = []
    for record in SeqIO.parse(fasta, "fasta"):
        header = record.description
        value = None

        # Extract value of input label parameter
        for part in header.split("::")[0].split("|"):
            if part.startswith(f"{label_key}="):
                value = float(part.split("=")[1])

        if value is None:
            raise ValueError("Parameter 'label_key' must be in the input FASTA header as |label_key=value|")

        if task == "classification": # cutoffs[i] stores the lowest value of the i-th bin
            found = False
            for i in range(len(cutoffs)-1):
                if value >= cutoffs[i] and value < cutoffs[i+1]:
                    label = i
                    found = True
                    break
            if found == False:
                label = len(cutoffs)-1
        elif task == "regression":
            label = float(value)
        else:
            raise ValueError("Parameter 'task' must be 'classification' or 'regression'")

        sequence = str(record.seq).upper()
        if "N" in sequence:
            continue  # skip ambiguous sequences

        records.append((sequence, float(label) if (int_regression == False and task == "regression") else int(label)))

    # Shuffle and split into input CSVs
    random.shuffle(records)
    total = len(records)
    n_train = int(split_ratio[0] * total)
    n_dev = int(split_ratio[1] * total)

    datasets = {
        "train.csv": records[:n_train],
        "dev.csv": records[n_train:n_train + n_dev],
        "test.csv": records[n_train + n_dev:]
    }

    os.makedirs(save_prefix, exist_ok=True)
    for name, data in datasets.items():
        df = pd.DataFrame(data, columns=["sequence", "label"])
        df.to_csv(f"{save_prefix}/{name}", index=False)

    print(f"Saved {len(records)} total examples with label key {label_key} to {save_prefix}/")
    print(pd.Series([r[1] for r in records]).value_counts())

''' For control experiments '''
def config_dnabert2_input_random_bins(fasta, save_prefix, labels, split_ratio=(0.8, 0.1, 0.1)):
    records = []
    for record in SeqIO.parse(fasta, "fasta"):
        header = record.description
        value = None

        sequence = str(record.seq).upper()
        if "N" in sequence:
            continue  # skip ambiguous sequences

        records.append((sequence, random.choice(labels))) # pick label randomly

    # Shuffle and split into input CSVs
    random.shuffle(records)
    total = len(records)
    n_train = int(split_ratio[0] * total)
    n_dev = int(split_ratio[1] * total)

    datasets = {
        "train.csv": records[:n_train],
        "dev.csv": records[n_train:n_train + n_dev],
        "test.csv": records[n_train + n_dev:]
    }

    os.makedirs(save_prefix, exist_ok=True)
    for name, data in datasets.items():
        df = pd.DataFrame(data, columns=["sequence", "label"])
        df.to_csv(f"{save_prefix}/{name}", index=False)

    print(f"Saved {len(records)} total examples with random labels {labels} to {save_prefix}/")
    print(pd.Series([r[1] for r in records]).value_counts())

