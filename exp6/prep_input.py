from entex_express.utils import convert_to_csv_input_labels

convert_to_csv_input_labels(
    fasta_file="promoters.fa",
    output_prefix="input",
    split_ratio=(0.8, 0.1, 0.1)
)
