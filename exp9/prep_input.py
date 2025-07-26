import pandas as pd
from entex_express.utils import generate_tissue_specific_bed

# config_dnabert2_input(fasta="liver_promoters.fa", label_key="TissueSpecific", save_prefix="input", task="regression")

generate_tissue_specific_bed(tissue="Liver", output_bed="liver_promoters.bed", downbp=1024, upbp=1024, ratio_cutoff=5, tpm_floor=1)


