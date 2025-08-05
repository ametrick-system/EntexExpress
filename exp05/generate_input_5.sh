#!/bin/bash

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp05/generate_input_5.log"
exec &> >(tee -a "$LOGFILE")

export RUN_TISSUE="Liver"
export RUN_OUTPUT_PREFIX="liver_promoters"
export RUN_DOWNBP=1024
export RUN_UPBP=1024
export RUN_RATIO_CUTOFF=5
export RUN_TPM_FLOOR=1
export RUN_LABEL_KEY="TissueSpecific"
export RUN_TASK="regression" # not really regression, just using 0/1 labels already in fasta
export RUN_INT_REGRESSION=True

module load miniconda
conda activate alphagenome-env

echo "Generating BED file..."
python3 -c "
from entex_express.utils import generate_tissue_specific_bed
generate_tissue_specific_bed(
    tissue='${RUN_TISSUE}',
    output_bed='${RUN_OUTPUT_PREFIX}.bed',
    downbp=${RUN_DOWNBP},
    upbp=${RUN_UPBP},
    ratio_cutoff=${RUN_RATIO_CUTOFF},
    tpm_floor=${RUN_TPM_FLOOR}
)
"

conda deactivate
conda activate dnabert2

echo "Generating FASTA file..."
bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed ${RUN_OUTPUT_PREFIX}.bed -s -name -fo ${RUN_OUTPUT_PREFIX}.fa


conda deactivate
conda activate alphagenome-env # to avoid module not found errors (this part will run under either env)

echo "Generating DNABERT2 input CSVs..."
python3 -c "
from entex_express.utils import config_dnabert2_input
config_dnabert2_input(
    fasta='${RUN_OUTPUT_PREFIX}.fa',
    label_key='${RUN_LABEL_KEY}',
    save_prefix='input',
    task='${RUN_TASK}',
	int_regression=${RUN_INT_REGRESSION}
)
"

echo "Generated DNABERT2 input successfully!"
