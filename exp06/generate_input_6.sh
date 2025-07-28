#!/bin/bash

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp6/generate_input_6.log"
exec &> >(tee -a "$LOGFILE")

export RUN_TISSUE="Brain_Substantia_nigra"
export RUN_OUTPUT_BED="brain_sn_promoters.bed"
export RUN_DOWNBP=1024
export RUN_UPBP=1024
export RUN_RATIO_CUTOFF=5
export RUN_TPM_FLOOR=1
export RUN_FASTA="brain_sn_promoters.fa"
export RUN_LABEL_KEY="TissueSpecific"
export RUN_TASK="regression"
export RUN_INT_REGRESSION=True

module load miniconda
conda activate alphagenome-env

echo "Generating BED file..."
python3 -c "
from entex_express.utils import generate_tissue_specific_bed
generate_tissue_specific_bed(
    tissue='${RUN_TISSUE}',
    output_bed='${RUN_OUTPUT_BED}',
    downbp=${RUN_DOWNBP},
    upbp=${RUN_UPBP},
    ratio_cutoff=${RUN_RATIO_CUTOFF},
    tpm_floor=${RUN_TPM_FLOOR}
)
"

conda deactivate
conda activate dnabert2

echo "Generating FASTA file..."
bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed ${RUN_OUTPUT_BED} -s -name -fo ${RUN_FASTA}


conda deactivate
conda activate alphagenome-env # to avoid module not found errors (this part will run under either env)

echo "Generating DNABERT2 input CSVs..."
python3 -c "
from entex_express.utils import config_dnabert2_input
config_dnabert2_input(
    fasta='${RUN_FASTA}',
    label_key='${RUN_LABEL_KEY}',
    save_prefix='input',
    task='${RUN_TASK}',
	int_regression='${RUN_INT_REGRESSION}'
)
"

echo "Generated DNABERT2 input successfully!"
