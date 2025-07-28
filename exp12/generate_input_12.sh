#!/bin/bash

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp12/generate_input_12.log"
exec &> >(tee -a "$LOGFILE")

export RUN_TISSUE="Liver"
export RUN_OUTPUT_PREFIX="liver_promoters"
export BED_TO_USE="${RUN_OUTPUT_PREFIX}.bed"
export RUN_DOWNBP=1024
export RUN_UPBP=1024
export RUN_RATIO_CUTOFF=5
export RUN_TPM_FLOOR=1
export RUN_FASTA="${RUN_OUTPUT_PREFIX}.fa"
export RUN_LABEL_KEY="SpecificityRatio"
export RUN_TASK="classification"
export RUN_INT_REGRESSION=False
export RUN_PLOT_COL="Specificity_Ratio"
export RUN_CUTOFFS="0 0.5 5 15"

module load miniconda
conda activate alphagenome-env

echo "Generating BED file..."
python3 -c "
from entex_express.utils import generate_tissue_specific_bed
generate_tissue_specific_bed(
    tissue='${RUN_TISSUE}',
    output_bed='${BED_TO_USE}',
    downbp=${RUN_DOWNBP},
    upbp=${RUN_UPBP},
    ratio_cutoff=${RUN_RATIO_CUTOFF},
    tpm_floor=${RUN_TPM_FLOOR}
)
"

conda deactivate
conda activate dnabert2

echo "Generating FASTA file..."
bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed ${BED_TO_USE} -s -name -fo ${RUN_FASTA}

conda deactivate
conda activate alphagenome-env # to avoid module not found errors (this part will run under either env)

echo "Visualizing the data..."
python3 -c "
from entex_express.visualizations import plot_histogram_from_bed
plot_histogram_from_bed(
    save_path='.',
    bed_path='${BED_TO_USE}',
    bed_cols=['chr', 'start', 'end', 'name', '${RUN_TISSUE}', 'log_tpm','Other_Median_TPM', 'Specificity_Ratio', 'strand'],
    plot_col='${RUN_PLOT_COL}',
    fig_name='${RUN_TISSUE}_${RUN_PLOT_COL}_histogram'
)
"

echo "Generating DNABERT2 input CSVs..."
python3 -c "
import os
from entex_express.utils import config_dnabert2_input
cutoffs = [float(x) for x in os.environ['RUN_CUTOFFS'].split()]
config_dnabert2_input(
    fasta='${RUN_FASTA}',
    label_key='${RUN_LABEL_KEY}',
    save_prefix='input',
    task='${RUN_TASK}',
    int_regression=${RUN_INT_REGRESSION},
    cutoffs=cutoffs
)
"

echo "Generated DNABERT2 input successfully!"
