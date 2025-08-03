#!/bin/bash

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp14/generate_input_14.log"
exec &> >(tee -a "$LOGFILE")

export RUN_TISSUE="Liver"
export RUN_OUTPUT_PREFIX="liver_promoters"
export RUN_DOWNBP=1024
export RUN_UPBP=1024
export RUN_RATIO_CUTOFF=5
export RUN_TPM_FLOOR=1
export RUN_ALL_GENES=True
export RUN_LABEL_KEY="log(TPM)"
export RUN_TASK="regression"
export RUN_INT_REGRESSION=False
export RUN_PLOT_COL="log_tpm"

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
    tpm_floor=${RUN_TPM_FLOOR},
    all_genes=${RUN_ALL_GENES}
)
"

conda deactivate
conda activate dnabert2

echo "Generating FASTA file..."
bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed ${RUN_OUTPUT_PREFIX}.bed -s -name -fo ${RUN_OUTPUT_PREFIX}.fa

echo "Visualizing the data..."
python3 -c "
from entex_express.visualizations import plot_histogram_from_bed
plot_histogram_from_bed(
    save_path='.',
    bed_path='${RUN_OUTPUT_PREFIX}.bed',
    bed_cols=['chr', 'start', 'end', 'name', '${RUN_TISSUE}', 'log_tpm','Other_Median_TPM', 'strand'],
    plot_col='${RUN_PLOT_COL}',
    fig_name='${RUN_TISSUE}_${RUN_PLOT_COL}_histogram.png'
)
"

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

