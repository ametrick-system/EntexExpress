#!/bin/bash

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp13/generate_input_13.log"
exec &> >(tee -a "$LOGFILE")

export RUN_TISSUE="thoracic_aorta"
export RUN_ASSAY="TF-ChIP-seq_CTCF"
export RUN_OUTPUT_PREFIX="hetSNVs_${RUN_TISSUE}_${RUN_ASSAY}"
export RUN_WINDOW=2048
export RUN_LABEL_KEY="ref_allele_ratio"
export RUN_TASK="regression"
export RUN_INT_REGRESSION=False
export RUN_PLOT_COL="ratio"

module load miniconda
conda activate alphagenome-env

echo "Generating BED file..."
python3 -c "
from entex_express.utils import generate_het_snvs_bed
generate_het_snvs_bed(
    tissue='${RUN_TISSUE}',
    assay='${RUN_ASSAY}',
    output_bed='${RUN_OUTPUT_PREFIX}.bed',
    window=${RUN_WINDOW}
)
"

conda deactivate
conda activate dnabert2

echo "Generating FASTA file..."
bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed ${RUN_OUTPUT_PREFIX}.bed -s -name -fo ${RUN_OUTPUT_PREFIX}.fa


conda deactivate
conda activate alphagenome-env # to avoid module not found errors (this part will run under either env)

echo "Visualizing the data..."
python3 -c "
from entex_express.visualizations import plot_histogram_from_bed
plot_histogram_from_bed(
    save_path='.',
    bed_path='${RUN_OUTPUT_PREFIX}.bed',
    bed_cols=['chr', 'start', 'end', 'name', 'ratio'],
    plot_col='${RUN_PLOT_COL}',
    fig_name='${RUN_TISSUE}_${RUN_ASSAY}_histogram'
)
"

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
