#!/bin/bash

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp02/generate_input_2.log"
exec &> >(tee -a "$LOGFILE")

export RUN_OUTPUT_PREFIX="promoters"
export RUN_PLOT_COL="num_tissues"
export RUN_LABEL_KEY="num_tissues"
export RUN_TASK="classification"

module load miniconda
conda activate dnabert2 # to avoid module not found errors (this part will run under either env)

echo "Visualizing the data..."
python3 -c "
from entex_express.visualizations import plot_histogram_from_bed
plot_histogram_from_bed(
    save_path='.',
    bed_path='${RUN_OUTPUT_PREFIX}.bed',
    bed_cols = ['chr', 'start', 'end', 'name', 'gene_name', 'strand', 'num_tissues', 'label'],
    plot_col='${RUN_PLOT_COL}',
    fig_name='${RUN_PLOT_COL}_histogram',
    bins_list=[0, 16, 29],
    labels=['Tissue Specific -> Moderate', 'Moderate -> Ubiquitous']
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
	cutoffs=[0, 16, 29]
)
"

echo "Generated DNABERT2 input successfully!"

