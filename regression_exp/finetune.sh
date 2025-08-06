#!/bin/bash

module load miniconda
conda activate dnabert2

cd ~/DNABERT_2/finetune

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/regression_exp/finetune.log"
exec &> >(tee -a "$LOGFILE")

export DATA_PATH=~/EntexExpress/regression_exp/input
export OUTPUT_PATH=~/EntexExpress/regression_exp/output
export RUN_NAME=entex_express_regression_exp

export RUN_OUTPUT_PREFIX="/home/asm242/EntexExpress/regression_exp/peaks"
export RUN_SAVE_PATH="/home/asm242/EntexExpress/regression_exp"
export RUN_PLOT_COL="log(signal)"
export RUN_LABEL_KEY="log(signal)"
export RUN_TASK="regression"


export MAX_LENGTH=114 # 0.25 * (max sequence length = 456)
export LR=2e-5

echo "Visualizing the data..."
python3 -c "
from entex_express.visualizations import plot_histogram_from_bed
plot_histogram_from_bed(
    save_path='${RUN_SAVE_PATH}',
    bed_path='${RUN_OUTPUT_PREFIX}.bed',
    bed_cols=['chr', 'start', 'end', 'name', 'signal', 'log(signal)'],
    plot_col='${RUN_PLOT_COL}',
    fig_name='${RUN_PLOT_COL}_histogram'
)
"

conda deactivate
conda activate alphagenome-env # to avoid module-not-found errors, even though no need for those packages here

echo "Generating DNABERT2 input CSVs..."
python3 -c "
from entex_express.utils import config_dnabert2_input
config_dnabert2_input(
    fasta='${RUN_OUTPUT_PREFIX}.fa',
    label_key='${RUN_LABEL_KEY}',
    save_prefix='${RUN_SAVE_PATH}/input',
    task='${RUN_TASK}'
)
"

echo "Generated DNABERT2 input successfully!"

conda deactivate
conda activate dnabert2

echo "Starting fine-tuning..."

# Training use DataParallel
python3 train_regression.py \
    --model_name_or_path zhihan1996/DNABERT-2-117M \
    --data_path  ${DATA_PATH} \
    --kmer -1 \
    --run_name ${RUN_NAME} \
    --model_max_length ${MAX_LENGTH} \
    --per_device_train_batch_size 8 \
    --per_device_eval_batch_size 16 \
    --gradient_accumulation_steps 1 \
    --learning_rate ${LR} \
    --num_train_epochs 20 \
    --fp16 \
    --save_steps 200 \
    --output_dir ${OUTPUT_PATH} \
    --evaluation_strategy steps \
    --eval_steps 200 \
    --warmup_steps 50 \
    --logging_steps 100 \
    --overwrite_output_dir True \
    --log_level info \
    --find_unused_parameters False \
	--task regression

