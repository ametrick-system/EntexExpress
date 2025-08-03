#!/bin/bash

module load miniconda
conda activate dnabert2

cd ~/DNABERT_2/finetune

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp14/finetune.log"
exec &> >(tee -a "$LOGFILE")

export DATA_PATH=~/EntexExpress/exp14/input
export OUTPUT_PATH=~/EntexExpress/exp14/output
export RUN_NAME=entex_express_exp14

export MAX_LENGTH=512 # 0.25 * (sequence length = 2048)
export LR=3e-5

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
    --num_train_epochs 5 \
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

