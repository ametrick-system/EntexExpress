#!/bin/bash

cd ~/DNABERT_2/finetune

# Write terminal output to log file
LOGFILE="~/EntexExpress/exp2/finetune.log"
exec &> >(tee -a "$LOGFILE")

# Change for each experiment
export DATA_PATH=~/EntexExpress/exp2/input
export OUTPUT_PATH=~/EntexExpress/exp2/output
export RUN_NAME=DNABERT2_exp2_promoters

export MAX_LENGTH=300 # set to 0.25 * your sequence length (1200)
export LR=3e-5

echo "Starting fine-tuning..."

# Training use DataParallel
python3 train.py \
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
    --find_unused_parameters False
    
# Training use DistributedDataParallel (more efficient)
# export num_gpu=4 # please change the value based on your setup

# torchrun --nproc_per_node=${num_gpu} train.py \
#     --model_name_or_path zhihan1996/DNABERT-2-117M \
#     --data_path  ${DATA_PATH} \
#     --kmer -1 \
#     --run_name DNABERT2_${DATA_PATH} \
#     --model_max_length ${MAX_LENGTH} \
#     --per_device_train_batch_size 8 \
#     --per_device_eval_batch_size 16 \
#     --gradient_accumulation_steps 1 \
#     --learning_rate ${LR} \
#     --num_train_epochs 5 \
#     --fp16 \
#     --save_steps 200 \
#     --output_dir output/dnabert2 \
#     --evaluation_strategy steps \
#     --eval_steps 200 \
#     --warmup_steps 50 \
#     --logging_steps 100 \
#     --overwrite_output_dir True \
#     --log_level info \
#     --find_unused_parameters False