#!/bin/bash

cd ~/DNABERT_2/finetune

# Write terminal output to log file
LOGFILE="$HOME/EntexExpress/exp02_control/finetune.log"
exec &> >(tee -a "$LOGFILE")

# CONFIG
export FASTA_PATH=~/EntexExpress/exp02_control/promoters.fa

export DATA_PATH=~/EntexExpress/exp02_control/input
export OUTPUT_PATH=~/EntexExpress/exp02_control/output
export RUN_NAME=entex_express_exp02_control

export MAX_LENGTH=512 # 0.25 * (sequence length = 2048)
export LR=3e-5

echo "Generating input CSVs with random labels..."

module load miniconda
conda activate alphagenome-env # avoid module-not-found errors from utils.py

python3 -c "
from entex_express.utils import config_dnabert2_input_random_bins
config_dnabert2_input_random_bins(
    fasta='${FASTA_PATH}',
    save_prefix='/home/asm242/EntexExpress/exp02_control/input',
    labels=[0, 1, 2]
)
"

echo "Starting fine-tuning..."

conda deactivate
conda activate dnabert2

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

