#!/bin/bash
#SBATCH --job-name=entex-express-dnabert2-exp5
#SBATCH --output=~/EntexExpress/exp5/finetune.log
#SBATCH --error=~/EntexExpress/exp5/finetune.err
#SBATCH --time=01:00:00
#SBATCH --partition=pi_gerstein_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --gres=gpu:1

# Activate your environment
source ~/.bashrc
conda activate dnabert2

# Move to script directory
cd ~/EntexExpress/exp5

# Run training script
./finetune.sh