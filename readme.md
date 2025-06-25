# EntexExpress: Using DNABERT2 to Search for Patterns in EN-TEx Gene Expression Data

**Research Question:** Can DNA sequence alone (specifically from promoter or other cis-regulatory regions) predict whether a gene is expressed in a given tissue — and possibly at what level?

- **Subquestion 1 (Experiments 1-2):** How well can DNABERT2 predict whether a protein-coding gene is broadly expressed (i.e., active in many tissues) based solely on its promoter sequence?

## Environment Setup

> Note: start in home directory

```bash
srun --partition=gpu_devel --gres=gpu:1 --cpus-per-task=2 --mem=8G --time=01:00:00 --pty bash

# Step 1: create virtual environment
module load miniconda
conda create -n dnabert2 python=3.8

# Step 2: activate virtual environment & install base software
conda activate dnabert2
git clone https://github.com/openai/triton.git
cd triton/python
pip install cmake

# Step 3: download DNABERT2 code and install other requirements/dependencies 
cd ~/
git clone https://github.com/MAGICS-LAB/DNABERT_2.git
cd DNABERT_2
python3 -m pip install -r requirements.txt
pip install biopython
```