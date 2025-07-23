# EntexExpress: Using AI Models to Search for and Predict Patterns in EN-TEx Gene Expression Data

**Research Question:** Can DNA sequence alone (specifically from promoter or other cis-regulatory regions) predict whether a gene is expressed in a given tissue — and possibly at what level?

- **Subquestion 1 (Experiments 1-4):** How well can DNABERT2 predict whether a protein-coding gene is broadly expressed (i.e., ubiquitous, moderate, tissue-specific) based solely on its promoter sequence?
- **Subquestion 2 (Experiment 5-6):** How well can DNABERT2 predict whether a protein-coding gene is broadly expressed in the Brain_Substantia_Nigra tissue based solely on its promoter sequence?
- **Subquestion 3 (Experiment 7):** How well can DNABERT2 *with an added regression head* predict the *level of expression* of a protein-coding gene in the Brain_Substantia_Nigra tissue based solely on its promoter sequence?

## Environment Setup (DNABERT-2)

```bash
# Step 1: create virtual environment (in home directory)
module load miniconda
conda create -n dnabert2 python=3.8

# Step 2: activate virtual environment & install base software
conda activate dnabert2
git clone https://github.com/openai/triton.git
cd triton/python
pip install cmake

# Step 3: download DNABERT2 code and install other requirements/dependencies (back in home directory)
git clone https://github.com/MAGICS-LAB/DNABERT_2.git
cd DNABERT_2
python3 -m pip install -r requirements.txt
pip install biopython
```

## Environment Setup (AlphaGenome)

```bash
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env
pip install -U alphagenome
```
