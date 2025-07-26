# EntexExpress: Using AI Models to Search for and Predict Patterns in EN-TEx Gene Expression Data

**Research Question:** Can DNA sequence alone (specifically from promoter or other cis-regulatory regions) predict whether a gene is expressed in a given tissue — and possibly at what level?

- **Subquestion 1 (Experiments 1-4):** How well can DNABERT2 predict whether a protein-coding gene is broadly expressed (i.e., ubiquitous, moderate, tissue-specific) based solely on its promoter sequence?
- **Subquestion 2 (Experiment 5)**: How accurately can AlphaGenome predict RNAseq for Brain_Substantia_Nigra?
- **Subquestion 3 (Experiments 6, 9):** How well can DNABERT2 classify a protein-coding gene as tissue-specific for a given tissue based solely on its promoter sequence?
- **Subquestion 4 (Experiments 7-8):** How well can DNABERT2 *with an added regression head* (DNABERT2-regression) predict the *level of expression* of *any* protein-coding gene in a given tissue based solely on its promoter sequence?
- **Subquestion 5 (Experiment 10):** How well can DNABERT2-regression predict the *level of expression* of a protein-coding gene that is *tissue-specific* to a given tissue based solely on its promoter sequence?

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
