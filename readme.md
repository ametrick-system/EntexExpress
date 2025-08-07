# EntexExpress: Using AI Models to Search for and Predict Patterns in EN-TEx Gene Expression Data

**Research Question:** What information can be predicted about a gene's level of expression across tissues and in general from DNA sequence alone (e.g. promoter region, start of gene, other cis-regulatory regions)?

1. **Subquestion 1 (Experiments 1-4):** Can the promoter sequence of a protein-coding gene in the human reference genome (hg38) alone predict whether this gene is broadly expressed (i.e. ubiquitous, moderate, tissue specific) in the human body?
2. **Subquestion 2 (Experiment 5):** Can the promoter sequence of a protein-coding gene in hg38 alone predict whether this gene is tissue-specific to a given tissue?
3. **Subquestion 3 (Experiment 6):** Can the promoter sequence of a protein-coding gene in the human reference genome (hg38) alone predict whether this gene is broadly expressed (i.e. ubiquitous, moderate, tissue specific) in the human body?
4. **Subquestion 4:** Can the promoter sequence of a protein-coding gene in an individual’s *personal genome* alone predict whether this gene is tissue-specific to a given tissue for this *individual*?


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
pip install biopython
```
