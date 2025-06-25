# Experiment 1

## Summary

- Protein coding genes
- Task: binary classification (1: expressed in ≥ 20 tissues, else 0)

For setup, simply activate the virtual environment: ``conda activate dnabert2``

## Data


bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed promoters.bed -s -name -fo promoters.fa