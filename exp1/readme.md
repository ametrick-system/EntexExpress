# Experiment 1

## Summary

- Protein coding genes
- Task: binary classification (1: expressed in ≥ 20 tissues, else 0)

## Environment Setup

pip install biopython
bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed promoters.bed -s -name -fo promoters.fa