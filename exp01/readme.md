# Experiment 1

## Summary
**Research Question:** Can the promoter sequence of a gene in the human reference genome (hg38) alone predict whether this gene is broadly expressed (i.e. ubiquitous, moderate, tissue specific) in the human body?

- **Input:** Promoter sequences of protein-coding genes, 2048 bp around the TSS
- **Label:** Binary (1 if gene expressed in ≥ 20 tissues, else 0)
- **Goal:** Test if general expression level correlates with promoter-intrinsic sequence features

## Raw Data
- Tissue data from Tissue_Specificity.zip ("the tissue specificity of gene expression and functional signals of cis-regulatory elements") downloaded from the EN-TEx portal. This experiment uses ``expressed_gene.tissue_specificity/pc.txt``, which contains a list of protein coding genes and how many tissues they are expressed in
- Human genome sequence from GENCODE Human Genome Release 48 primary assembly (GRCh38)

## Input
- ``extract_pc_promoters.py`` extracts the promoter region for each gene in pc.txt (1024 nucleotides upstream and 1024 nucleotides downstream) and stores results in ``promoters.bed`` with binary labels indicating whether that gene is expressed in ≥ 20 tissues with label counts 9945 0's and 6055 1's.
- ``promoters.fa`` contains the actual promoter sequences from GRCh38 based on windows in ``promoters.bed`` (generated with bedtools command ``bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed promoters.bed -s -name -fo promoters.fa``)
- ``generate_input.py`` creates the 3 input CSV files: ``train.csv`` (80%), ``test.csv`` (10%), and ``dev.csv`` (10%) to be passed to DNABERT2 in finetuning, where each line contains a promoter sequence and a 1 or 0 indicating whether that gene is expressed in ≥ 20 tissues
-  ``train.csv``, ``test.csv``, and ``dev.csv`` are all stored in the ``input`` directory

## Fine-tuning DNABERT2
The script ``finetune.sh`` runs the DNABERT2 fine tuning with the default weights, steps, epochs, and learning rate. To run, first make the file executable (``chmod +x finetune.sh``) then run (``./finetune.sh``). All the terminal output from the fine-tuning process is saved in ``finetune.log``.

## Results
Here is a table that presents the results stored in ``output/results/entex_express_exp01/eval_results.json`` (note that each value is rounded to the nearest hundredth). Additionally, ``plot_roc.py`` plots the Receiver Operating Characteristic (ROC) curve in ``roc_curve.png`` and computes the AUROC.

| Metric                      | Value  | Notes                                                                 |
|-----------------------------|--------|-----------------------------------------------------------------------|
| `Loss`                 | 0.600  | Slightly better than random (~0.693)                                  |
| `Accuracy`             | 0.696  | Significantly better than random (~0.5)                                        |
| `F1`                   | 0.660  | Good measure of accuracy with input inbalances removed     |
| `Precision`            | 0.677  | Equal to 1-(False Positives)                               |
| `Recall`               | 0.655  | True Positive Rate                                     |
| `Matthews Correlation` | 0.331  | Indicates weak predictive power – better than random prediction (0), but far from perfect (1) |
| `AUROC`                     | 0.713  | Moderately good predictions (perfect predictions = 1)   |

## Conclusion
While 20 may make biological sense as a cutoff from tissue specific to ubiquitous if we insist on a binary split, it would also make sense that there wouldn't be a huge different in sequence motifs in the promoter region for genes that are expressed in 20 tissues and genes that are expressed in 21 tissues. Rather, since the distribution has values at every integer between 1 and 29, wherever we place the 0/1 cutoff, there is likely to not be detectible differences across the closest values to this split.

