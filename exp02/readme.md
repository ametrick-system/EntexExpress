# Experiment 2

## Summary
**Research Question:** Can the promoter sequence of a gene in the human reference genome (hg38) alone predict whether this gene is broadly expressed (i.e. ubiquitous, moderate, tissue specific) in the human body?

- **Input:** Promoter sequences of protein-coding genes, 2048 bp around the TSS
- **Label:** 0 (tissue-specific), 1 (moderately expressed), or 2 (ubiquitous)
- **Goal:** Test if general expression level correlates with promoter-intrinsic sequence features

## Raw Data
- Tissue data from Tissue_Specificity.zip ("the tissue specificity of gene expression and functional signals of cis-regulatory elements") downloaded from the EN-TEx portal. This experiment uses ``expressed_gene.tissue_specificity/pc.txt``, which contains a list of protein coding genes and how many tissues they are expressed in
- Human genome sequence from GENCODE Human Genome Release 48 primary assembly (GRCh38)

## Input
- ``extract_pc_promoters.py`` extracts the promoter region for each gene in pc.txt (1024 nucleotides upstream and 1024 nucleotides downstream) and stores results in ``promoters.bed`` with ternary labels indicating whether that gene is expressed in ≤ 15 tissues (0), > 15 and < 29 tissues (1), and = 29 tissues (2), with label counts 5183 0's, 4137 1's, and 6680 2's
- ``promoters.fa`` contains the actual promoter sequences from GRCh38 based on windows in ``promoters.bed`` (generated with bedtools command ``bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed promoters.bed -s -name -fo promoters.fa``)
- ``generate_input.py`` creates the 3 input CSV files: ``train.csv`` (80%), ``test.csv`` (10%), and ``dev.csv`` (10%) to be passed to DNABERT2 in finetuning, where each line contains a promoter sequence and a 0, 1 or 2 indicating tissue specificity of that gene as described above
-  ``train.csv``, ``test.csv``, and ``dev.csv`` are all stored in the ``input`` directory

## Fine-tuning DNABERT2
The script ``finetune.sh`` runs the DNABERT2 fine tuning with the default weights, steps, epochs, and learning rate. To run, first make the file executable (``chmod +x finetune.sh``) then run (``./finetune.sh``). All the terminal output from the fine-tuning process is saved in ``finetune.log``.

## Results
Here is a table that presents the results stored in ``output/results/entex_express_exp02/eval_results.json`` (note that each value is rounded to the nearest hundredth). Note that we also include the results from the control version of Experiment 2 here in order to compare.

| Metric     | exp02 | exp02_control |
|------------|--------|----------------|
| Loss       | 1.068  | 1.098          |
| Accuracy   | 0.483  | 0.339          |
| F1         | 0.348  | 0.169          |
| Precision  | 0.317  | 0.113          |
| Recall     | 0.404  | 0.333          |
| MC         | 0.152  | 0.000          |

## Conclusion
While certainly better than random, splitting into 3 bins with these cutoffs yields low predictive power. 

