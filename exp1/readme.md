# Experiment 1

## Summary
**Research Question:** How well can DNABERT2 predict whether a protein-coding gene is broadly expressed (i.e., active in many tissues) based solely on its promoter sequence?

- **Input:** Promoter sequences of protein-coding genes
- **Label:** Binary (1 if gene expressed in ≥ 20 tissues, else 0)
- **Goal:** Test if general expression level correlates with promoter-intrinsic sequence features

For setup, simply activate the virtual environment: ``conda activate dnabert2``

## Raw Data
- Tissue data from Tissue_Specificity.zip ("the tissue specificity of gene expression and functional signals of cis-regulatory elements") downloaded from the EN-TEx portal. This experiment uses ``expressed_gene.tissue_specificity/pc.txt``, which contains a list of protein coding genes and how many tissues they are expressed in
- Human genome sequence from GENCODE Human Genome Release 48 primary assembly (GRCh38)

## Input
- ``extract_pc_promoters.py`` extracts the promoter region for each gene in pc.txt and stores results in ``promoters.bed``
- ``promoters.fa`` contains the actual promoter sequences from GRCh38 based on windows in ``promoters.bed`` (generated with bedtools command ``bedtools getfasta -fi ~/LargeFiles/GRCh38.primary_assembly.genome.fa -bed promoters.bed -s -name -fo promoters.fa``)
- ``prep_input.py`` creates the 3 input CSV files: ``train.csv`` (80%), ``test.csv`` (10%), and ``dev.csv`` (10%) to be passed to DNABERT2 in finetuning, each line contains a promoter sequence and a 1 or 0 indicating whether that gene is expressed in ≥ 20 tissues
-  ``train.csv``, ``test.csv``, and ``dev.csv`` are all stored in the ``input`` directory

## Fine-tuning DNABERT2
The script ``finetune.sh`` runs the DNABERT2 fine tuning with the default weights, steps, epochs, and learning rate. To run, first make the file executable (``chmod +x finetune.sh``) then run (``./finetune.sh``). All the output from the fine-tuning process is saved in ``finetune.log``.

## Results
Here is a table that presents the results stored in ``output/results/DNABERT2_exp1_promoters/eval_results.json`` (note that each value is rounded to the nearest hundredth).

| Metric                      | Value  | Notes                                                                 |
|-----------------------------|--------|-----------------------------------------------------------------------|
| `eval_loss`                 | 0.658  | Slightly better than random (~0.693)                                  |
| `eval_accuracy`             | 0.632  | Also better than random (~0.5)                                        |
| `eval_f1`                   | 0.387  | Low — indicates class imbalance or poor precision/recall balance      |
| `eval_precision`            | 0.316  | Many predicted positives were incorrect                               |
| `eval_recall`               | 0.500  | Half of true positives were found                                     |
| `eval_matthews_correlation` | 0.000  | Suggests predictions are not better than random chance from a correlation perspective (more sensitive than accuracy) |

## Conclusion
The F1 result suggested that the data might not have been balanced. Sure enough (as discovered in ``interpret_results.py``), there were 7946 1's and 4854 0's in ``train.csv``, indicating that 20 is not a good threshold for separating between 1/0. 


