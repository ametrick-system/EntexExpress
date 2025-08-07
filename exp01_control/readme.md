# Experiment 1 CONTROL

## Summary
**Research Question:** Can the promoter sequence of a gene in the human reference genome (hg38) alone predict whether this gene is broadly expressed (i.e. ubiquitous, moderate, tissue specific) in the human body?

- **Input:** Promoter sequences of protein-coding genes, 2048 bp around the TSS
- **Label:** Assigned random binary labels with equal probability of 0/1
- **Goal:** Test if general expression level correlates with promoter-intrinsic sequence features

## Raw Data
- ``promoters.fa`` exactly from Experiment 1

## Input and Finetuning
- ``finetune.sh`` first generates the 3 input CSV files ``train.csv`` (80%), ``test.csv`` (10%), and ``dev.csv`` (10%) to be passed to DNABERT2 in finetuning, where each line contains a promoter sequence and a 1 or 0 label assigned randomly with equal probably
- ``finetune.sh`` then runs the DNABERT2 fine tuning with the default weights, steps, epochs, and learning rate
- To run, first make the file executable (``chmod +x finetune.sh``) then run (``./finetune.sh``). All the terminal output from the fine-tuning process is saved in ``finetune.log``.

## Results
Here is a table that presents the results stored in ``output/results/entex_express_exp01_control/eval_results.json`` (note that each value is rounded to the nearest hundredth). Additionally, ``plot_roc.py`` plots the Receiver Operating Characteristic (ROC) curve in ``roc_curve.png`` and computes the AUROC.

| Metric                      | Value  |                                                               |
|-----------------------------|--------|
| `Loss`                 | 0.694  |
| `Accuracy`             | 0.485  |
| `F1`                   | 0.327  |
| `Precision`            | 0.243  |
| `Recall`               | 0.500  | True Positive Rate                                     |
| `Matthews Correlation` | 0.000  |
| `AUROC`                     | 0.501  |

## Conclusion
The control experiment performed as expected, indicating the interpretation of the results from Experiment 1 are valid.

