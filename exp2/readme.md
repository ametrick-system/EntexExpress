# Experiment 2

## Summary
For Experiment 2, we simply replicated Experiment 1 but with balanced 1's and 0's in the input. This meant a threshold of 27, which led to 6501 0's and 6299 1's in ``train.csv``.

## Results
Here is a table that presents the results stored in ``output/results/DNABERT2_exp2_promoters/eval_results.json`` (note that each value is rounded to the nearest hundredth).

| Metric                      | Value  | Notes                                                                 |
|-----------------------------|--------|-----------------------------------------------------------------------|
| `eval_loss`                 | 0.662  | Slightly better than random (~0.693) but slightly worse than Experiment 1 (0.658)                                  |
| `eval_accuracy`             | 0.615  | Better than random (~0.5), worse than Experiment 1 (0.632)                                        |
| `eval_f1`                   | 0.606  | Moderate performance, much better than Experiment 1 (0.387) |
| `eval_precision`            | 0.636  | Moderate performance, much better than Experiment 1 (0.316)                               |
| `eval_recall`               | 0.621  | Moderate performance, more than half of TPs found, better than Experiment 1 (0.5)                                     |
| `eval_matthews_correlation` | 0.256  | Moderate signal |

## Conclusion
DNABERT2 learned moderately well. However, the results reflect the fact that binary classification is not a natural measure for level of gene expression. Experiment 3 will divide the input into closer bins.