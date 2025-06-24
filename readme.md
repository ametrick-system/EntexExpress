**Research Question:** Can DNA sequence alone (promoter or other cis-regulatory regions) predict whether a gene is expressed in a given tissue, and at what level?

## Environment Setup 

```bash
srun --partition=gpu_devel --gres=gpu:1 --cpus-per-task=2 --mem=8G --time=01:00:00 --pty bash
conda activate dnabert2
```