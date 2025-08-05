# Must run from alphagenome-env conda virtual environment

from entex_express.utils import config_dnabert2_input

config_dnabert2_input(
    fasta='promoters.fa',
    label_key='label',
    save_prefix='input',
    task='classification',
    cutoffs=[0, 1, 2]
)

