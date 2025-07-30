from entex_express.utils import config_dnabert2_input

config_dnabert2_input(
    fasta='promoters.fa',
    label_key='num_tissues',
    save_prefix='input',
    task='regression'
)
