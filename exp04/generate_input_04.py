from entex_express.utils import config_dnabert2_input

config_dnabert2_input(
    fasta='promoters.fa',
    label_key='num_tissues',
    save_prefix='input',
    task='classification',
    cutoffs=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29] # so gene expressed in 1 tissue will have label "0"
)
