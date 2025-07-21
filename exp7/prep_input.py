import pandas as pd
from entex_express.utils import convert_to_csv_input_tpms

convert_to_csv_input_tpms("promoters_brain_substantia_nigra.fa", output_prefix="input", task="regression")
