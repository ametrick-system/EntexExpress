import pandas as pd
from entex_express.utils import convert_to_csv_input

convert_to_csv_input("promoters.fa", "promoters.bed", output_prefix="input", task="binary", threshold=27)

# ensure balanced input
df = pd.read_csv("input/train.csv")
print(df["label"].value_counts())

# result: 6501 0's, 6299 1's