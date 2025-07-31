import pandas as pd
from entex_express.utils import convert_to_csv_input_random_labels

convert_to_csv_input_random_labels("promoters.fa", "promoters.bed", output_prefix="input", labels=[0, 1, 2])

# ensure input is reasonably balanced
df = pd.read_csv("input/train.csv")
print(df["label"].value_counts())

'''
2    4400
1    4215
0    4185
'''