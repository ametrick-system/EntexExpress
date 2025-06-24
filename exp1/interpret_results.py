import pandas as pd

df = pd.read_csv("input/train.csv")
print(df["label"].value_counts())

# results: 7946 1's, 4854 0's