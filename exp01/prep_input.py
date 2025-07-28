from entex_express.utils import convert_to_csv_input

convert_to_csv_input("promoters.fa", "promoters.bed", output_prefix="input", task="binary", threshold=20)