import pandas as pd
import matplotlib.pyplot as plt
import os

bed_path = "promoters_brain_substantia_nigra.bed"
bed_cols = ["chr", "start", "end", "name", "tpm", "log_tpm", "strand"]

if not os.path.exists(bed_path):
    raise FileNotFoundError(f"Cannot find {bed_path}")

bed = pd.read_csv(bed_path, sep="\t", names=bed_cols)

# Plot histogram of log(tpm)
plt.figure(figsize=(8, 5))
plt.hist(bed["log_tpm"], bins=30, edgecolor="black")
plt.title("Distribution of log(TPM) per Gene")
plt.xlabel("log(TPM)")
plt.ylabel("Number of Genes")
plt.grid(True)
plt.tight_layout()
plt.savefig("histogram_log_tpm.png")
plt.close()

# Assign bin cutoffs
bins = [0, 1, 2.5, 12]
labels = [0, 1, 2]
bed["expression_bin"] = pd.cut(bed["log_tpm"], bins=bins, labels=labels, include_lowest=True)
print(bed["expression_bin"].value_counts(sort=False))

# Plot histogram again with bin cutoffs
plt.figure(figsize=(8, 5))
plt.hist(bed["log_tpm"], bins=30, edgecolor="black")
plt.axvline(x=1, color='red', linestyle='--')
plt.axvline(x=2.5, color='green', linestyle='--')
plt.title("Distribution of log(TPM) with Bins")
plt.xlabel("log(TPM)")
plt.ylabel("Number of Genes")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("histogram_with_bins.png")
plt.close()

print("Done! Plots saved:")
print("- histogram_log_tpm.png")
print("- histogram_with_bins.png")