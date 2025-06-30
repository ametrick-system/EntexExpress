import pandas as pd
import matplotlib.pyplot as plt
import os

bed_path = "promoters.bed"
bed_cols = ["chrom", "start", "end", "label", "name", "strand", "tissue_count"]

if not os.path.exists(bed_path):
    raise FileNotFoundError(f"Cannot find {bed_path}")

bed = pd.read_csv(bed_path, sep="\t", names=bed_cols)

# Plot histogram of tissue counts
plt.figure(figsize=(8, 5))
plt.hist(bed["tissue_count"], bins=30, edgecolor="black")
plt.title("Distribution of Tissue Counts per Gene")
plt.xlabel("Number of Tissues")
plt.ylabel("Number of Genes")
plt.grid(True)
plt.tight_layout()
plt.savefig("histogram_tissue_counts.png")
plt.close()

# Assign bin cutoffs
bins = [0, 15, 28, 29]
labels = [0, 1, 2]
bed["expression_bin"] = pd.cut(bed["tissue_count"], bins=bins, labels=labels, include_lowest=True)
print(bed["expression_bin"].value_counts(sort=False))


# Plot histogram again with bin cutoffs
plt.figure(figsize=(8, 5))
plt.hist(bed["tissue_count"], bins=30, edgecolor="black")
plt.axvline(x=15, color='red', linestyle='--', label='Tissue-specific ↔ Intermediate')
plt.axvline(x=28, color='green', linestyle='--', label='Intermediate ↔ Ubiquitous')
plt.title("Distribution of Tissue Counts with Bins")
plt.xlabel("Number of Tissues")
plt.ylabel("Number of Genes")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("histogram_with_bins.png")
plt.close()

print("Done! Plots saved:")
print("- histogram_tissue_counts.png")
print("- histogram_with_bins.png")