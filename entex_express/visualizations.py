# Must be run while in dnabert2 conda virtual environment

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import torch
import numpy as np
from sklearn.metrics import roc_curve, auc
from transformers import AutoTokenizer, AutoModelForSequenceClassification, Trainer
from .dnabert2_train_regression import SupervisedDataset, DataCollatorForSupervisedDataset

# For binary classification tasks
def plot_roc_curve(checkpoint_path, data_path):
    tokenizer = AutoTokenizer.from_pretrained(checkpoint_path, trust_remote_code=True)
    model = AutoModelForSequenceClassification.from_pretrained(checkpoint_path, trust_remote_code=True)

    test_dataset = SupervisedDataset(
        data_path=f"{data_path}/test.csv",
        tokenizer=tokenizer,
        kmer=-1,
        task="classification"
    )

    data_collator = DataCollatorForSupervisedDataset(
        tokenizer=tokenizer,
        regression=False
    )

    # Make predictions
    trainer = Trainer(model=model, tokenizer=tokenizer, data_collator=data_collator)
    predictions_output = trainer.predict(test_dataset)

    logits = predictions_output[0] if isinstance(predictions_output, tuple) else predictions_output.predictions
    if isinstance(logits, tuple):
        logits = logits[0]

    # Convert to numpy array
    if isinstance(logits, list):
        logits = np.stack([np.array(x) for x in logits])
    else:
        logits = np.array(logits)

    if logits.ndim != 2 or logits.shape[1] != 2:
        raise ValueError(f"Expected logits to be shape (N, 2), got shape {logits.shape}")

    # Softmax
    exps = np.exp(logits - np.max(logits, axis=1, keepdims=True))  # stability
    probs = exps / np.sum(exps, axis=1, keepdims=True)
    positive_probs = probs[:, 1]

    labels = np.array(test_dataset.labels)

    # ROC & AUROC
    fpr, tpr, _ = roc_curve(labels, positive_probs)
    roc_auc = auc(fpr, tpr)

    # Plot
    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, label=f"AUROC = {roc_auc:.3f}", color="blue")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("roc_curve.png", dpi=300)
    print("Saved ROC curve as roc_curve.png")


def plot_regression_results_scatter(checkpoint_path, data_path):
    tokenizer = AutoTokenizer.from_pretrained(checkpoint_path, trust_remote_code=True)
    model = AutoModelForSequenceClassification.from_pretrained(checkpoint_path, trust_remote_code=True)

    test_dataset = SupervisedDataset(
        data_path=f"{data_path}/test.csv",
        tokenizer=tokenizer,
        kmer=-1,
        task="regression"
    )

    data_collator = DataCollatorForSupervisedDataset(
        tokenizer=tokenizer,
        regression=True
    )

    trainer = Trainer(model=model, tokenizer=tokenizer, data_collator=data_collator)

    # Make predictions
    predictions_output = trainer.predict(test_dataset)

    if isinstance(predictions_output, tuple):
        predictions = predictions_output[0]  # first element = logits
    else:
        predictions = predictions_output.predictions

    # If model outputs a tuple (logits, extra)
    if isinstance(predictions, tuple):
        predictions = predictions[0]

    preds = np.array(predictions).squeeze()
    labels = np.array(test_dataset.labels)

    # Plot predicted vs actual
    plt.figure(figsize=(6, 6))
    plt.scatter(labels, preds, alpha=0.5)
    plt.plot([labels.min(), labels.max()], [labels.min(), labels.max()], 'r--', label="Perfect Prediction")
    plt.xlabel("Actual Values")
    plt.ylabel("Predicted Values")
    plt.title("Predicted vs Actual (Regression)")
    plt.legend()

    plt.savefig("pred_vs_actual.png", dpi=300)
    print("Saved plot as pred_vs_actual.png")

def plot_histogram_from_bed(save_path, bed_path, bed_cols, plot_col, fig_name, bins_list=None, labels=None):
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"Cannot find {bed_path}")

    bed = pd.read_csv(bed_path, sep="\t", names=bed_cols)

    # Plot histogram of plot_col
    plt.figure(figsize=(8, 5))
    plt.hist(bed[plot_col], bins=30, edgecolor="black")
    plt.title(f"Distribution of {plot_col}")
    plt.xlabel(f"{plot_col}")
    plt.ylabel("Number of Genes")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{fig_name}.png")
    plt.close()

    print(f"{fig_name}.png saved in {save_path}")

    # Also plot histogram with bin cutoffs marked, if given
    if bins_list is not None:
        bed["bin"] = pd.cut(bed[plot_col], bins=bins_list, labels=labels, include_lowest=True)

        plt.figure(figsize=(8, 5))
        plt.hist(bed[plot_col], bins=30, edgecolor="black")

        colors = list(mcolors.CSS4_COLORS.keys())

        for i in range(1,len(bins_list)):
            bin_x = bins_list[i]
            plt.axvline(x=bin_x, color=colors[i % len(colors)], linestyle='--', label=labels[i-1])
       
        plt.title(f"Distribution of {plot_col} with Bins")
        plt.xlabel(f"{plot_col}")
        plt.ylabel("Number of Genes")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{fig_name}_with_bins.png")
        plt.close()

        print(f"{fig_name}_with_bins.png saved in {save_path}")
