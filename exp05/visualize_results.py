import torch
import numpy as np
import matplotlib.pyplot as plt
from transformers import AutoTokenizer, AutoModelForSequenceClassification, Trainer
from dnabert2_train_regression import SupervisedDataset, DataCollatorForSupervisedDataset

# ---------------------------
# 1. Load from checkpoint
# ---------------------------
checkpoint_path = "/home/asm242/EntexExpress/exp05/output/checkpoint-16000"
tokenizer = AutoTokenizer.from_pretrained(checkpoint_path, trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained(checkpoint_path, trust_remote_code=True)

# ---------------------------
# 2. Load the test dataset
# ---------------------------
data_path = "/home/asm242/EntexExpress/exp05/input"
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

# ---------------------------
# 3. Make predictions
# ---------------------------
predictions_output = trainer.predict(test_dataset)

# Unpack depending on Hugging Face version
if isinstance(predictions_output, tuple):
    predictions = predictions_output[0]  # first element = logits
else:
    predictions = predictions_output.predictions

# If model outputs a tuple (logits, extra)
if isinstance(predictions, tuple):
    predictions = predictions[0]

# Now squeeze properly
preds = np.array(predictions).squeeze()

labels = np.array(test_dataset.labels)

# ---------------------------
# 4. Plot Predicted vs Actual
# ---------------------------
plt.figure(figsize=(6, 6))
plt.scatter(labels, preds, alpha=0.5)
plt.plot([labels.min(), labels.max()], [labels.min(), labels.max()], 'r--', label="Perfect Prediction")
plt.xlabel("Actual Values")
plt.ylabel("Predicted Values")
plt.title("Predicted vs Actual (Regression)")
plt.legend()

# Save to file
plt.savefig("pred_vs_actual.png", dpi=300)
plt.show()

print("Saved plot as pred_vs_actual.png")

