from pathlib import Path
import pandas as pd

#Look into the folder that contains all the 44 HTSeq files
counts_dir = Path("data/raw/GSE103001/counts")

#Find all the .txt.gz files
count_files = sorted(counts_dir.glob("*.txt.gz"))
print("Number of count files:", len(count_files))
if len(count_files) != 44:
    raise ValueError(
        f"Expected 44 count files, found {len(count_files)}"
    )

samples = []

#Clean up the names on file
for file in count_files:
    GSM = file.name.split("_")[0]
    sample = pd.read_csv(
        file,
        sep=r"\s+",
        header=None,
        names=["GeneID", GSM]
    )

    #Remove the HTSeq summary rows:
    sample = sample[
        ~sample["GeneID"].str.startswith("__")
    ]

    # Make Ensembl IDs the row names
    sample = sample.set_index("GeneID")
    samples.append(sample)

#Join all the samples by GeneID
count_matrix = pd.concat(
    samples,
    axis=1
)

print("\nCount matrix dimensions:")
print(count_matrix.shape)
print("\nMissing values:")
print(count_matrix.isna().sum().sum())

#Validate the count matrix against the metadata
metadata = pd.read_csv("data/processed/GSE103001/sample_metadata.csv")

matrix_samples = set(count_matrix.columns)
metadata_samples = set(metadata["GSM"])

print ("\nSamples missing from matrix:")
print(metadata_samples - matrix_samples)

print("\nSamples missing from metadata:")
print(matrix_samples - metadata_samples)

#Save
count_matrix.to_csv(
    "data/processed/GSE103001/count_matrix.csv"
)
print ("\nCount matrix saved.")