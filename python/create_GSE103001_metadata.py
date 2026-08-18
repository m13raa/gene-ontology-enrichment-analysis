import pandas as pd

# Patient identifiers from GSE103001
patients = [
    "12-02",
    "12-03",
    "13-02",
    "13-03",
    "13-05",
    "13-07",
    "13-10",
    "13-11",
    "13-13",
    "13-17",
    "14-25",
    "14-31",
    "14-33",
    "14-34",
    "14-35",
    "14-36",
    "14-38",
    "14-41",
    "15-01",
    "15-02",
    "15-03",
    "15-04",
]

normal_gsm = [f"GSM{i}" for i in range(2752350, 2752372)]
tumor_gsm = [f"GSM{i}" for i in range(2752372, 2752394)]

metadata = pd.DataFrame({
    "GSM": normal_gsm + tumor_gsm,
    "patient": patients + patients,
    "condition": ["normal"] * 22 + ["tumor"] * 22
})

print(metadata)
print("\nSample counts:")
print(metadata["condition"].value_counts())

print("\nNumber of patients:", metadata["patient"].nunique())

metadata.to_csv(
    "data/processed/GSE103001/sample_metadata.csv",
    index=False
)

print(
    "\nSaved metadata to "
    "data/processed/GSE103001/sample_metadata.csv"
)