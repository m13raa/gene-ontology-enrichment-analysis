import pandas as pd
# Load the gene2go annotation file from NCBI
gene2go = pd.read_csv("data/raw/gene2go", sep="\t")
# Isolate and keep the human annotations
human_gene2go = gene2go[gene2go["#tax_id"] == 9606]
# Get list of human Gene ID
background_genes = human_gene2go["GeneID"].unique()
print("Number of human genes in the background:", len(background_genes))
#Save the background genes to a text file
pd.Series(background_genes).to_csv("data/processed/background_genes.txt", index=False, header=False)
print("Background genes saved")
