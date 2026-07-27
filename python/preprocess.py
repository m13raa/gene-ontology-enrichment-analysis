import pandas as pd
print ("start go preprocessing data")
# Load gene2go data by creating a variable pointing to the relative file path
# Seperate the .gz file by tab character
gene2go = pd.read_csv("data/raw/gene2go.gz", sep="\t")
#Analyze the dataframe
print(gene2go.head())
print(gene2go.shape)
print(gene2go.columns)
print(gene2go.info())
#Count the number of different species in the dataframe
print(gene2go["#tax_id"].nunique())
#How many of the species are the human species (homo sapiens?) Taxonomy ID is 9606 for human species
homo_sapiens_gene2go = gene2go[gene2go["#tax_id"]==9606]
print(homo_sapiens_gene2go.shape)
# print(homo_sapiens_gene2go ["#tax_id"].unique()) Check the previous filtering
#How many human genes have GO annotations within the gene2go dataframe?
print(homo_sapiens_gene2go["GeneID"].nunique())
#Create a dictionary mapping GeneIDs to their GO annotations
gene_2_go_ditionary=(homo_sapiens_gene2go.groupby("GeneID")["GO_ID"].apply(list).to_dict())
print(len(gene_2_go_ditionary))
print(list(gene_2_go_ditionary.items())[:5])
#Clean up the GO terms, repeated GO ID's not necessary as they prove the same gene-to-go connection through different forms of evidence
gene_2_go_ditionary=(homo_sapiens_gene2go.groupby("GeneID")["GO_ID"].apply(lambda x:sorted(set(x)))).to_dict()
