import pandas as pd 
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
# DATASET = "demo" #to run the demo 
DATASET = "GSE103001" #to run the real analysis

# Load the study genes
with open(f"data/processed/{DATASET}/study_genes.txt") as f:
    study_genes = [
        int(line.strip())
        for line in f
        if line.strip()
    ]
# Load the background genes
with open(f"data/processed/{DATASET}/background_genes.txt") as f:
    background_genes = [
        int(line.strip())
        for line in f
        if line.strip()         
    ]   

print("Number of study genes:", len(study_genes))
print("Number of background genes:", len(background_genes))

#See if every study gene is in the background gene list
missing_genes = set(study_genes) - set(background_genes)
print("Study genes missing from background", missing_genes)
#missing_genes should return an empty set if all study genes are in the background gene list

# Load the NCBI Gene Ontology annotations
gene2go = pd.read_csv("data/raw/gene2go", sep="\t")
# Filter to keep the Homo sapiens annotations
human_gene2go = gene2go[gene2go["#tax_id"]==9606]
print("Human GO annotations:", human_gene2go.shape)
#Build GeneID to GO mapping
gene_to_go=(human_gene2go.groupby("GeneID")["GO_ID"].apply(set).to_dict())
# Load GO DAG
obo_file = "data/raw/go-basic.obo"
go_dag = GODag(obo_file)
goea = GOEnrichmentStudy(
    background_genes,
    gene_to_go,
    go_dag,
    propagate_counts=True,
    alpha=0.05,
    methods=["fdr_bh"]
)
# Run enrichment analysis
go_results = goea.run_study(study_genes)
print("\nNumber of GO terms test:")
print(len(go_results))

# Keep the statistically significant GO terms per FDR
significant_results=[
    result
    for result in go_results
    if result.p_fdr_bh < 0.05
]
print("\nNumber of significant GO terms (FDR < 0.05):")
print(len(significant_results))

# Look at the first few significant terms
for result in significant_results[:10]:
    print(
        "\nGO term:", result.GO, #GO ID
        "\nName:", result.name, #biological term name
        "\nNamespace:", result.NS, #ontology namespace (BP, CC, MF)
        "\nEnrichment:",result.enrichment, #enriched or purified
        "\nStudy count:", result.study_count, #number of study genes annotated to this GO term
        "\nStudy total:", result.study_n, #total number of study genes
        "\nPopulation count:", result.pop_count, #number of background genes annotated to this GO term
        "\nPopulation total:", result.pop_n, #total number of background genes
        "\nRaw p-value:",result.p_uncorrected, #raw p-value
        "\nFDR:",result.p_fdr_bh #Benjamini-Hochberg corrected p-value
    )  

# Export the results to a CSV file

results_df = pd.DataFrame([
    {
        "GO_ID": result.GO,
        "Name": result.name,
        "Namespace": result.NS,
        "Enrichment": result.enrichment,
        "Study_Count": result.study_count,
        "Study_Total": result.study_n,
        "Population_Count": result.pop_count,
        "Population_Total": result.pop_n,
        "P_Value": result.p_uncorrected,
        "FDR": result.p_fdr_bh
    }
    for result in significant_results
])

# Check that the columns were created correctly
print("\nResult columns:")
print(results_df.columns)

# Sort strongest results first
results_df = results_df.sort_values(
    by="FDR",
    ascending=True
)

# Save results
results_df.to_csv(
    f"results/{DATASET}/goatools_enrichment_results.csv",
    index=False
)

print(f"\nResults saved to results/{DATASET}/goatools_enrichment_results.csv")