import pandas as pd

# Load GOATools results
python_results = pd.read_csv("results/goatools_enrichment_results.csv")
# Load clusterProfiler results
r_results = pd.read_csv("results/R_GO_enrichment.csv")
# Print shape of results
print("GOATools results:", python_results.shape)
print("clusterProfiler results:", r_results.shape)

# Make sets out of the python and r GO ID's
python_go_ids = set(python_results["GO_ID"])
r_go_ids = set(r_results["ID"])
print ("\nUnique GOATools terms:")
print(len(python_go_ids))
print("\nUnique clusterProfiler terms:")
print(len(r_go_ids))

#Find the overlapping GO terms
shared_go_ids = python_go_ids & r_go_ids
print("\nShared significant GO terms:")
print(len(shared_go_ids))

#Find terms found using python only
python_only = python_go_ids - r_go_ids
#Find terms found using r only
r_only = r_go_ids - python_go_ids
print("\nGOATools-only terms:")
print(len(python_only))
print("\nclusterProfiler-only terms")
print(len(r_only))

#Calculate the percentages of overlap
python_overlap_pct = (len(shared_go_ids)/len(python_go_ids)*100)
r_overlap_pct = (len(shared_go_ids)/len(r_go_ids)*100)

print("\nPercentage of GOATools terms also found by clusterProfiler",
      round(python_overlap_pct,2),"%")
print("\nPercentage of clusterProfiler terms also found by GOATools:",
      round (r_overlap_pct, 2),"%")

#clusterProfiler returned more items, hence the smaller percentage compared to the GOATools perspective

#Create a table
shared_python = python_results[python_results["GO_ID"].isin(shared_go_ids)].copy()
shared_r = r_results[r_results["ID"].isin(shared_go_ids)].copy()
#isin  creates true/false values, asking if the GO ID is in the shared set

#Rename the columns to minimize confusion
shared_r = shared_r.rename(
    columns={
        "ID":"GO_ID",
        "Description": "R_Name",
        "p.adjust": "R_Adjusted_P"
    }
)
shared_python = shared_python.rename(
    columns={
        "Name": "Python_Name",
        "FDR": "Python_FDR"
    }
)

#Merge to create a table and compare
comparison = pd.merge(
    shared_python,
    shared_r[
        [
            "GO_ID",
            "R_Name",
            "R_Adjusted_P"
        ]
    ],
    on="GO_ID",
    how="inner"
)

#Save the comparison table
comparison = comparison.sort_values(
    "Python_FDR"
)
comparison.to_csv(
    "results/python_vs_r_shared_go_terms.csv",
    index=False
)
print("\nShared comparison table saved.")

#Save GO terms specific to the package
python_only_df = python_results[python_results["GO_ID"].isin(python_only)]
python_only_df.to_csv("results/goatools_only_terms.csv", index=False)
r_only_df = r_results[r_results["ID"].isin(r_only)]
r_only_df.to_csv("results/clusterProfiler_only_terms.csv", index=False)

#Create a comparison visual
import matplotlib.pyplot as plt

comparison_counts = {
    "Shared": len(shared_go_ids),
    "GOATools only": len(python_only),
    "clusterProfiler only": len(r_only)
}

plt.figure(figsize=(8,5))
plt.bar(
    comparison_counts.keys(),
    comparison_counts.values()
)
plt.ylabel("Number of Significant GO Terms")
plt.title("Comparison of Significant GO Terms")

plt.tight_layout()

plt.savefig(
    "figures/python_vs_r_go_term_comparison.png", 
    dpi=300,
    bbox_inches="tight"
)
plt.show()