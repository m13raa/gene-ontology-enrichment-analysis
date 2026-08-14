import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

results = pd.read_csv("results/goatools_enrichment_results.csv")
print(results.head())
print("\nShape:", results.shape)

#log the FDR values for better visualization
results["neg_log10_FDR"] = -np.log10(results["FDR"])
print(results[["Name", "FDR", "neg_log10_FDR"]].head())

#select the top 10 significant GO terms based on FDR
top_terms = (results.sort_values("FDR").head(10).copy())
print("\nTop 10 GO terms:")
print(
    top_terms[["GO_ID", "Name", "FDR", "neg_log10_FDR"]
    ]
)
#create enrichment plot
plt.figure(figsize=(10,7))

plt.barh(top_terms["Name"], top_terms["neg_log10_FDR"])

plt.xlabel("-log10(FDR)")
plt.ylabel("GO Term")

plt.title("Top 10 Significantly Enriched GO Terms")
plt.gca().invert_yaxis()  # Invert y-axis to have the most significant term at the top
plt.tight_layout() #Don't have the labels be cut off

#save figure before plotting
plt.savefig("figures/top_10_go_terms_enrichment.png", dpi=300, bbox_inches='tight')

plt.show()

#calculate the fold enrichment
results["Fold_Enrichment"] = (results["Study_Count"]/results["Study_Total"])/(results["Population_Count"]/results["Population_Total"])

#check
print(results[["Name", "FDR", "Fold_Enrichment"]].head())

#select top 10 significant GO terms based on fold enrichment

top_enriched_terms = (results.sort_values("Fold_Enrichment", ascending=False).head(10).copy())

#create fold enrichment plot
plt.figure(figsize=(10,7))
plt.barh(top_enriched_terms["Name"], top_enriched_terms["Fold_Enrichment"])

plt.xlabel("Fold Enrichment")
plt.ylabel("GO Term")

plt.title("Top 10 GO Terms by Fold Enrichment")
plt.gca().invert_yaxis()  # Invert y-axis to have the most enriched term at the top
plt.tight_layout() #Don't have the labels be cut off")")

#save figure before plotting
plt.savefig("figures/top_10_go_terms_fold_enrichment.png", dpi=300, bbox_inches='tight')

plt.show()