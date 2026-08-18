library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

#DATASET <- "demo" # to run the demo
DATASET <- "GSE103001" # to run the analysis

#Load study genes as an R vector
study_genes <-scan(
    paste0("data/processed/", DATASET, "/study_genes.txt"),
    what = integer()
)

print (length(study_genes))

#Load background genes as an R vector
background_genes <-scan(
    paste0("data/processed/", DATASET, "/background_genes.txt"),
    what = integer()
)
print(length(background_genes))

#Run enrichGO()

EGO <- enrichGO(
    gene = study_genes,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)
print(EGO)

head(as.data.frame(EGO))

#Export results
EGO_df <- as.data.frame(EGO)
#Check if enrichment actually returned
if (nrow(EGO_df) == 0) {
    stop("No significant GO enrichment results were found")
}

write.csv(
    EGO_df,
    paste0("results/", DATASET, "/R_GO_enrichment.csv"),
    row.names = FALSE
)

print(dim(EGO_df))

#Create a dot plot
print("Creating dot plot")
png(
    paste0("figures/", DATASET, "/r_go_dotplot.png"),
    width = 1800,
    height = 1200,
    res = 150
)

dotplot(EGO, showCategory = 10, x = "GeneRatio", color = "p.adjust", font.size = 10)

dev.off()

print("Dot plot complete")

#Create a bar plot
print("Creating bar plot")

png(
    paste0("figures/", DATASET, "/r_go_barplot.png"),
    width = 1800,
    height = 1200,
    res = 150
)

barplot(EGO, showCategory = 10, font.size = 10)

dev.off()

print("Bar plot complete")

#Calculate term similarity for the enrichment map
print("Creating enrichment map")

EGO_sim <- pairwise_termsim(EGO)

png(
    paste0("figures/", DATASET, "/r_go_enrichment_map.png"),
    width = 1800,
    height = 1400,
    res = 150
)

emapplot(EGO_sim, showCategory = 10)

dev.off()

print ("Enrichment map complete")