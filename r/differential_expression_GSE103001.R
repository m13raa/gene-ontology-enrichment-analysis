#Install neccesary libraries to conduct DEA
#BiocManager::install("DESeq2")

library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

counts_matrix <- read.csv(
    "data/processed/GSE103001/count_matrix.csv",
    row.names = 1,
    check.names = FALSE
)
metadata <- read.csv(
    "data/processed/GSE103001/sample_metadata.csv"
)
rownames(metadata) <- metadata$GSM

metadata <- metadata[
    colnames(counts_matrix),
]
stopifnot(
    all(
        rownames(metadata) ==
        colnames(counts_matrix)
    )
)

#Define the experimental variables

metadata$patient <- factor(
    metadata$patient
)

metadata$condition <- factor(
    metadata$condition,
    levels = c(
        "normal",
        "tumor"
    )
)

#Paired DESeq2 design
dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = metadata,
    design = ~ patient + condition #Each patient has tumor + normal data, so use patient + data to look for patient differences within before estimating the effect of the tumor
)

#Remove the low count genes to reduce noise and make the analysis more stable
keep <- rowSums(
    counts(dds) >= 10
) >=5 #Keep the gene if it has at least 10 counts in at least 5 samples
dds <- dds[keep,]
cat("Genes after low count filtering:", nrow(dds),"\n")

#Run DESeq2. Asking for tumor relativeness to normal
dds <-DESeq(dds)
res <-results(
    dds,
    contrast = c(
        "condition",
        "tumor",
        "normal"
    ), alpha = 0.05
)

#Results table
res_df <-as.data.frame(res)

res_df$ENSEMBL <-rownames(
    res_df 
)
res_df <- res_df[
    !is.na(res_df$padj),
]
res_df <-res_df[
    order(res_df$padj),
]

# Select the significantly differentially expressed genes
sig <- res_df[
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
]
cat(
    "Significant DE genes:",
    nrow(sig),
    "\n"
)
cat(
    "Upregulated:",
    sum(sig$log2FoldChange >= 1),
    "\n"
)
cat(
    "Downregulated:",
    sum(sig$log2FoldChange <= -1),
    "\n"
)
#Save the results table
write.csv(
    res_df,
    "results/GSE103001/differential_expression_all.csv",
    row.names = FALSE
)

#Count files use ENSEMBL IDs, so convert them to ENTREZ ID's
mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(
        res_df$ENSEMBL
    ),
    keytype = "ENSEMBL",
    columns = c(
        "ENTREZID",
        "SYMBOL"
    )
)
#Remove the unmapped entries
mapping <- mapping[
    !is.na(mapping$ENTREZID),
]
#Map the sig genes
sig_mapped <- merge(
    sig,
    mapping,
    by = "ENSEMBL"
)

#Save the mapped table
write.csv(
    sig_mapped,
    "results/GSE103001/differential_expression_significant_mapped.csv",
    row.names = FALSE
)

#Generate the gene set to conduct the enrichment analysis on
study_genes <- unique(
    sig_mapped$ENTREZID
)

writeLines(
    study_genes,
    "data/processed/GSE103001/study_genes.txt"
)

cat(
    "Study genes:",
    length(study_genes),
    "\n"
)

#Generate the RNA-seq background. Sharpen it so that all genes included in the background are ones that were eligible for DE testing # nolint
background_ensembl <- rownames(
    dds
)
background_mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = background_ensembl,
    keytype = "ENSEMBL",
    columns = "ENTREZID"
)
background_mapping <- background_mapping[
    !is.na(
        background_mapping$ENTREZID
    ),
]
background_genes <- unique(
    background_mapping$ENTREZID
)
writeLines(
    background_genes,
    "data/processed/GSE103001/background_genes.txt"
)
#Validate the background txt
cat(
    "Background genes:",
    length(background_genes),
    "\n"
)
cat(
    "Study genes outside background:",
    length(
        setdiff(
            study_genes,
            background_genes
        )
    ),
    "\n"
)

#Plots

#Create a volcano plot to see which individual genes differ between tumor and normal

volcano_df <- res_df |> subset (!is.na(padj))#Copy the DESeq2 results so you can modify it without changing the original

volcano_df$significant <- (
    volcano_df$padj <0.05 & abs (
        volcano_df$log2FoldChange
    ) >= 1
) #Creates a true/false column called significant, tests for a significant gene. gene is significant when padj <0.05 & log2FoldChange = 1

volcano_plot <- ggplot(
    volcano_df,
    aes(
        x = log2FoldChange,
        y = -log10(padj)
    )
) +
    geom_point(
        aes(
            alpha = significant
        )
    ) +
    labs (
        title = "GSE103001: Tumor vs Normal Differential Expression",
        x = "log2 Fold Change",
        y = "-log10 Adjusted p-value"
    ) +
    theme_minimal()

#Save volcano plot
ggsave(
    "figures/GSE103001/volcano_plot.png",
    volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
)
#Create a PCA plot
#Does tumor and normal expression profiles differ globally?
 
 vsd <- vst( #variance stabilizing transformation
    dds,
    blind = FALSE
 )
pca_plot <- plotPCA(
    vsd,
    intgroup = "condition"
)

ggsave(
    "figures/GSE103001/pca_plot.png",
    pca_plot,
    width = 8,
    height =6,
    dpi = 300
)
