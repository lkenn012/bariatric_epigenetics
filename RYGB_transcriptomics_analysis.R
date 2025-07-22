### Luke Kennedy (July, 2025)
### Differential gene expression analysis for Pileggi, Mohottalage et al. (2025)
##
## This code imports transcriptomics data quantified via Salmon and conducts
## Differential expression analysis and several visualization of the RNASeq 
## analysis outputs.

# LOAD LIBRARIES
library(tximeta)
library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(ComplexHeatmap)
library(tidyverse)

#####
#  Use tximeta to get salmon RNAseq quantification for DEG analysis
#####
sampleTab_path <- "salmonRNAseq_samples.csv"  # REPLACE WITH CORRECT PATH
colData <- read.csv(sampleTab_path, stringsAsFactors = TRUE)

# load quantification files with tximeta according to colData
se <- tximeta(colData)  # NOTE: this may take some time

# Convert transcript-level data to gene-level
gse <- summarizeToGene(se)
gse$condition <- factor(gse$condition)  # Convert to factor if not already
gse$patient_IDs <- factor(gse$patient_IDs)

gse$condition <- relevel(gse$condition, "Baseline")   # Set correct reference

gse <- addIds(gse, "SYMBOL", gene=TRUE) # add gene symbols

# Construct DESeq dataset object
dds <- DESeqDataSet(gse, design = ~ patient_IDs + condition)
print(dds)

# Only want genes with expression in >25% of samples (n=18)
keep <- rowSums(counts(dds)) >= 18
dds <- dds[keep,]
print(dds)

# Perform DEG analysis
dds <- DESeq(dds)
deg_res <- results(dds, contrast = c("condition", "Post Sx", "Baseline"))

deg_df <- as.data.frame(deg_res)
deg_df$Symbols <- rowData(dds)$SYMBOL

#####
## Save RNAseq data as csv files
#####
write.csv(deg_df, "DESeq2_results_fromSalmon.csv", row.names = TRUE)

# Also want gene x sample values before and after normalization
raw_counts <- assays(gse)$counts
raw_df <- as.data.frame(raw_counts)
write.csv(raw_df, "DESeq2_RNA_rawCounts.csv", row.names = TRUE)

raw_tpm <- assays(gse)$abundance
raw_df <- as.data.frame(raw_tpm)
write.csv(raw_df, "DESeq2_RNA_rawTPMs.csv", row.names = TRUE)

norm_counts <- estimateSizeFactors(dds)
norm_counts <- counts(norm_counts, normalized=TRUE)
norm_df <- as.data.frame(norm_counts)
norm_df$Symbols <- rowData(dds)$SYMBOL
write.csv(norm_df, "DESeq2_RNA_normCounts_wSymbols.csv", row.names = TRUE)

#####
## END-OF: Save RNAseq data as csv files
#####

#####
## Visualize Post-surgery vs. Baseline changes
#####

# Visualize overall changes via PCA, and specific gene changes as volcano

# Use variance stabilize transformed data that is blind to sample and conditions
vsd <- vst(dds)

# use ggplot to visualize paired PCA
plotPCA(vsd, intgroup=c("condition"), ntop=10000)

pca_data <- plotPCA(vsd, intgroup=c("condition","patient_IDs"), ntop=10000, returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
# Define colors for consistency
group_colours <- c("Baseline" = rgb(255, 203, 106, maxColorValue=255),
                   "Post Sx" = rgb(56, 72, 97, maxColorValue=255))

# Plot data
tiff("PCA_rnaVST_CI95.tiff", units="in", width=6, height=5, res=600)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(aes(fill = condition), shape=21, size = 2.5, alpha = 0.8, 
             color="black", stroke=1) + # Points colored by group
  stat_ellipse(aes(fill = condition), type="norm", geom="polygon", alpha=0.1,
               level=0.95, show.legend=FALSE) +  # Ellipses around groups
  scale_color_manual(values = group_colours) +
  scale_fill_manual(values = group_colours) +
  labs(
    title = "PCA of RNAseq (VST)",
    x = paste("Explained Variance PC1 (",percentVar[1],"%)", sep=""),
    y = paste("Explained Variance PC2 (", percentVar[2],"%)", sep=""),
    color = "Sample Group"
  ) +
  theme_classic() +
  theme(legend.position = "right")
dev.off()

# For plotting, first remove NA values
res_df <- na.omit(deg_df)

# Add a column for significance criteria
res_df <- res_df %>%
  mutate(
  Significance = case_when(
  res_df$padj < 0.05 & res_df$log2FoldChange <= -0.5 ~ "Negative DEG",
  res_df$padj < 0.05 & res_df$log2FoldChange >= 0.5 ~ "Positive DEG",
  TRUE ~ "Non-Significant"
  )
  )

# Specify colours for plotting
group_colours <- c("Non-Significant" = "darkgrey",
                  "Negative DEG" = rgb(255, 203, 106, maxColorValue=255), 
                  "Positive DEG" = rgb(56, 72, 97, maxColorValue=255))

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = group_colours) +
  labs(title = "Differential Gene Expression (Post- vs. Pre-)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  theme_classic()

# "Zoomed-in" plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = group_colours) +
  labs(title = "Differential Gene Expression (Post- vs. Pre-)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "slategrey", linewidth = 0.5) +

  xlim(-2.55, 2.55) +
  ylim(0, 25) +
  theme_classic()
#####
## END-OF: Visualize Post-surgery vs. Baseline changes
#####

#####
## Heatmap of select gene expression sets pre- and post surgery
#####
allGenes_of_interest <- read.csv("ribo_geneGroups.csv", header=TRUE) # REPLACE WITH CORRECT PATH

# We are interested in plotting a few different gene sets, specify here
genes_of_interest <- allGenes_of_interest[grep(
  "Transcription|POL I|POL III|45S processing/rRNA modification", 
  allGenes_of_interest$Group), ]
res_df_tidy <- res_df %>%
  rownames_to_column(var = "Ensembl_ID")

mapped_ids_df <- genes_of_interest %>%
  inner_join(res_df_tidy, by = "Symbols") %>%
  distinct(Symbols, Ensembl_ID, Group)

# Sort rows so ordering is the same across objects
mapped_ids_df <- arrange(mapped_ids_df, Symbols)
genes_of_interest <- arrange(genes_of_interest, Symbols)

# Z-score and subset data for plotting
scaled_normDF <- t(scale(t(norm_df[mapped_ids_df$Ensembl_ID,1:70])))

# Formatting for plotting
condition_colours <- list(Group=c("Baseline" = rgb(255, 203, 106, maxColorValue=255), 
                       "Post Sx" = rgb(56, 72, 97, maxColorValue=255))
)
row.names(scaled_normDF) <- mapped_ids_df$Symbols

col_order = c(seq(from=1, to=70, by=2), seq(from=2, to=70, by=2))
sorted_df = scaled_normDF[, col_order]
reorder_cond <- as.character(gse$condition)[col_order]

tiff("Zscore_heat_biogenesisSelect.tiff", units="in", width=4, height=6, res=600)
Heatmap(
  matrix = sorted_df,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  column_split = reorder_cond,
  row_split = mapped_ids_df$Group,  # NOTE: For Mitoribosome plots, comment  
  row_gap = unit(1, "mm"),         #  out row_split,row_gap
  border = TRUE,
  row_title_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=4),
  column_title_gp = gpar(fontsize=8),
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(Group=reorder_cond, name= "Condition", 
                                     col=condition_colours, show_legend=FALSE, 
                                     show_annotation_name=FALSE, 
                                     simple_anno_size=unit(2, "mm")),
  name = "Z-score"
)
dev.off()

# Create Fold change heatmap for each individual to better visualize differences
temp_df <- norm_df[mapped_ids_df$Ensembl_ID,1:70]
logFC_df <- log2(temp_df[, col_order[36:70]] / temp_df[, col_order[1:35]])
row.names(logFC_df) <- mapped_ids_df$Symbols
logFC_df <- logFC_df[is.finite(rowSums(logFC_df)),]  # removes inf and nan rows
logFC <- as.matrix(logFC_df)

# Since genes may have been removed, get genes for plotting
logFC_idsDF <- mapped_ids_df[mapped_ids_df$Symbols %in% row.names(logFC_df), ]


# Plot heatmap
tiff("logFC_heat_biogenesisSelect.tiff", units="in", width=6, height=6, res=600)
Heatmap(
  matrix = logFC,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  row_split = logFC_idsDF$Group,  # NOTE: For Mitoribosome plots, comment
  row_gap = unit(1, "mm"),         #  out row_split,row_gap
  border = TRUE,
  row_title_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=5),
  column_title_gp = gpar(fontsize=8),
  show_column_names = FALSE,
  name = "log2FC Post-Sx/Baseline",
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 4))
)

dev.off()
#####
## END-OF: Heatmap of select gene expression sets pre- and post surgery
#####
