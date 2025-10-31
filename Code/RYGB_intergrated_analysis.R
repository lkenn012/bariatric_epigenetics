### Luke Kennedy (Oct, 2025)
### Integrated analysis of transcriptomics & methylation data for Pileggi, 
### Mohottalage et al. (2025)
##
## This code imports processed data produced from `RYGB_methylation_analysis.R` 
## & `RYGB_transcriptomics_analysis.R` for integration, analysis, and 
## visualizations of genes which potentially undergo epigenetic reprogramming

# LOAD LIBRARIES
library(tidyverse)
library(clusterProfiler)
library(ggplot2)

# We first need to load the relevant data for our integrated analysis
# NOTE: This assumes that both `RYGB_methylation_analysis.R` & 
# `RYGB_transcriptomics_analysis.R` have already been run

# load our data here for convenience
dmr_path <- "C:/Users/User/R/DMRs_Annots_noY_22-10-2025.csv"  # REPLACE WITH CORRECT PATH
dmr_df <- read.csv(dmr_path, header = TRUE, stringsAsFactors = TRUE)

dmp_path <- "DMPs_Annots_noY_19-10-2025.csv"  # REPLACE WITH CORRECT PATH
dmp_df <- read.csv(dmp_path, header = TRUE, stringsAsFactors = TRUE)

normMeth_path <- "Minfi_normM_noY_Luke_29-10-2025.csv"  # REPLACE WITH CORRECT PATH
normMeth_df <- read.csv(normMeth_path, header = TRUE, stringsAsFactors = TRUE)

deg_path <- "C:/Users/User/R/DESeq2_results_fromSalmon_all.csv"  # REPLACE WITH CORRECT PATH
degs_df <- read.csv(deg_path, header = TRUE, stringsAsFactors = TRUE)

norm_path <- "DESeq2_RNA_normCounts.csv"  # REPLACE WITH CORRECT PATH
normRNA_df <- read.csv(norm_path, header = TRUE, stringsAsFactors = TRUE)

# Add unique row identifiers to each dataframe
dmr_df <- dmr_df %>% mutate(dmr_idx = row_number())
degs_df <- degs_df %>% mutate(deg_idx = row_number())

#####
##   MEDEG ANALYSIS
#####

# Finally, identify the genes that are DMRs and DEGs in our methylation and gene
#  expression data, this identifies differentially expressed genes with 
# methylation regulation - meDEGs

# identify DEGs
# Selects for any gene that passes padj < 0.05
sig_DEGs <- degs_df %>%
  mutate(
    # Define groups of genes
    group = case_when(
      padj < 0.05 ~ "sig",
      TRUE ~ "not_sig"
    )
  ) %>%
  filter(!grepl("not_sig", group))

sig_DEGs$orig_deg_idx <- row_number(sig_DEGs)


# Repeat for DMRs
sig_DMRs <- dmr_df %>%
  mutate(
    # Define groups of genes
    group = case_when(
      min_smoothed_fdr < 0.05 & no.cpgs > 1 ~ "sig",
      TRUE ~ "not_sig"
    )
  ) %>%
  filter(!grepl("not_sig", group))

sig_DMRs$orig_deg_idx <- row_number(sig_DMRs)

# Get all gene IDs associated with DMRs
DMR_symbols <- sig_DMRs %>%
  dplyr::select("overlapping.genes") %>%
  separate_rows("overlapping.genes", sep = ",\\s*") %>%
  dplyr::pull("overlapping.genes") %>%
  unique()

# Also get unique gene IDs from DEGs so that intersection
DEG_symbols <- sig_DEGs %>%
  dplyr::select("Symbols") %>%
  separate_rows("Symbols", sep = ",\\s*") %>%
  dplyr::pull("Symbols") %>%
  unique()

intersection_genes <- intersect(DMR_symbols, DEG_symbols)

meth_set <- setdiff(DMR_symbols, intersection_genes)
expr_set <- setdiff(DEG_symbols, intersection_genes)

gene_sets <- c(
  "DMR" = length(meth_set),
  "DEG" = length(expr_set),
  "DMR&DEG" = length(intersection_genes)
)
geneInfo_sets <- c(
  "DMR" = meth_set,
  "DEG" = expr_set,
  "DMR&DEG" = intersection_genes
)
genes_df <- data.frame(name=names(gene_sets), value=gene_sets)
paste0("Numbers of DEGs, DMRs, and meDEGs:")
paste0(gene_sets)
geneInfo_df <- data.frame(name=names(geneInfo_sets), value=geneInfo_sets)
write.csv(geneInfo_df, "DMR_DEG_loose.csv", row.names = FALSE)

## Enrichment analysis of meDEGs

# GO BP enrichment
ego <- enrichGO(gene= intersection_genes,
                OrgDb= org.Hs.eg.db,
                keyType = "SYMBOL",
                ont= "BP",
                universe=as.character(degs_df$Symbols)
)

cluster_summary <- as.data.frame(ego)
write.csv(cluster_summary, file = "meDEG_GObp_enrich.csv", row.names=FALSE)

tiff("meDEG_GOBPenrich.tiff", units="in", width=7, height=5, res=600)
dotplot(ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of meDEGs (FDR<0.05)")
dev.off()

# KEGG enrichment
intersect_keggList <- bitr(intersection_genes, fromType="SYMBOL", toType="ENTREZID", 
                           OrgDb=org.Hs.eg.db)

eKEGG <- enrichKEGG(gene= intersect_keggList$ENTREZID,
                    organism= "hsa",
                    keyType= "ncbi-geneid",
                    universe=all_keggIDs$ENTREZID
)

cluster_summary <- as.data.frame(ego)
write.csv(cluster_summary, file = "meDEG_KEGG_enrich.csv", row.names=FALSE)

tiff("meDEG_KEGG_enrich.tiff", units="in", width=7, height=5, res=600)
dotplot(eKEGG, showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("KEGG Enrichment of meDEGs (FDR<0.05)")
dev.off()

# Plot the genes which has both methylated regions & expression 
# data. Here, plot all DMR values with DEG value on scatter (one point per DMR)

# For each dataset, match on individual gene symbols to identify overlapping genes
dmr_split <- dmr_df %>%
  drop_na("overlapping.genes") %>%
  mutate(
    gene_list = str_split(overlapping.genes, pattern = ",\\s*")
  ) %>%
  unnest(gene_list) %>%
  rename(gene = gene_list) 

degs_split <- degs_df %>%
  drop_na("Symbols") %>%
  mutate(
    gene_list = str_split(Symbols, pattern = ",\\s*")
  ) %>%
  unnest(gene_list) %>%
  rename(gene = gene_list) 

matching_pairs <- inner_join(dmr_split, degs_split, by = "gene", 
                             relationship = "many-to-many")

# Selects for any gene that passes padj < 0.05
grouped_pairData <- matching_pairs %>%
  mutate(
    # Define groups of genes
    group = case_when(
      padj < 0.05 & min_smoothed_fdr < 0.05 ~ "meDEG",
      TRUE ~ "non-sig"
    ),
    subgroup = case_when(
      group == "meDEG" & log2FoldChange > 0 & maxdiff > 0 ~ "hyper-up",
      group == "meDEG" & log2FoldChange > 0 & maxdiff < 0 ~ "hypo-up",
      group == "meDEG" & log2FoldChange < 0 & maxdiff > 0 ~ "hyper-down",
      group == "meDEG" & log2FoldChange < 0 & maxdiff < 0 ~ "hypo-down"
    ),
    most_diff = case_when(
      log2FoldChange > 0.5 & (maxdiff > 0.03 | maxdiff < -0.03) ~ "label",
      log2FoldChange < -0.5 & (maxdiff > 0.03 | maxdiff < -0.03) ~ "label",
      TRUE ~ "no_label"
    )
  )

meDegs <- grouped_pairData[grepl("meDEG", grouped_pairData$group), ]
write.csv(meDegs, 'meDEG_groups.csv', row.names=TRUE)

# plot these data
tiff("diff_expr_meth_plot.tiff", units="in", width=6, height=4, res=600)
ggplot(mapping = aes(x = maxdiff, y = log2FoldChange)) +
  geom_hline(yintercept = 0, color = "gray28", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "gray28", linewidth = 0.25) +
  geom_point(data = grouped_pairData, color = "grey60", size=0.7, alpha=0.2) +
  geom_point(data = meDegs, 
             aes(x = maxdiff, y = log2FoldChange, color = subgroup), 
             size=0.7, alpha=0.8) +
  geom_label_repel(data = subset(meDegs, most_diff == "label"),
                   seed = 42,
                   aes(x = maxdiff, y = log2FoldChange),
                   label = subset(meDegs, most_diff == "label")$Symbols,
                   max.overlaps = 10,
                   min.segment.length = 0.4,
                   nudge_y = 0.01,
                   size = 2,
                   segment.size = 0.15,
                   label.size = 0.1,
                   fill = "white",
                   alpha = 0.8) +
  labs(title = "Differential expression & methylation (FDR <0.05)") +
  theme_classic()
dev.off()

# zoomed plot
tiff("diff_expr_meth_plot_zoom.tiff", units="in", width=6, height=4, res=600)
ggplot(mapping = aes(x = maxdiff, y = log2FoldChange)) +
  geom_hline(yintercept = 0, color = "gray28", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "gray28", linewidth = 0.25) +
  
  geom_point(data = grouped_pairData, color = "grey60", size=0.7, alpha=0.2) +
  geom_point(data = meDegs, 
             aes(x = maxdiff, y = log2FoldChange, color = subgroup), 
             size=0.7, alpha=0.8) +
  geom_label_repel(data = subset(meDegs, most_diff == "label"),
                   seed = 42,
                   aes(x = maxdiff, y = log2FoldChange),
                   label = subset(meDegs, most_diff == "label")$Symbols,
                   max.overlaps = 10,
                   min.segment.length = 0.4,
                   nudge_y = 0.01,
                   size = 2,
                   segment.size = 0.15,
                   label.size = 0.1,
                   fill = "white",
                   alpha = 0.8) +  
  xlim(-0.1, 0.1) +
  ylim(-2, 2) +
  labs(title = "Differential expression & methylation (FDR <0.05)") +
  theme_classic()
dev.off()


#####
##   CORRELATION ANALYSIS
#####

# For each differentially expressed gene, identify associated CpGs
# Correlate normalized gene counts and M-values via Spearman


# Get DEGs and pair with gene expression and methylation

# First get gene associations from methylation site data
# Get relevant methylation data as dataframe
meth_cols <- c("Name", "UCSC_RefGene_Group", 
               "UCSC_RefGene_Name", "GencodeV41_Group",
               "GencodeV41_Name")

subset_meth <- dmp_df[, meth_cols]

cpg_UCSC_long <- subset_meth %>%
  mutate(orig_cpg_idx = row_number()) %>% # original row index for matching
  filter(!is.na(UCSC_RefGene_Name) & str_trim(UCSC_RefGene_Name) != "") %>%
  separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  rename(cpg_UCSC_RefGene = UCSC_RefGene_Name) %>%
  mutate(matched_on_col = "UCSC_RefGene_Name")

cpg_Gencode_long <- subset_meth %>%
  mutate(orig_cpg_idx = row_number()) %>%
  filter(!is.na(GencodeV41_Name) & str_trim(GencodeV41_Name) != "") %>%
  separate_rows(GencodeV41_Name, sep = ";") %>%
  rename(cpg_Gencode_Name = GencodeV41_Name) %>%
  mutate(matched_on_col = "GencodeV41_Name")

# Now get all relevant DEG information
degs_df$ENSG <- sub("\\.\\d+$", "", degs_df$X) # remove version ids

# Get all ids
deg_ids <- sig_DEGs %>%
  pivot_longer(
    cols = c("X", "Symbols", "ENSG"), # "X" is placeholder given to row name (ENSG IDs)
    names_to = "DEG_ID_type",
    values_to = "DEG_IDs"
  ) %>%
  distinct(orig_deg_idx, DEG_IDs, DEG_ID_type) # Get unique IDs

# Now get all matched IDs found in methylation data
# First on UCSC identifier column
matched_genes_UCSC <- cpg_UCSC_long %>%
  dplyr::inner_join(deg_ids, by = c("cpg_UCSC_RefGene" = "DEG_IDs")) %>%
  dplyr::select(orig_deg_idx, orig_cpg_idx, matched_on_col, cpg_UCSC_RefGene, 
                DEG_ID_type, everything())

# Repeat for GENCODE IDs
matched_genes_GENCODE <- cpg_Gencode_long %>%
  dplyr::inner_join(deg_ids, by = c("cpg_Gencode_Name" = "DEG_IDs")) %>%
  dplyr::select(orig_deg_idx, orig_cpg_idx, matched_on_col, cpg_Gencode_Name, 
                DEG_ID_type, everything())

# Combine UCSC and GENCODE matches, and remove any duplicate
all_matches_detailed <- bind_rows(matched_genes_UCSC, matched_genes_GENCODE) %>%
  distinct(orig_deg_idx, orig_cpg_idx, .keep_all = TRUE)

all_matches_detailed$Matched_id <- coalesce(all_matches_detailed$cpg_UCSC_RefGene, 
                                            all_matches_detailed$cpg_Gencode_Name)

# get gene expression IDs that were matched (ENSG with version #) for getting
# RNAseq expression data
all_matches_detailed$matched_ENSG <- sig_DEGs[all_matches_detailed$orig_deg_idx,"X"]
all_matches_detailed$matched_Symbol <- sig_DEGs[all_matches_detailed$orig_deg_idx,"Symbols"]

# Combine all data into a single dataframe
matched_dataDF <- all_matches_detailed[, c("Name", "matched_ENSG", "matched_Symbol")]

# Include columns with consistent headers across data for joining
matched_dataDF <- matched_dataDF %>%
  rename(
    X = matched_ENSG,
    Y = Name
  )

normMeth_df$Y <- rownames(normMeth_df) # create column of site names

# Get data corresponding to each matched ID
matched_dataDF <- matched_dataDF %>%
  left_join(normRNA_df, by = "X") %>%
  left_join(d3, by = "Y")

# Get column names for each dataset for correlations
rna_samples <- colnames(normRNA_df)[-1] # skips column with gene IDs
meth_samples <- colnames(normMeth_df)[-length(colnames(normMeth_df))]

# Compute correlations and p-values across data
corr_df <- matched_dataDF %>%
  rowwise() %>%
  mutate(
    # get columns for each set of data
    rna_vec = list(c_across(all_of(rna_samples))),
    meth_vec = list(c_across(all_of(meth_samples))),
    pearson_corr = cor(unlist(rna_vec), unlist(meth_vec), method = "pearson"),
    pearson_p = cor.test(unlist(rna_vec), unlist(meth_vec), method = "pearson",
                         alternative = "two.sided")$p.value,
    spearman_corr = cor(unlist(rna_vec), unlist(meth_vec), method = "spearman"),
    spearman_p = cor.test(unlist(rna_vec), unlist(meth_vec), method = "spearman",
                          alternative = "two.sided")$p.value
  ) %>%
  ungroup() %>% # Remove rowwise grouping
  dplyr::select(X, Y, matched_Symbol, pearson_corr, pearson_p, spearman_corr, 
                spearman_p, dist_corr, dist_p, rna_samples, meth_samples)

write.csv(corr_df, "paired_sigDEGs_wCpG_corr.csv", row.names=TRUE)

# To visualize the results of correlation analysis, create a manhattan plot

# Find the minimum non-zero p-value to replace 0 with
min_nonzero_p <- min(corr_df$spearman_p[corr_df$spearman_p > 0], na.rm = TRUE)

# Replace p-values of 0 with min_nonzero_p / 10
corr_df$spearman_p[corr_df$spearman_p == 0] <- min_nonzero_p / 10

# Define some genes of interest which have multiple sig corrs for labeling
corrs_of_interest <- c("MYOC1", "MYOM2", "ATP2A1", "SLC16A3", "NME4", "ZNF205",
                       "SELENOW","TPM1", "PITX1", "S100A1", "CDKN1C", "ALDH5A1",
                       "HOMER3", "KANSL2", "FEZ2", "ZNF385A")

# subset and format correlations
corPlt_df <- corr_df %>%
  group_by(matched_Symbol) %>%
  arrange(position, .by_group=TRUE) %>%
  ungroup() %>%
  mutate(gene_order = row_number()) %>%
  mutate(log10_p = -log10(spearman_p)) %>%
  mutate(is_highlight = ifelse(log10_p > 5, "yes", "no")) %>%
  mutate(is_annotate = ifelse(matched_Symbol %in% corrs_of_interest, "yes", "no")) %>%
  mutate(abs_spear = abs(spearman_corr)) %>%   # For plotting corr values
  mutate(is_highlightCorr = case_when(spearman_corr <= -0.5 ~ "Negative",
                                      spearman_corr >= 0.5 ~ "Positive")
  )

# Further format to visualize and annotate relvant genes
significant_genes_summary_df <- corPlt_df %>%
  filter(is_highlight == "yes") %>%
  group_by(matched_Symbol) %>%
  summarise(
    max_log10_p = max(log10_p, na.rm = TRUE),
    max_corr = max(abs_spear, na.rm = TRUE),
    label_gene_order = gene_order[which.max(log10_p)]
  ) %>%
  ungroup() %>%
  mutate(is_annotate_for_plot = ifelse(matched_Symbol %in% corrs_of_interest, "yes", "no"))

axisdf <- corPlt_df %>%
  group_by(matched_Symbol) %>%
  summarize(center=(max(gene_order) + min(gene_order)) / 2)

# get unique genes for manual color scale
plt_genes <- as.factor(corPlt_df$matched_Symbol)

# Plot correlations
tiff("manh_corr_DEG_CpG_allSig.tiff", units="in", width=8, height=6, res=600)
ggplot(corPlt_df, aes(x=gene_order, y=abs_spear)) + # Use pre-calculated log10_p
  
  # Show all points (non-highlighted)
  geom_point(aes(color=as.factor(matched_Symbol)), alpha=0.8, size=0.2) +
  scale_color_manual(values = rep(c("azure3", "azure4"), nlevels(plt_genes))) +
  
  geom_segment(data = subset(significant_genes_summary_df, is_annotate_for_plot == "yes"),
               aes(x = label_gene_order, xend = label_gene_order,
                   y = 0, yend = max_corr),
               color = "darkgreen",
               linewidth = 1,
               alpha = 0.2) +
  
  geom_point(data=subset(corPlt_df, is_highlightCorr=="Positive"), color="goldenrod", size=0.7) + 
  geom_point(data=subset(corPlt_df, is_highlightCorr=="Negative"), color="darkmagenta", size=0.7) + 
  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue", linewidth = 0.75) +
  
  geom_label_repel(data = subset(significant_genes_summary_df, is_annotate_for_plot == "yes"),
                   aes(x = label_gene_order, y = max_corr), # y position is corr
                   label = unique(subset(significant_genes_summary_df, is_annotate_for_plot == "yes")$matched_Symbol), # Use unique labels
                   nudge_y = 0.025,
                   size = 3,
                   label.padding = unit(0.15, "lines"),
                   label.size = 0.1,
                   color = "black",
                   fill = "white",
                   alpha = 0.8) +
  
  # Plot and axis formatting
  scale_x_continuous(expand = c(0.005, 10)) +
  scale_y_continuous(expand = c(-0.1, 0.0)) +
  theme_classic() +
  theme(legend.position="none") +
  ylim(0,1) +
  labs(x = "CpG site (relative position)", y = "Spearman rho (absolute)") # Add meaningful axis labels
dev.off()
