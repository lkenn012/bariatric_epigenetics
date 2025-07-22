### Luke Kennedy (July, 2025)
### Differential gene expression analysis for Pileggi, Mohottalage et al. (2025)
##
## This code imports transcriptomics data quantified via Salmon and conducts
## Differential expression analysis and several visualization of the RNASeq 
## analysis outputs.## Differential DNA methylation analysis for Mohottalage, Pileggi et al. (2025)


# LOAD LIBRARIES
library(tidyverse)
library(minfi)
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
library("org.Hs.eg.db")
library(limma)
library(DMRcate)
library(ggplot2)
library(ggrepel)
library(missMethyl)

### The following code is credited to Majid Nikpay

# Read in methylation sample sheet file
basedir="C:/methylation_data"    # REPLACE WITH CORRECT PATH
targets <- read.metharray.sheet(basedir)

RGset <- read.metharray.exp(targets = targets)
qcReport(RGset, pdf= "qcReport_minfi.pdf")
GRGset <- mapToGenome(RGset)
GRGset <- dropLociWithSnps(GRGset, snps=c("SBE","CpG"), maf=0) #whether to keep SNPs
d1=preprocessQuantile(GRGset)

# Remove Chromosome Y
probe_locations_d1 <- granges(d1)
chromosome_names_d1 <- seqnames(probe_locations_d1)
probes_to_keep <- !(chromosome_names_d1 %in% c("chrY"))
d1 <- d1[probes_to_keep, ]

# Get Beta and M-values
d2 = getBeta(d1)
d3 = getM(d1)
write.csv(d2, "Minfi_normB_noY.csv", row.names = TRUE)
write.csv(d3, "Minfi_normM_noY.csv", row.names = TRUE)

### The following code is credited to Luke Kennedy

# get annotation names
names <- getAnnotation(d1)
write.csv(names, "Minfi_CpGAnnots.csv", row.names = TRUE)

# Consider only gene-associated sites
ucsc_annot <- !is.na(names$UCSC_RefGene_Name) & (names$UCSC_RefGene_Name != "")
gencode_annot <- !is.na(names$GencodeV41_Name) & (names$GencodeV41_Name != "")
probes_to_keep_both_annotated <- ucsc_annot | gencode_annot

d1_genes <- d1[probes_to_keep_both_annotated, ]

# Get Beta and M-values
d2_genes = getBeta(d1_genes)
d3_genes = getM(d1_genes)
write.csv(d2_genes, "Minfi_normB_noY_genesOnly.csv", row.names = TRUE)
write.csv(d3_genes, "Minfi_normM_noY_genesOnly.csv", row.names = TRUE)

# 1. Map genomic regions to corresponding genes
gset <- makeGenomicRatioSetFromMatrix(
  mat = d3_genes,
  array= "IlluminaHumanMethylationEPICv2",
  annotation = "20a1.hg38",
  what = "M"
)

#####
## 2. Determine differentially methylated probes (DMPs) - minfi
#####

# Get Sample IDs and conditions as columns for paired differential analysis
targets$Sample_Group <- ifelse(grepl("PSx", targets$Sample_Name), "Post Sx", "Baseline")
id_match <- str_extract(targets$Sample_Name, "BAR-\\d+")
targets$Sample_ID <- id_match
head(targets)

# Define strings as factors and construct design matrix, including sample pairs
patient_IDs <- factor(targets$Sample_ID)
condition <- factor(targets$Sample_Group)

design <- model.matrix(~ patient_IDs + condition)

m_values_matrix <- getM(gset)

# Fit data with limma for DEG analysis
fit <- lmFit(m_values_matrix, design)
eBayesfit <- eBayes(fit)
summary(decideTests(eBayesfit))

# Get formatted results
temp_table <- topTable(eBayesfit, 
                       coef="conditionPost Sx", 
                       adjust.method="BH", 
                       number=1000000
                       )

names$Probe <- rownames(names)
temp_table$Probe <- rownames(temp_table)
results_annotated <- merge(
  temp_table,
  names,
  by = "Probe",
  all.x = TRUE # Keep all rows from dmp_results
)

rownames(results_annotated) <- results_annotated$Probe
results_annotated$Probe <- NULL # Removed duplicated Probes
results_annotated <- results_annotated[order(results_annotated$adj.P.Val), ]

write.csv(results_annotated, "DMPs_Annots_genesOnly_noY.csv", row.names = TRUE)

#####
## 3. Determine differentially methylated regions (DMRs)
#####
myannotation <- cpg.annotate("array", object=d3_genes, what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type="differential",
                             design=design, coef="conditionPost Sx", fdr=1)

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg38")

write.csv(results.ranges, "DMRs_Annots_noY_noGenes.csv", row.names = TRUE)

# Plot DMRs as a simple volcano plot
DMRplot_data <- as.data.frame(results.ranges) %>%
  mutate(
    neg_log10_FDR = -log10(min_smoothed_fdr),
    # Define significant DMRs
    Significance = case_when(
      min_smoothed_fdr < 0.05 & maxdiff > 0.025 ~ "Significant Up",
      min_smoothed_fdr < 0.05 & maxdiff < -0.025 ~ "Significant Down",
      TRUE ~ "Not Significant"
    )
  )

# Specify group colours for plotting
group_colours <- c("Non-Significant" = "darkgrey",
                   "Significant Down" = rgb(255, 203, 106, maxColorValue=255), 
                   "Significant Up" = rgb(56, 72, 97, maxColorValue=255))
# Volcano plot
tiff("DMR_volcano_all.tiff", units="in", width=8, height=6, res=600)
ggplot(DMRplot_data, aes(x = maxdiff, y = -log10(min_smoothed_fdr), color = Significance)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = group_colours) +
  labs(title = "Differential Methylation regions (Post- vs. Pre-)",
       x = "Max CpG site \u0394 M-value",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = 0.025, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = -0.025, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  xlim(-0.15, 0.15) +
  theme_classic()
dev.off()

# "Zoomed-in" plot
tiff("DMR_volcano_zoom.tiff", units="in", width=8, height=6, res=600)
ggplot(DMRplot_data, aes(x = maxdiff, y = -log10(min_smoothed_fdr), color = Significance)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = group_colours) +
  labs(title = "Differential Methylation regions (Post- vs. Pre-)",
       x = "Maximum \u0394 M-value",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = 0.025, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  geom_vline(xintercept = -0.025, linetype = "dashed", color = "slategrey", linewidth = 0.5) +
  xlim(-0.11, 0.11) +
  ylim(0, 17) +
  theme_classic()
dev.off()

#####
## 4. use missMethyl to get enriched gene sets in DMRs
#####

# Define significant DMRs
sig_DMRs <- results.ranges[results.ranges$no.cpgs > 1 & results.ranges$min_smoothed_fdr < 0.05]
#                            & (results.ranges$maxdiff > 0.025 | results.ranges$maxdiff < -0.025)], strict critera

gst.region <- goregion(sig_DMRs, all.cpg=rownames(d3), 
                       collection="GO", array.type="EPIC_V2", plot.bias=TRUE,
                       genomic.features="ALL", sig.genes="TRUE")
enrich_sets <- topGSA(gst.region, n=250)
write.csv(enrich_sets, "EnrichGO_misMethyl.csv", row.names=TRUE)

# Generate a simple dotplot to visualize some of these results
enrich_sets$logFDR <- -1*log10(enrich_sets$FDR)

tiff("misMethyl_GOenrich_noY.tiff", units="in", width=8, height=4, res=600)
ggplot(data=enrich_sets[1:20,], aes(x=logFDR, y=reorder(TERM,logFDR), size=DE)) +
  geom_point(alpha=0.5, shape=21, color="black", fill="darkviolet") +
  scale_size(range = c(1, 8), name="Gene count") +
  theme(legend.position="bottom") +
  ylab("GO term") +
  xlab("-log10 FDR") +
  ggtitle("GO Enrichment DMRs (FDR < 0.05)") +
  theme_classic()

dev.off()

# Repeat for KEGG
gst.region <- goregion(sig_DMRs, all.cpg=rownames(d3), 
                       collection="KEGG", array.type="EPIC_V2", plot.bias=TRUE,
                       genomic.features="ALL", sig.genes="TRUE")
enrichKEGG_set <- topGSA(gst.region, n=250)

# Generate a simple dotplot to visualize some of these results
enrichKEGG_set$logFDR <- -1*log10(enrichKEGG_set$FDR)

tiff("misMethyl_KEGGenrich_noY.tiff", units="in", width=8, height=4, res=400)
ggplot(data=enrich_sets[1:20,], aes(x=logFDR, y=reorder(Description,logFDR), size=DE)) +
  geom_point(alpha=0.5, shape=21, color="black", fill="deeppink3") +
  scale_size(range = c(1, 8), name="Gene count") +
  theme(legend.position="bottom") +
  ylab("GO term") +
  xlab("-log10 FDR") +
  ggtitle("GO Enrichment DMRs (FDR < 0.05)") +
  theme_classic()

dev.off()

# Visualize distribution of methylated regions for genes through a simple barplot
sig_DMRs <- sig_DMRs %>%
  direction = ifelse(maxdiff > 0,"pos", "neg")

up_positions <- sig_DMRs[sig_DMRs[,direction == "pos"], c("start", "end")]

# Identify CpG sites constituting DMRs

# ~~~~~ NOTE: THIS STEP IS VERY SLOW ~~~~~
locs <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations
locs.ranges <- GRanges(locs$chr, IRanges(locs$pos, locs$pos))
names(locs.ranges) <- rownames(locs)
constituent.cpgs <- list()
for (i in 1:length(results.ranges)){
  constituent.cpgs <- c(constituent.cpgs, locs.ranges[locs.ranges %over% results.ranges[i]])
}

# Get all CpG site positions for each DMR as a new column in our DMR output
list_of_start_positions <- lapply(constituent.cpgs, start)
cpg_positions_aggregated_string <- sapply(list_of_start_positions, function(positions_vector) {
  paste(positions_vector, collapse = ", ")
})

# Get the positions for significant DMRs
sig_idxs <- which(results.ranges$min_smoothed_fdr < 0.05)
sig_positions <- list_of_start_positions[sig_idxs]

# Get directionality for later plotting
sig_results <- results.ranges[sig_idxs,]
up_idxs <- which(sig_results$maxdiff > 0)

# Can save results for future reference
results.ranges$cpg_positions <- cpg_positions_aggregated_string
write.csv(results.ranges, "DMRs_wSites_noY.csv", row.names = TRUE)

# For each DMR site, get the relative location to the gene for plotting
results.ranges$All_groups <- paste(results.ranges$UCSC_RefGene_Group, results.ranges$GencodeV41_Group, sep=";")
value_lookup_table <- setNames(results.ranges$All_groups, results.ranges$pos)

# Function for extracting and aggregating list of locations
get_loc <- function(value_table, loc_vector) {
  found_values <- value_table[as.character(loc_vector)] 
  found_values <- na.omit(found_values)
  
  # Aggregate into a single comma-separated string
  if (length(found_values) > 0) {
    return(paste(found_values, collapse = ";"))
  } else {
    return(NA_character_) # Return NA if no values were found
  }
}

# Function for get unique values from list
extract_unique <- function(s, split_regex) {
  # Split string into list of locations
  split_values <- unlist(strsplit(s, split = split_regex))
  split_values <- split_values[split_values != ""] # remove empty strings
  
  unique_vals <- unique(split_values)
  
  return(sort(unique_vals))
}

# Get unique locations for each DMR
upDMR_positions <- lapply(sig_positions[up_idxs], get_loc, 
                          value_table = value_lookup_table) 

delimiter_regex <- ",\\s*|;\\s*" # for splitting string of locations into list
unique_upDMR_pos <- lapply(upDMR_positions, extract_unique, 
                           split_regex=delimiter_regex) 

# Repeat for down DMRs
downDMR_positions <- lapply(sig_positions[-up_idxs], get_loc, 
                          value_table = value_lookup_table) 
unique_downDMR_pos <- lapply(downDMR_positions, extract_unique, 
                           split_regex=delimiter_regex) 

# Define function for getting counts for the locations of interest
get_annotation_counts<- function(annotations_list, group_name) {
  all_annotations <- unlist(annotations_list)
  # Distinguish 1st Exon and remainder of Gene Body
  all_annotations <- gsub("\\bexon_1\\b", "1st Exon", all_annotations)
  all_annotations <- gsub("exon_[0-9]+", "Gene Body", all_annotations)
  counts_df <- as.data.frame(table(all_annotations))
  names(counts_df) <- c("Category", "Count")
  counts_df$Change_Type <- group_name # Add a column to identify the group
  return(counts_df)
}

# Get counts & combine
pos_counts_df <- get_annotation_counts(unique_upDMR_pos, "Positive DMR")
neg_counts_df <- get_annotation_counts(unique_downDMR_pos, "Negative DMR")

combined_counts_df <- rbind(pos_counts_df, neg_counts_df)

# Reorder the 'Category' factor levels in the combined dataframe
# so that the bars appear in a consistent order (e.g., by total frequency)
combined_counts_df$Category <- factor(combined_counts_df$Category,
                                      levels = total_category_counts$Category)

# Get percentage for plotting as well
total_pos <- sum(combined_counts_df[combined_counts_df$Change_Type == "Positive DMR", "Count"])
total_neg <- sum(combined_counts_df[combined_counts_df$Change_Type == "Negative DMR", "Count"])

combined_counts_df$Percent = case_when(combined_counts_df$Change_Type == "Positive DMR" ~ (combined_counts_df$Count/total_pos)*100,
                                       combined_counts_df$Change_Type == "Negative DMR" ~ (combined_counts_df$Count/total_neg)*100
)

# Plot these occurences
tiff("DMR_regionCounts.tiff", units="in", width=8, height=4, res=600)
ggplot(combined_counts_df, aes(x = Category, y = Count, fill = Change_Type)) +
  geom_bar(stat = "identity", width=0.5, position = position_dodge(width=-0.5),
            color = "white") +
  scale_fill_manual(values = c("Negative DMR" = rgb(255, 203, 106, maxColorValue=255), 
                               "Positive DMR" = rgb(56, 72, 97, maxColorValue=255))) +
  labs(title = "Frequency of DMR Genomic overlap (Gene-associated DMRs)",
       y = "Number of Occurrences",
       fill = "Methylation change (Post/Pre)") + # Legend title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate x-axis labels
    panel.grid.major.x = element_blank(), # Remove vertical grid lines for clarity
    plot.title = element_text(face = "bold", hjust = 0.5) # Center and bold title
  )
dev.off()

# Plot these %
tiff("DMR_regionPercent.tiff", units="in", width=8, height=4, res=600)
ggplot(combined_counts_df, aes(x = Category, y = Percent, fill = Change_Type)) +
  geom_bar(stat = "identity", width=0.5, position = position_dodge(width=-0.5),
           color = "white") +
  scale_fill_manual(values = c("Negative DMR" = rgb(255, 203, 106, maxColorValue=255), 
                               "Positive DMR" = rgb(56, 72, 97, maxColorValue=255))) +
  labs(title = "Frequency of DMR Genomic overlap (Gene-associated DMRs)",
       y = "% Occurences",
       fill = "Methylation change (Post/Pre)") + # Legend title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate x-axis labels
    panel.grid.major.x = element_blank(), # Remove vertical grid lines for clarity
    plot.title = element_text(face = "bold", hjust = 0.5) # Center and bold title
  )
dev.off()
#####
##  Plot methylation values for CpG sites in DMRs of interest
#####

# A list of CpG site positions corresponding to DMRs in `results.ranges`
DMR_cpg_locs <- list("ACACB" = c(109131279, 109131375, 109131325, 109131311,
                              109131254, 109130433),
                  "PTPRE" = c(127998673, 127999496, 127999576, 127999513,
                              127997935, 127998755, 127997965),
                  "PTPRE, AL390236.1" = c(128030456, 128030506, 128030534,
                                          128029175, 128030233, 128030732,
                                          128029588),
                  "PPARGC1A" = c(23890170, 23890175, 23890460, 23890998,
                                 23891266, 23890917, 23890212, 23889973,
                                 23890892),
                  "AC002451.1, PDK4" = c(95592780, 95592792) # AC002451.1, PDK4
)

## Define functions for getting CpG data from positions in DMR_cpg_locs and then
## plotting this data.

# function for processing CpG data for plots
get_cpgData <- function(cpg_pos, methylation_data, cpg_info, sample_info) {
  
  # Get CpG names and corresponding values
  DMR_sites <- cpg_info[cpg_info$pos %in% cpg_pos,]
  DMR_sites <- as.data.frame(DMR_sites)
  DMR_sites <- DMR_sites %>%   # Sort by descending pos since on reverse strand
    dplyr::arrange(desc(pos))
  
  DMRsite_vals <- methylation_data[DMR_sites$Name,]
  print(DMR_sites)
  # Format data for plotting
  DMRsite_vals <- as.data.frame(DMRsite_vals) %>%
    tibble::rownames_to_column("CpG_site") %>%
    mutate(GencodeV41_Group = DMR_sites$GencodeV41_Group, 
           UCSC_RefGene_Group = DMR_sites$UCSC_RefGene_Group,
           cpg_number = row_number()
           )
  print(unique(DMRsite_vals$GencodeV41_Group))
  print(unique(DMRsite_vals$UCSC_RefGene_Group))
  
  DMRsite_vals_long <- DMRsite_vals %>%
    pivot_longer(
      cols = -c(CpG_site, GencodeV41_Group, UCSC_RefGene_Group, cpg_number),
      names_to = "Sample_ID",
      values_to = "Norm_Beta" 
    )

  # Get sample identifiers equivalent to those in methylation data
  meth_sampleData <- sample_info[, c("Sample_Group","Basename")] # Get group and meth sample ID
  meth_sampleData$Sample_ID <- substr(meth_sampleData$Basename, 
                                      nchar(meth_sampleData$Basename)-18, 
                                      nchar(meth_sampleData$Basename)
                                      )
  
  # Add Group and CpG regions
  DMRsite_vals_final <- DMRsite_vals_long %>%
    mutate(Region = case_when(
      grepl("TSS|exon_1", GencodeV41_Group) ~ "Promoter",
      grepl("TSS|exon_1", UCSC_RefGene_Group) ~ "Promoter",
      TRUE ~ "Gene Body"
      )
      ) %>%
    left_join(meth_sampleData, by = "Sample_ID")
  
  # For plotting remove probe IDs from CpG sites
  DMRsite_vals_final$CpG_name <- substr(DMRsite_vals_final$CpG_site, 1,
                                        nchar(DMRsite_vals_final$CpG_site)-5
  )
  
  # Get some summary stats for plotting
  summary_data <- DMRsite_vals_final %>%
    group_by(CpG_name, Sample_Group) %>%
    summarise(
      Mean_meth = mean(Norm_Beta),
      Median_meth = median(Norm_Beta),
      SD_meth = sd(Norm_Beta),
      gene_region = Region,
      rel_cpgPos = cpg_number,
      .groups = 'drop' # Drop grouping after summarizing
    )

  # For labeling plots, get CpG sites for each region
  region_locs <-list()
  region_types <- unique(summary_data$gene_region)
  for (region in region_types) {
    region_sites <- summary_data[summary_data$gene_region == region, ]
    region_locs[[region]] <- c(min(region_sites$rel_cpgPos), max(region_sites$rel_cpgPos))
  }
  
  # Format region locations for plotting
  region_df <- data.frame(region_locs)
  region_df <- as.data.frame(t(region_df))
  colnames(region_df) <- c("min", "max")
  region_df$gene_region <- rownames(region_df)
  
  return(list(summary = as.data.frame(summary_data), 
              region_pos = region_df)
         )
}

# function for plotting
plot_cpgs <- function(gene_name, summary_meth, meth_regions, group_formats) {
  # Plot the values organized by position
  plt <- ggplot(summary_meth, aes(x = CpG_name, y = Mean_meth, 
                                  group = Sample_Group, color = Sample_Group)) +
    geom_line(linewidth=0.75, alpha=0.9) +
    geom_point() +
    geom_pointrange(aes(ymin=Mean_meth-SD_meth, ymax=Mean_meth+SD_meth),
                    size=0.2) +
    scale_color_manual(values = group_formats) +
        # NOTE: this does not properly format - redo manually
    # Add annotation labels for CpG site regions
    geom_rect(data = meth_regions,
              aes(xmin = min, xmax = max, ymin = 0.0, ymax = 0.05),
              fill = "gray", alpha = 0.7, inherit.aes = FALSE) + # inherit.aes=FALSE is important
    # Add region text labels
    geom_text(data = meth_regions,
              aes(x = (max+min)/2, y = 0.042, label = gene_region),
              vjust = 1, # Align top of text with y position
              size = 3, fontface = "italic", inherit.aes = FALSE) +
    scale_color_manual(values = group_colours) +
    labs(
      title = paste0("DMR methylation across CpG sites (",gene_name, ")"),
      x = "CpG Site",
      y = "CpG β-value (normalized)"
    ) +
    ylim(0,1) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 65,     # Rotate labels by 90 degrees
      hjust = 1,      # Horizontal justification: 1 means align right (bottom of rotated label)
      size = 10,      # Adjust font size if needed
      face = "plain"  # Font face (e.g., "bold", "italic", "plain")
    )
    )
  return(plt)
}

# Define colours for plotting
group_colours <- c("Baseline" = rgb(255, 203, 106, maxColorValue=255), 
                   "Post Sx" = rgb(56, 72, 97, maxColorValue=255)
                   )
all_gene_plots <- list()

# Create plots for cpg sites of each gene
for (gene_name in names(DMR_cpg_locs)) {
  current_gene_cpgs <- DMR_cpg_locs[[gene_name]]
  current_gene_summary <- get_cpgData(current_gene_cpgs, d2, names, targets)
  gene_plt <- plot_cpgs(gene_name, current_gene_summary$summary, current_gene_summary$region_pos, group_colours)
  all_gene_plots[[gene_name]] <- gene_plt
  ggsave(filename=paste0(gene_name,"DMR_cpg_line_wSD.tiff"),
         plot=gene_plt,
         width=5,
         height=4)
}

#####
## PCA of M-/B-values between groups
#####

# Define colors for sample groups
group_levels <- levels(condition)
color_palette <- c(rgb(255, 203, 106, maxColorValue=255), rgb(56, 72, 97, maxColorValue=255))

# Create a vector of colors for each sample based on its group
sample_colors <- color_palette[as.numeric(condition)]

mds_data <- plotMDS(d3, top=nrow(d3), gene.selection="common", plot="FALSE")
plotMDS(d3, top=1000, gene.selection="common", plot="TRUE",
        labels = NULL,          # Remove text labels
        pch = 16,               # dots
        col = sample_colors,    # Vector of colors for each point
        cex = 1.5,              # Adjust point size
        main = "MDS Plot of Sample Methylation Profiles (Pre- & Post)"
)

# Can use ggplot to create MDS plot showing paired points
plot_df <- data.frame(
  MDS1 = mds_data$x,
  MDS2 = mds_data$y,
  Sample_Group = condition,
  Patient_ID = patient_IDs
)

# Define colors for consistency
group_colors <- c("Baseline" = rgb(255, 203, 106, maxColorValue=255), "Post Sx" = rgb(56, 72, 97, maxColorValue=255))

# Paired plot
tiff("PCA_methNormM_paired.tiff", units="in", width=6, height=5, res=600)
ggplot(plot_df, aes(x = MDS1, y = MDS2, group = Patient_ID)) +
  geom_line(color = "grey70", linewidth = 1, alpha = 0.9) + # Lines connecting pairs
  geom_point(aes(fill = Sample_Group), colour="black", pch=21, size=4, alpha=1) + # Points colored by group
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colours) +
  labs(
    title = "PCA Pre- & Post-RYGB methylation profiles (norm. M-values)",
    x = paste("Principal Component 1 (", round(mds_data$var.explained[1]*100,1),"%)", sep=""),
    y = paste("Principal Component 2 (", round(mds_data$var.explained[2]*100,1),"%)", sep=""),
    color = "Sample Group"
  ) +
  xlim(-0.22, 0.15) +
  ylim(-0.18, 0.22) +
  theme_classic() +
  theme(legend.position = "right")
dev.off()

#####
##   CORRELATION ANALYSIS
#####

# For each differentially expressed gene, identify associated CpGs
# Correlate normalized gene counts and M-values via Spearman

# Load transcriptomics data
##### NOTE: This assumes "RYGB_transcriptomics_analysis" code has already been
##### run. Outputs from that analysis are necessary for correlation analysis
deg_path <- "DESeq2_results_fromSalmon.csv"  # REPLACE WITH CORRECT PATH
degs_df <- read.csv(deg_path, header=TRUE, stringsAsFactors = TRUE)

norm_path <- "DESeq2_RNA_normCounts.csv"  # REPLACE WITH CORRECT PATH
normRNA_df <- read.csv(norm_path, header=TRUE, stringsAsFactors = TRUE)

# Get DEGs and pair with gene expression and methylation

# First get gene associations from methylation site data
# Get relevant methylation data as dataframe
meth_cols <- c("Name", "UCSC_RefGene_Group", 
               "UCSC_RefGene_Name", "GencodeV41_Group",
                "GencodeV41_Name")

subset_meth <- results_annotated[, meth_cols]
subset_meth <- as.data.frame(subset_meth)

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

# Selects for any gene that passes padj < 0.05
sig_DEGs <- degs_df %>%
  mutate(
    # Define groups of genes
    group = case_when(
      padj < 0.05 ~ "sig",
#      padj < 0.05 & log2FoldChange <= -0.5 ~ "neg_DEG",
#      padj < 0.05 & log2FoldChange > 0.5 ~ "pos_DEG",
      TRUE ~ "not_sig"
    )
  ) %>%
  filter(!grepl("not_sig", group))

sig_DEGs$orig_deg_idx <- row_number(sig_DEGs)

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

d3 <- as.data.frame(d3)
d3$Y <- rownames(d3) # create column of site names

# Get data corresponding to each matched ID
matched_dataDF <- matched_dataDF %>%
  left_join(normRNA_df, by = "X") %>%
  left_join(d3, by = "Y")

# Get column names for each dataset for correlations
rna_samples <- colnames(normRNA_df)[-1] # skips column with gene IDs
meth_samples <- colnames(d3)[-length(colnames(d3))]

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
# Order genes according to the median position of methylation sites
# Points should be visually distinct between genes (e.g., grey - black)
# And points/genes of interest labeled/highlighted
corr_df <- read.csv("sigRNA_CpG_corr_summary_04-06-2025.csv")

# Find the minimum non-zero p-value to replace 0 with
min_nonzero_p <- min(corr_df$spearman_p[corr_df$spearman_p > 0], na.rm = TRUE)

# Replace p-values of 0 with min_nonzero_p / 10
corr_df$spearman_p[corr_df$spearman_p == 0] <- min_nonzero_p / 10

# Define some genes of interest which have multiple sig corrs
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

# Make the plot
tiff("manh_corrP_allSig.tiff", units="in", width=8, height=6, res=600)
ggplot(corPlt_df, aes(x=gene_order, y=log10_p)) +
  
  # Show all points (non-highlighted)
  geom_point(aes(color=as.factor(matched_Symbol)), alpha=0.7, size=0.4) +
  scale_color_manual(values = rep(c("azure3", "azure4"), nlevels(plt_genes))) +
  
  # Highlight region of plots with highly correlated genes
  geom_segment(data = subset(significant_genes_summary_df, is_annotate_for_plot == "yes"),
               aes(x = label_gene_order, xend = label_gene_order,
                   y = 0, yend = max_log10_p),
               color = "darkgoldenrod",
               linewidth = 1,
               alpha = 0.2) +
  
  # Add highlighted points
  geom_point(data=subset(corPlt_df, is_highlight=="yes"), color="orange", size=0.7) + # Slightly larger size
  
  # Significance threshold
  geom_hline(yintercept = 5, linetype = "dashed", color = "blue", linewidth = 0.8) +

  
  # Add label over region for a gene
  geom_label(data = subset(significant_genes_summary_df, is_annotate_for_plot == "yes"),
             aes(x = label_gene_order, y = max_log10_p), # y position is max_log10_p
             label = unique(subset(significant_genes_summary_df, is_annotate_for_plot == "yes")$matched_Symbol), # Use unique labels
             position = position_jitter(1000, 0.3), # Nudge label slightly above the line/point
             size = 3,
             label.padding = unit(0.15, "lines"),
             label.size = 0.1,
             color = "black",
             fill = "white",
             alpha = 0.8) +
  
  scale_x_continuous(expand = c(0.005,0.01)) +
  scale_y_continuous(expand = c(-0.01, 0.25)) +
  theme_classic() +
  theme(legend.position="none") +
  labs(x = "Gene CpG sites", y = "-log10(p-value)") # Add meaningful axis labels
dev.off()

# Repeat plotting using correlation values instead of p-values
tiff("manh_corrRho_allSig_test.tiff", units="in", width=8, height=6, res=600)
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

#####
## Correlation analysis and plot between baseline CpG methylation and PWL
#####

# Correlate normalized beta values pre-Surgery with age-adjusted percent weight loss
age_pwl <- read.csv("age_adjust_pwl.csv", header=TRUE)

# Get baseline CpG methylation
# Data is in the same order as pwl, with Pre-, post- consecutive
pre_Bmeth <- d2[, c(seq(from=1, to=70, by=2))]

# Compute correlations, p-values, and FDR
pwl_corrs <- apply(pre_Bmeth, MARGIN=1, FUN=cor.test, 
                   y=age_pwl$Age.adjusted.PWL, method="pearson", 
                   alternative="two.sided")
# Compute FDRs
pwlCorr_p <- sapply(pwl_corrs, function(x) x$p.value)
pwlCorr_FDR <- p.adjust(pwlCorr_p, method = "BH")

pwlCorr_df <- data.frame(names(pwl_corrs), 
                         names$chr,
                         names$pos,
                         names$UCSC_RefGene_Name,
                         names$GencodeV41_Name,
                         sapply(pwl_corrs, function(x) x$estimate), 
                         pwlCorr_p, 
                         -1*log10(pwlCorr_p),   # -log10*p for plotting
                         pwlCorr_FDR)
colnames(pwlCorr_df) <- c("CpG", "chr", "pos", "UCSC_RefGene","GencodeV41_Gene", 
                          "Pearson_rho", "Pearson_p", "log10_p", "Pearson_FDR")
# Save correlations
write.csv(pwlCorr_df, "ageAdjPWL_CpGnormB_corrs.csv", row.names = TRUE)

# Modify dataframe for plotting
pwlCorr_df$chr <- substr(pwlCorr_df$chr, 4, nchar(pwlCorr_df$chr))
pwlCorr_df$chr <- factor(pwlCorr_df$chr, 
                       c(as.character(1:22), "X", "Y")) # for proper ordering
  
pwlCorPlt_df <- pwlCorr_df %>% 
  mutate(
    is_highlight = ifelse(log10_p > 5, "yes", "no"),
    is_annotate = ifelse(log10_p > 6, "yes", "no"),
    CpG_name = substr(CpG, 1, nchar(CpG)-5),
    chr_numeric = as.numeric(chr),
    chr_jittered = jitter(chr_numeric, amount = 0.4), # For consistent labeling
    absRho = abs(Pearson_rho),
    is_highlightCorr = case_when(Pearson_rho <= -0.6 ~ "Negative",
                                 Pearson_rho >= 0.6 ~ "Positive"),
    is_annotateCorr = ifelse(absRho > 0.7, "yes", "no")
    )

axis_labels <- levels(pwlCorPlt_df$chr)
axis_breaks <- as.numeric(seq_along(axis_labels))


# TEST using relative positions rather than absolute
pwlCorPlt_df$abs_pos <- pwlCorPlt_df$pos
pwlCorPlt_df$pos <- 1:nrow(pwlCorPlt_df)

##### NEW
chr_max_pos <- pwlCorPlt_df %>%
  group_by(chr) %>%
  summarise(max_pos = max(as.numeric(pos), na.rm = TRUE)) %>% # Make sure 'pos' is your genomic position column
  ungroup() %>%
  arrange(chr) # Arrange by the factor levels

# Calculate offsets (cumulative sum of previous chromosome lengths)
data_cum <- chr_max_pos %>%
  mutate(offset = lag(cumsum(max_pos), default = 0)) %>%
  select(chr, offset)

# Join the offset back to the original data and calculate cumulative BP
pwlCorPlt_df_manhattan <- pwlCorPlt_df %>%
  left_join(data_cum, by = "chr") %>%
  mutate(cumulative_bp = pos + offset) # 'pos' is your original position column

# Calculate chromosome midpoints for x-axis labels
axis_df <- pwlCorPlt_df_manhattan %>%
  group_by(chr) %>%
  summarise(
    chr_start = min(pos, na.rm = TRUE),
    chr_end = max(pos, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(chr_mid = (chr_start + chr_end) / 2)
#####

# Make the plot
tiff("manh_corrPWL_allSig.tiff", units="in", width=8, height=6, res=600)
ggplot(pwlCorPlt_df, aes(x=chr_jittered, y=log10_p)) + # Use pre-calculated log10_p
  
  # Show all points (non-highlighted)
  geom_jitter(aes(color=as.factor(chr)), alpha=0.7, size=0.4, 
              position=position_jitter(seed = 42)) +
  scale_color_manual(values = rep(c("grey", "lightgoldenrod3"), 
                                  length(levels(pwlCorr_df$chr)))) +
  geom_label(data = subset(pwlCorPlt_df, is_annotate == "yes"),
             aes(x = chr_jittered, y = log10_p, label=CpG_name),
             position = position_jitter(0.5, 0.3, seed=42), # Nudge labels
             size = 2.5,
             label.padding = unit(0.15, "lines"),
             label.size = 0.1,
             color = "black",
             fill = "white",
             alpha = 0.8) +
  scale_x_continuous(
    breaks = axis_breaks, # Numeric positions for breaks
    labels = axis_labels, # Original chromosome names for labels
    expand = c(0.01, 0) # Reduce space at the ends
  ) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange", linewidth = 0.8) +
  geom_point(data=subset(pwlCorPlt_df, is_highlight=="yes"), 
             aes(x=chr_jittered, y=log10_p), 
             color="blue", size=0.7) +

  theme_classic() +
  theme(legend.position="none") +
  labs(x = "Chromosome", y = "-log10(p-value)")
dev.off()

# Repeat using correlation values
tiff("manh_corrPWL_Rho06_test.tiff", units="in", width=8, height=6, res=600)
ggplot(pwlCorPlt_df_manhattan, aes(x=pos, y=absRho)) +
  
  # Show all points (non-highlighted)
  geom_point(aes(color=chr), alpha=0.7, size=0.4) +
  scale_color_manual(values = rep(c("grey", "lightgoldenrod3"), 
                                  length(levels(pwlCorPlt_df_manhattan$chr)))) +
  
  geom_label_repel(data = subset(pwlCorPlt_df_manhattan, is_annotateCorr == "yes"),
             aes(x = pos, y = absRho, label=CpG_name),
             nudge_y = 0.04,
             size = 3,
             label.padding = unit(0.15, "lines"),
             label.size = 0.1,
             color = "black",
             fill = "white",
             alpha = 0.8) +
  scale_x_continuous(
    breaks = axis_df$chr_mid, # Midpoints for labels
    labels = axis_df$chr,     # Chromosome labels
    expand = c(0.01, 0)     # Reduce empty space at ends
    ) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "orange", linewidth = 0.8) +
  geom_point(data=subset(pwlCorPlt_df_manhattan, is_highlightCorr=="Positive"), 
             aes(x=pos, y=absRho), 
             color="dodgerblue3", size=0.7) +
  geom_point(data=subset(pwlCorPlt_df_manhattan, is_highlightCorr=="Negative"), 
             aes(x=pos, y=absRho), 
             color="darkseagreen3", size=0.7) +
  
  theme_classic() +
  ylim(0,1) +
  theme(axis.text.x=element_text(size=6), 
        legend.position="none") +
  labs(x = "Chromosome", y = "Pearson rho (absolute)")
dev.off()

#####
##   END-OF: CORRELATION ANALYSIS
#####


# For one more visulization, identify the genes that are DMRs and DEGs in our
# methylation and expression data, this identifies differentially expressed
# genes with methylation regulation - mDEGs
library(eulerr)

sig_DMR_df <- as.data.frame(sig_DMRs)
# Get all gene IDs associated with DMRs
DMR_symbols <- sig_DMR_df %>%
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
geneInfo_df <- data.frame(name=names(geneInfo_sets), value=geneInfo_sets)
write.csv(geneInfo_df, "DMR_DEG_loose_NoGene.csv", row.names = FALSE)

# Plot as venn diagram
fit <- euler(gene_sets)
tiff("DMR_DEG_venn_noYonlyGene_CpG.tiff", units="in", width=8, height=6, res=600)
plot(fit,
     quantities = TRUE, # Show the counts in each region
     labels = c("DMR Genes", "DEG Genes"), # Labels for your circles
     fills = list(fill = c("blue", "red"), alpha = 0.2), # Fill colors for areas
     edges = list(col = c("cornflowerblue", "lightcoral"), lwd = 3, alpha = 0.7), # Outline colors and line width
     main = "Overlap of DMR and DEG Genes" # Plot title
)
dev.off()
####
# Enrichment analysis to DEGs and mDEGs
####
library(clusterProfiler)

genes_df <- read.csv("DMR_DEG_loose.csv", header = TRUE)
intersection_genes <- genes_df[grep("DMR&DEG", genes_df$name), "value"]

# ENSG IDs are prefered over Symbols due to duplicates
DEG_IDlist <- sig_DEGs %>%
  dplyr::select("ENSG") %>%
  dplyr::pull("ENSG") %>%
  unique()

ego <- enrichGO(gene= DEG_IDlist,
                OrgDb= org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "ALL",
                universe=as.character(degs_df$ENSG)
                )

cluster_summary <- as.data.frame(ego)
write.csv(cluster_summary, file = "DEG_allGO_enrich.csv", row.names=FALSE)

s_ego <- clusterProfiler::simplify(ego)
dotplot(ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of differentially expressed genes (FDR<0.05)")

# KEGG enrichment
DEG_keggList <- bitr(DEG_IDlist, fromType="ENSEMBL", toType="ENTREZID", 
                     OrgDb=org.Hs.eg.db)
all_keggIDs <- bitr(as.character(degs_df$ENSG), fromType="ENSEMBL", 
                    toType="ENTREZID", OrgDb=org.Hs.eg.db)

eKEGG <- enrichKEGG(gene= DEG_keggList$ENTREZID,
                organism= "hsa",
                keyType= "ncbi-geneid",
                universe=all_keggIDs$ENTREZID
                )

cluster_summary <- as.data.frame(eKEGG)
write.csv(cluster_summary, file = "DEG_KEGG_enrich.csv", row.names=FALSE)
tiff("DEG_KEGGenrich.tiff", units="in", width=7, height=4, res=600)

dotplot(eKEGG,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_classic() +
  ggtitle("KEGG Enrichment of differentially expressed genes (FDR<0.05)")
dev.off()
# Up-down gene set overlap plot
# Identify up- and down-regulated genes
upDown_DEGs <- degs_df %>%
  mutate(
    # Define groups of genes
    group = case_when(
      padj < 0.05 & log2FoldChange < 0 ~ "negative DEG",
      padj < 0.05 & log2FoldChange > 0 ~ "positive DEG",
      TRUE ~ "not_sig"
    )
  ) %>%
  filter(!grepl("not_sig", group))

upDown_list <- list(upDown_DEGs[upDown_DEGs$group == "positive DEG", "ENSG"], 
                    upDown_DEGs[upDown_DEGs$group == "negative DEG", "ENSG"])
names(upDown_list)<-c("Up-regulated","Down-regulated")

cclust<-compareCluster(ENSG~group, data=upDown_DEGs,
                       fun = enrichGO,
                       OrgDb= org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont= "BP",
                       universe=as.character(degs_df$ENSG))

tiff("DEG_splitBPenrich.tiff", units="in", width=7, height=7, res=600)
dotplot(cclust,showCategory=15, label_format=60) +
  theme_minimal()
dev.off()

cluster_summary <- as.data.frame(cclust)
write.csv(cluster_summary, file = "DEG_splitBPenrich.csv", row.names=FALSE)

# Repeat for KEGG
upDown_keggList <- bitr(upDown_DEGs$ENSG, fromType="ENSEMBL", toType="ENTREZID", 
                     OrgDb=org.Hs.eg.db)

# some ensembl have multiple NCBI for the same protein, remove any duplicates
upDown_keggList <- upDown_keggList %>%
  distinct(ENSEMBL, .keep_all = TRUE)

upDown_DEGs$group <- as.factor(upDown_DEGs$group) # Explicitly convert 'group' to factor
upDown_DEGs$ENTREZ <- as.character(upDown_keggList$ENTREZID)

cclust<-compareCluster(ENTREZ~group, data=upDown_KEGG,
                       fun = enrichKEGG,
                       organism = "hsa",
                       universe=all_keggIDs$ENTREZID)
tiff("DEG_splitKEGGenrich.tiff", units="in", width=7, height=7, res=600)
dotplot(cclust,showCategory=15, label_format=60) +
  theme_minimal()
dev.off()
cluster_summary <- as.data.frame(cclust)
write.csv(cluster_summary, file = "DEG_splitKEGGenrich.csv", row.names=FALSE)

# DMR-DEG enrichments
ego <- enrichGO(gene= intersection_genes,
                OrgDb= org.Hs.eg.db,
                keyType = "SYMBOL",
                ont= "BP",
                universe=as.character(degs_df$Symbols)
)

cluster_summary <- as.data.frame(ego)
write.csv(cluster_summary, file = "meDEG_allGO_enrich.csv", row.names=FALSE)

# s_ego <- clusterProfiler::simplify(ego)   # Condensing some terms
tiff("meDEG_GOBPenrich.tiff", units="in", width=7, height=5, res=600)
dotplot(ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of meDEGs (FDR<0.05)")
dev.off()

intersect_keggList <- bitr(intersection_genes, fromType="SYMBOL", toType="ENTREZID", 
                     OrgDb=org.Hs.eg.db)

eKEGG <- enrichKEGG(gene= intersect_keggList$ENTREZID,
                    organism= "hsa",
                    keyType= "ncbi-geneid",
                    universe=all_keggIDs$ENTREZID
)

cluster_summary <- as.data.frame(ego)
write.csv(cluster_summary, file = "meDEG_KEGG_enrich.csv", row.names=FALSE)

dotplot(eKEGG, showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("KEGG Enrichment of meDEGs (FDR<0.05)")

temp <- data.frame(eKEGG)
