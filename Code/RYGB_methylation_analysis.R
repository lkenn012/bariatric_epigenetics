### Luke Kennedy (Oct, 2025)
### Differential methylation analysis for Pileggi, Mohottalage et al. (2025)
##
## This code imports methylation data (EPIC array, iDAT files) for analysis and 
## visualization.


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
library(stats)

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

# 1. Map genomic regions for corresponding genes
gset <- makeGenomicRatioSetFromMatrix(
  mat = d3,
  array= "IlluminaHumanMethylationEPICv2",
  annotation = "20a1.hg38",
  what = "M"
)

#####
##    Determine differentially methylated probes (DMPs) - minfi
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

write.csv(results_annotated, "DMPs_Annots_noY.csv", row.names = TRUE)

# Additionally, we want to evaluate overall methylation status between the two
# groups.

# For simplicity, we know columns are ordered in pre,post for each patient
preM_matrix <- d3[, seq(1, ncol(d3), 2)]
postM_matrix <- d3[, seq(2, ncol(d3), 2)]

preB_matrix <- d2[, seq(1, ncol(d2), 2)]
postB_matrix <- d2[, seq(2, ncol(d2), 2)]

# Compute individual-level methylation levels
indiv_meth <-data.frame()
for (i in 1:ncol(preM_matrix)) {
  pre_M <- preM_matrix[,i]
  post_M <- postM_matrix[,i]
  pre_B <- preB_matrix[,i]
  post_B <- postB_matrix[,i]
  
  indiv_meth <- rbind(indiv_meth, c(mean(pre_M), mean(post_M), mean(pre_B), 
                                    mean(post_B), median(pre_M), 
                                    median(post_M), median(pre_B),
                                    median(post_B))
  )
}

print(indiv_meth)
rownames(indiv_meth) <- levels(patient_IDs)
colnames(indiv_meth) <- c("mean pre-M", "mean post-M", "mean pre-B", 
                          "mean post-B", "median pre-M", "median post-M", 
                          "median pre-B", "median post-B"
                          )
write.csv(indiv_meth, "overallMeth_noY_perPatient.csv", row.names = TRUE)



#####
## 3. Determine differentially methylated regions (DMRs)
#####
myannotation <- cpg.annotate("array", object=d3, what = "M",
                             arraytype = "EPICv2", epicv2Filter = "mean",
                             epicv2Remap = TRUE, analysis.type="differential",
                             design=design, coef="conditionPost Sx", fdr=1)

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg38")

write.csv(results.ranges, "DMRs_Annots_noY.csv", row.names = TRUE)

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

# NOTE: missMethyl enriochment analysis is quite time-consuming

# Define significant DMRs
sig_DMRs <- results.ranges[results.ranges$no.cpgs > 1 & results.ranges$min_smoothed_fdr < 0.05]

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

# Repeat using only promoter region DMRs
gst.region <- goregion(sig_DMRs, all.cpg=rownames(d3), 
                       collection="GO", array.type="EPIC_V2", plot.bias=TRUE,
                       genomic.features=c("TSS200","TSS1500", "1stExon"), 
                       sig.genes="TRUE")
enrich_sets <- topGSA(gst.region, n=250)
enrich_sets$logFDR <- -1*log10(enrich_sets$FDR)

tiff("misMethyl_GOenrich_promoter_noY.tiff", units="in", width=8, height=4, res=600)
ggplot(data=enrich_sets[1:20,], aes(x=logFDR, y=reorder(TERM,logFDR), size=DE)) +
  geom_point(alpha=0.5, shape=21, color="black", fill="darkslateblue") +
  scale_size(range = c(1, 8), name="Gene count") +
  theme(legend.position="bottom") +
  ylab("GO term") +
  xlab("-log10 FDR") +
  ggtitle("GO Enrichment promoter DMRs (FDR < 0.05)") +
  theme_classic()

dev.off()

# Visualize distribution of methylated regions for genes through a simple barplot
sig_DMRs <- sig_DMRs %>%
  direction = ifelse(maxdiff > 0,"pos", "neg")

up_positions <- sig_DMRs[sig_DMRs[,direction == "pos"], c("start", "end")]

# Identify CpG sites constituting DMRs
locs <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations
locs.ranges <- GRanges(locs$chr, IRanges(locs$pos, locs$pos))
names(locs.ranges) <- rownames(locs)

# Find the overlaps for all gene-associated DMRs
gene_ranges <- results.ranges[!is.na(results.ranges$overlapping.genes), ]
overlaps <- findOverlaps(locs.ranges, gene_ranges)

# Get CpG & DMR overlaps
cpg_indices <- queryHits(overlaps)
dmr_indices <- subjectHits(overlaps)

# Split overlapping CpGs (locs.ranges[cpg_indices]) by DMR index (dmr_indices)
constituent_cpgs <- split(locs.ranges[cpg_indices], dmr_indices)

# Get all CpG site positions for each DMR as a new column in our DMR output
list_of_start_positions <- lapply(constituent.cpgs, start)
cpg_positions_aggregated_string <- sapply(
  list_of_start_positions, function(positions_vector) {
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
  all_annotations <- gsub("\\bexon_1\\b", "1st Exon", all_annotations) # exon 1
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

print("--- Combined Annotation Counts Data ---")
print(combined_counts_df)  # Plotting of these data was done in GraphPad

# save the results
# save the results
dmr_regionsDF <- sig_results %>%
  as_tibble() %>%
  mutate(
    upDMR = (row_number() %in% up_idxs),
    gene_regions = vector("list", length = n()),
    gene_regions = replace(gene_regions, upDMR, unique_upDMR_pos),
    gene_regions = replace(gene_regions, !upDMR, unique_downDMR_pos)
  ) %>%
  mutate(
    gene_regions = map_chr(gene_regions, ~paste(.x, collapse = "; "))
  ) %>%
  select(-upDMR)
write.csv(dmr_regionsDF, "sigDMR_noY_positions.csv")

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

pwlCorPlt_df$abs_pos <- pwlCorPlt_df$pos
pwlCorPlt_df$pos <- 1:nrow(pwlCorPlt_df)

# Calculate chromosome midpoints for x-axis labels
axis_df <- pwlCorPlt_df %>%
  group_by(chr) %>%
  summarise(
    chr_start = min(pos, na.rm = TRUE),
    chr_end = max(pos, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(chr_mid = (chr_start + chr_end) / 2)


# Plot correlation values
tiff("manh_corrPWL.tiff", units="in", width=8, height=6, res=600)
ggplot(pwlCorPlt_df, aes(x=pos, y=absRho)) +
  
  # Show all points (non-highlighted)
  geom_point(aes(color=chr), alpha=0.7, size=0.4) +
  scale_color_manual(values = rep(c("grey", "lightgoldenrod3"), 
                                  length(levels(pwlCorPlt_df$chr)))) +
  
  geom_label_repel(data = subset(pwlCorPlt_df, is_annotateCorr == "yes"),
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
  geom_point(data=subset(pwlCorPlt_df, is_highlightCorr=="Positive"), 
             aes(x=pos, y=absRho), 
             color="dodgerblue3", size=0.7) +
  geom_point(data=subset(pwlCorPlt_df, is_highlightCorr=="Negative"), 
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
