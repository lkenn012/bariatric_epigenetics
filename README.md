# Description - Epigenetic reprogramming of skeletal muscle by bariatric surgery
R code and data for analysis and visualization of integrated transcriptomics & DNA methylation analysis in bariatric surgery
Forthcoming paper. Harper lab, submitted.
🧬📄 Coming soon

## Code
Contains all R code used in the analysis of Salmon-quantified RNAseq read count and EPIC array IDAT data.

## Data
All files used in the analysis (Coming soon)

## Environment information
The following version of R along with all accompanying libraries are used in the code provided here. All code was implemented on Windows with an installation time of approximately 15 minutes, and a total runtime of ~1.5 hours (+5 hours if including all enrichment analyses). Times dependent on workstation specifications.

R  --4.4.3
RStudio  --2024.12.1+563

### R libraries
- tidyverse --2.0.0
- minfi --1.52.1
- IlluminaHumanMethylationEPICv2anno.20a1.hg38 --1.0.0
- limma --3.62.2
- DMRcate --3.2.1
- missMethyl --1.40.3
- tximeta --1.24.0
- DESeq2 --1.46.0
- org.Hs.eg.db --3.20.0
- ggplot2 --3.5.2
- ggrepel --0.9.6
- tidyverse --2.0.0
- clusterProfiler --4.14.6
- ComplexHeatmap --2.22.0
- BiocVersion --3.20.0 (for Bioconductor, version-specific package installation)
- remotes -- 2.5.0 (for non-Bioconductor, version-specific package installation)

Installation information:
RStudio installation download from https://posit.co/download/rstudio-desktop/
```
install.packages("remotes")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("BiocVersion", version = "3.20")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
BiocManager::install("limma")
BiocManager::install("DMRcate")
BiocManager::install("missMethyl")
BiocManager::install("tximeta")
BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db")
install_version("ggplot2", version = "3.5.2")
install_version("ggrepel", version = "0.9.6")
install_version("tidyverse", version = "2.0.0")
BiocManager::install("clusterProfiler")
BiocManager::install("ComplexHeatmap")
```
