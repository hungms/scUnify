# Installation


## Install CRAN packages
options(timeout = 500)
cran <- c(
    "ggplot2", "dplyr", "ggrepel", "patchwork", "cowplot", "scales", 
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony")
install.packages(cran)

## Install Bioconductor packages
bioc <- c("biomaRt", "BiocParallel", "UCell", "celda", "scDblFinder")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioc)

## Install Github packages
git <- c("chris-mcginnis-ucsf/DoubletFinder", "mojaveazure/seurat-disk")
devtools::install_github(git)

## Optional Modules
### Under Development (as of 31 Aug 2024)