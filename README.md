# scUnify
scUnify is an R package to simplify single-cell omic pipelines and provide customisable visualization functions for multi-omic datasets. Please refer to the [scUnify vignette](https://mshung229.github.io/scUnify/) for example usage and explanations of each function and module. For comprehensive applications of the functions/modules in real single-cell data analysis please refer to the [SCWorkBook documentation](https://github.com/mshung229/scworkbook) for more details.

## Installation
```{r}
remotes::install_github("mshung229/scUnify")
```

## Dependencies
```{r}
## Install CRAN packages
options(timeout = 500)
cran <- c(
    "ggplot2", "dplyr", "ggrepel", "patchwork", "cowplot", "scales", "ggh4x",
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony")
install.packages(cran)

## Install Bioconductor packages
bioc <- c("biomaRt", "BiocParallel", "UCell", "celda", "scDblFinder", "ComplexHeatmap")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioc)

## Install Github packages
git <- c("chris-mcginnis-ucsf/DoubletFinder", "mojaveazure/seurat-disk")
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(git)
```

## Optional Modules
```{r}
## Optional Modules
### Under Development (as of 31 Aug 2024)
```