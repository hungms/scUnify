# scUnify
scUnify is an R package to simplify single-cell omic pipelines and provide customisable visualization functions for multi-omic datasets. Please refer to the [scUnify : Quick Start](https://hungms.github.io/scUnify/) for example usage and explanations of each function and module.  
  
For comprehensive applications of the functions/modules in real single-cell data analysis please refer to the [SCWorkBook](https://github.com/hungms/scworkbook) documentation for more details.

## Installation
```ruby
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
remotes::install_github("hungms/scUnify", dependencies=FALSE)
```

## Core dependencies
```ruby
### Install CRAN packages
options(timeout = 500)
cran <- c(
    "tidyverse", "ggrepel", "patchwork", "cowplot", "scales", "ggh4x",
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony", "dsb")
install.packages(cran)

### Install Bioconductor packages
bioc <- c("biomaRt", "BiocParallel", "UCell", "celda", "scDblFinder", "ComplexHeatmap")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioc)

### Install Github packages
git <- c("chris-mcginnis-ucsf/DoubletFinder", "hungms/MuDataSeurat")
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(git)
```

## Dependencies for additional modules :   

---
### VDJ
```ruby
## Install CRAN packages
cran <- c()
install.packages(cran)

## Install Bioconductor packages
bioc <- c()
BiocManager::install(bioc)

## Install Github packages
git <- c()
devtools::install_github(git)
```

### Cell Classification
```ruby
## Install CRAN packages
cran <- c()
install.packages(cran)

## Install Bioconductor packages
bioc <- c()
BiocManager::install(bioc)

## Install Github packages
git <- c()
devtools::install_github(git)
```

### Differential Expression & Pathway Analysis
```ruby
## Install CRAN packages
cran <- c()
install.packages(cran)

## Install Bioconductor packages
bioc <- c()
BiocManager::install(bioc)

## Install Github packages
git <- c()
devtools::install_github(git)
```