pkgs <- c(
    "tidyverse", "ggrepel", "patchwork", "cowplot", "scales", 
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony", "ggh4x", 
    "biomaRt", "BiocParallel", "UCell", "celda", "scDblFinder", 
    "DoubletFinder", "ComplexHeatmap")

for(x in pkgs){
    usethis::use_package(x)}