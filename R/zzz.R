.onLoad <- function(libname, pkgname) {
  packages <- c(
    "ggplot2", "dplyr", "ggrepel", "patchwork", "cowplot", "scales", 
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony", "ggh4x", 
    "biomaRt", "BiocParallel", "UCell", "celda", "scDblFinder",
    "DoubletFinder", "SeuratDisk", "ComplexHeatmap")
  invisible(lapply(packages, library, character.only = TRUE))
  dir <- system.file("essentials", package = pkgname)
  scripts <- list.files(dir, full.names = TRUE)
  for(script in scripts) {
    source(script)
  }
}
