.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
    "tidyverse", "ggrepel", "patchwork", "cowplot", "scales", 
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony", "ggh4x", 
    "biomaRt", "BiocParallel", "UCell", "scDblFinder",
    "DoubletFinder", "ComplexHeatmap")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  options(Seurat.object.assay.version = "v5")
  options(future.globals.maxSize = 1e12)
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
  errorMessage <- NULL

  ### load essential scripts
  #dir <- system.file("essentials", package = pkgname)
  #scripts <- list.files(dir, full.names = TRUE)
  #for(script in scripts){
  #  source(script)}

}
