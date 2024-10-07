#' dsb_normalize
#' 
#' Normalize CITEseq expressions by DSB normalization
#' 
#' @param x a list of seurat objects
#' @param dir 
#' @param denoise.counts
#' @param use.isotype.control
#' @param isotype.control.name.vec
#' @return a list of seurat objects with DSB normalised "data" in ADT assay
#' @export
dsb_normalize <- function(x, dir, denoise.counts = T, use.isotype.control = F, isotype.control.name.vec = NULL){
    stopifnot(length(x) == length(dir))
    raw.dir <- gsub("/outs/.*", "/outs/multi/count/raw_feature_bc_matrix", dir)

    for(i in seq_along(x)){
        adt.cell <- Read10X(dir[i])[["Antibody Capture"]][rownames(x[["ADT"]]),]
        adt.raw <- Read10X(raw.dir[i])[["Antibody Capture"]][rownames(x[["ADT"]]),]
        adt.raw <- adt.raw[, -c(which(colnames(adt.raw) %in% colnames(adt.cell)))]
        x[[i]][["ADT"]]$data = dsb::DSBNormalizeProtein(
            cell_protein_matrix = adt.cell, 
            empty_drop_matrix = adt.raw, 
            denoise.counts = denoise.counts, 
            use.isotype.control = use.isotype.control,
            isotype.control.name.vec = isotype.control.name.vec)}
    
    return(x)
    }


