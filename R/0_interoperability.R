#' convert_seurat_to_anndata
#'
#' Convert Seurat to AnnData with SeuratDisk::Convert()
#' @param x Seurat object
#' @param h5ad path for output .h5ad file
#' @param columns a vector of metadata columns to keep. Defaults to NULL to keep all columns
#' @param pca reduction name for X_pca. Defaults to "pca"
#' @param umap reduction name for X_umap. Defaults to "umap"
#' @param assay assay name for gene expression. Defaults to "RNA"
#' @param adt_assay assay name for ADT counts if any
#' @param return_genes if TRUE, return genes from assays named c("CC", "BCR", "TCR", "MHC") to current assay. Defaults to TRUE
#' @export
convert_seurat_to_anndata <- function(x, h5ad, columns = NULL, umap = "umap", pca = NULL, snn = NULL, nn = NULL, assay = "RNA", adt_assay = NULL, return_genes = T){
    stopifnot(all(c(umap, pca) %in% names(x@reductions)))
    stopifnot(all(c(snn, nn) %in% names(x@graphs)))

    # downgrade to v4
    options(Seurat.object.assay.version = "v3")
    if(!is.object(x)){
        if(file.exists(x)){
            x <- qread(paste0(x))}}
    
    # preselect columns if necessary
    if(!length(columns) > 0){
        message(paste0("selecting all columns from metadata"))
        columns <- colnames(x@meta.data)}
    else{
	    message(paste0("selecting columns from metadata: ", paste0(columns, collapse = ", ")))
        x@meta.data <- x@meta.data[,columns]}
    
    # convert non-numeric columns to character class
    classlist <- c()
    for(i in seq_along(colnames(x@meta.data))){
        if(is.integer(x@meta.data[[i]]) | is.numeric(x@meta.data[[i]])){
            x@meta.data[[i]] <- ifelse(is.na(x@meta.data[[i]]), 0, x@meta.data[[i]])
            x@meta.data[[i]] <- as.numeric(x@meta.data[[i]])}
        if(!is.numeric(x@meta.data[[i]])){
            x@meta.data[[i]] <- ifelse(is.na(x@meta.data[[i]]), "NA", x@meta.data[[i]])
            x@meta.data[[i]] <- as.character(x@meta.data[[i]])}}

    
    # add adt features to metadata
    if(length(adt_assay) > 0){
        stopifnot(adt_assay %in% names(x@assays))
        adt_features <- paste0(tolower(adt_assay), "_", rownames(x[[adt_assay]]))
        adt <- FetchData(x, vars = adt_features, slot = "counts")
        x@meta.data[colnames(adt)] <- NULL
        x@meta.data <- cbind(x@meta.data, adt)}
    
    # return CC, BCR, TCR, MHC genes back to assay slot
    if(return_genes & !str_detect(assay, "SCT|integrated")){
        return_assays <- intersect(c("CC", "BCR", "TCR", "MHC"), names(x@assays))
        if(length(return_assays) > 0){
            message(paste0("returning assays (", paste0(return_assays, collapse = ", "), ") to ", assay, " assay"))
            for(i in return_assays){
                x <- return_genes(x, from.assay = i, to.assay = assay)}}}

    # convert assay to V3
    for(i in names(x@assays)){
        if(!str_detect(i, "SCT|integrated")){
            x[[i]] <- as(object = x[[i]], Class = "Assay")}}

    # store PCA/UMAP/neighbours
    #message("storing PCA, UMAP to the appropriate reduction slot")
    DefaultAssay(x) <- assay
    x@reductions[["umap"]] <- x@reductions[[umap]]
    x@reductions[["pca"]] <- x@reductions[[pca]]
    x@reductions <- x@reductions[c("umap", "pca")]
    for(i in c("umap", "pca")){
        x@reductions[[i]]@assay.used <- assay}

    if(length(snn) > 0 & length(nn) > 0){
        x@graphs[[paste0(assay, "_snn")]] <- x@graphs[[snn]]
        x@graphs[[paste0(assay, "_nn")]] <- x@graphs[[nn]]}
    
    MuDataSeurat::WriteH5AD(x, h5ad, assay=assay, scale.data = F) #https://github.com/zqfang/MuDataSeurat #remotes::install_github("zqfang/MuDataSeurat", dependencies = F)
    options(Seurat.object.assay.version = "v5")
    }

#' convert_seurat_to_sce
#'
#' Convert Seurat to SingleCellExperiment object with SeuratWrapper::as.SingleCellExperiment()
#' @param x Seurat object
#' @param assay assay name for gene expression. Defaults to "RNA"
#' @param return_genes if TRUE, return genes from assays named c("CC", "BCR", "TCR", "MHC") to current assay. Defaults to TRUE
#' @export
convert_seurat_to_sce <- function(x, assay = "RNA", return_genes = T){

    if(return_genes & !str_detect(assay, "SCT")){
        return_assays <- intersect(c("CC", "BCR", "TCR", "MHC"), names(x@assays))
        if(length(return_assays) > 0){
            message(paste0("returning assays (", paste0(return_assays, collapse = ", "), ") to ", assay, " assay"))
            for(i in return_assays){
                x <- return_genes(x, from.assay = i, to.assay = assay)}}}

    x <- DietSeurat(x, assay = assay)
    x[[assay]] <- as(object = x[[assay]], Class = "Assay")
    sce <- as.SingleCellExperiment(x)
    return(sce)}
