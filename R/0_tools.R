#' create_seurat_object
#'
#' Create Seurat object from cellranger-multi outputs
#' @param dir a vector of directories
#' @param samples a vector of sample names for each directory
#' @param hto_str prefix to identify HTO tag names if HTO library is present. Defaults to NULL
#' @export
create_seurat_object <- function(dir, samples, hto_str = NULL, adt_normalize = T){
    stopifnot(length(samples) == length(dir))

    x <- list()
    for(i in seq_along(samples)){

       	# read 10x output
        message(paste0(samples[i], " --- Loading Sample ", i))
        message("Step 1 : Adding RNA counts")
        mtx <- Read10X(dir[i])

        # add rna counts
        any.adt <- "Antibody Capture" %in% names(mtx)
        if(any.adt){rna <- mtx$`Gene Expression`}
        else{rna <- mtx}
        seu_obj <- suppressWarnings(CreateSeuratObject(as(rna, "sparseMatrix"), assay="RNA", min.cells = 3, min.features = 200)) # import gene expression/adt to seurat, remove cells with less than 200 features and features expressed in less than 3 cells

        # check for modalities
        message("Step 2 : Check for Modalities")
        if(any.adt){
            any.hto <- any(str_detect(rownames(mtx$`Antibody Capture`), hto_str))
            any.adt <- any(!(str_detect(rownames(mtx$`Antibody Capture`), hto_str)))}
        else{
            any.hto <- F}
        adt.skip <- ifelse(any.adt, "", "--- SKIPPED")
        hto.skip <- ifelse(any.hto, "", " --- SKIPPED")

       	# add hto counts
        message(paste0("Step 3 : Adding HTO counts", hto.skip))
        if(any.hto){
            hash <- mtx$`Antibody Capture`[c(which(str_detect(rownames(mtx$`Antibody Capture`), hto_str))),]
            seu_obj[["HTO"]] <- CreateAssay5Object(as(hash[,colnames(seu_obj[["RNA"]]$counts)], "sparseMatrix"), min.cells = 0, min.features = 0)
            }

       	# add adt counts
        message(paste0("Step 4 : Adding ADT counts", adt.skip))
        if(any.adt){
            prot <- mtx$`Antibody Capture`
            if(any.hto){
                prot <- prot[-c(which(str_detect(rownames(prot), hto_str))),]}
            seu_obj[["ADT"]] <- CreateAssay5Object(as(prot[,colnames(seu_obj[["RNA"]]$counts)], "sparseMatrix"), min.cells = 0, min.features = 0)
            if(adt_normalize){
                seu_obj <- run_dsb(seu_obj, dir = dir[i], denoise.counts = T, use.isotype.control = F, isotype.control.name.vec = NULL)}
            }

       	# add sample id in metadata
        message("Step 5 : Creating Seurat Object")
        seu_obj <- set_assay_keys(seu_obj)
        seu_obj@meta.data$samples <- rep(samples[i], ncol(seu_obj))
        seu_obj <- RenameCells(seu_obj, new.names = paste0(samples[i], "_", colnames(seu_obj)))
        x[[i]] <- seu_obj}

    names(x) <- paste0(samples)
    return(x)}

#' remove_genes
#'
#' Remove genes from current assay to new assay
#' @param x Seurat object
#' @param features a vector genes/features to remove
#' @param orig.assay current assay name, defaults to "RNA"
#' @param new.assay new assay name to store removed genes
#' @export
remove_genes <- function(x, features = NULL, orig.assay = "RNA", new.assay){

    if(orig.assay == "RNA"){
        stopifnot("data" %in% Layers(x, assay = "RNA"))}

    DefaultAssay(x) <- orig.assay
    keep <- which(rownames(x) %in% c(features))
    message(paste0("Detected genes = ", length(keep)))

    if(length(keep) <= 1){
        stop("Detected genes <= 1, stop subsetting genes")}

    feature.counts <- x[[orig.assay]]$counts[keep,]
    new.counts <- x[[orig.assay]]$counts[-c(keep),]

    if("data" %in% Layers(x, assay = orig.assay)){
        feature.data <- x[[orig.assay]]$data[keep,]
        new.data <- x[[orig.assay]]$data[-c(keep),]
        x[[new.assay]] <- CreateAssay5Object(counts = as(feature.counts, "sparseMatrix"), data = as(feature.data, "sparseMatrix"), min.features = 0, min.cells = 0)
        x[[orig.assay]] <- CreateAssay5Object(counts = as(new.counts, "sparseMatrix"), data = as(new.data, "sparseMatrix"), min.features = 0, min.cells = 0)}
    else{
        x[[new.assay]] <- CreateAssay5Object(counts = as(feature.counts, "sparseMatrix"), min.features = 0, min.cells = 0)
        x[[orig.assay]] <- CreateAssay5Object(counts = as(new.counts, "sparseMatrix"), min.features = 0, min.cells = 0)}

    return(x)}

#' return_genes
#'
#' Return genes from an assay to a specific assay
#' @param x Seurat object
#' @param from.assay assay name from which genes/features are taken
#' @param to.assay assay name to which genes/features are stored, defaults to "RNA"
#' @export
return_genes <- function(x, from.assay, to.assay = "RNA"){
    counts <- rbind(x[[from.assay]]$counts, x[[to.assay]]$counts)
    data <- rbind(x[[from.assay]]$data, x[[to.assay]]$data)
    x[[to.assay]] <- CreateAssay5Object(counts = counts, data = data, min.cells = 0, min.features = 0)
    x[[from.assay]] <- NULL
    return(x)}

#' remove_vdj_genes
#'
#' Remove VDJ genes from an assay
#' @param x Seurat object
#' @param bcr if TRUE, remove BCR-VDJ genes. Defaults to FALSE
#' @param tcr if TRUE, remove TCR-VDJ genes. Defaults to FALSE
#' @param orig.assay current assay name, defaults to "RNA"
#' @export
remove_vdj_genes <- function(x, bcr = T, tcr = T, orig.assay = "RNA"){
    if(bcr){
        bcr.genes <- rownames(x[[orig.assay]])[which(str_detect(rownames(x[[orig.assay]]), bcr.string))]
        x <- remove_genes(x, features = bcr.genes, orig.assay = orig.assay, new.assay = "BCR")}
    if(tcr){
        tcr.genes <- rownames(x[[orig.assay]])[which(str_detect(rownames(x[[orig.assay]]), tcr.string))]
        x <- remove_genes(x, features = tcr.genes, orig.assay = orig.assay, new.assay = "TCR")}
    return(x)
        }
        
#' set_assay_keys
#'
#' Bugfix for setting assay keys when creating Seurat object
#' @param x Seurat object
#' @export
set_assay_keys <- function(x){
    for(assay in names(x@assays)){
        x[[assay]]@key <- paste0(tolower(assay), "_")}
    return(x)}

#' join_layers
#'
#' Join layers for multiple assays 
#' @param x Seurat object
#' @param assays A vector of assay names to join. Defaults to {assays = c("RNA", "BCR", "TCR", "CC", "ADT", "HTO")}
#' @return Seurat object
#' @export
join_layers <- function(x, assays = c("RNA", "BCR", "TCR", "CC", "ADT", "HTO")){
    assays.to.join <- names(x@assays)[which(names(x@assays) %in% assays)]
    for(i in assays.to.join){
        x[[i]] <- JoinLayers(x[[i]])}
    return(x)}

#' list_to_df
#'
#' Convert list to dataframe
#' @param list A list of vectors
#' @return A dataframe where each column is a vector from the list
#' @export
list_to_df <- function(list){
    max_length <- max(sapply(list, length))
    padded.list <- lapply(list, function(v) {
        c(v, rep("", max_length - length(v)))})
    df <- as.data.frame(padded.list, stringsAsFactors = FALSE)
    return(df)}

#' metadata_checkpoint
#'
#' add checkpoint for metadata
#' @param x Seurat Object
#' @param dir directory storing metadata
#' @param save save metadata as tsv
#' @param load specific metadata.tsv to load
#' @return Seurat Object with updated metadata file
#' @export
metadata_checkpoint <- function(x, dir = "seurat/0_metadata/", save = F, load = NULL){
    if(save){
        date <- gsub("-", "", Sys.Date())
        write.table(x@meta.data, paste0(dir, "/", date, "_metadata.tsv"), quote = F, row.names = T, sep = "\t")}

    if(length(load) == 0){
        file <- list.files(dir)[length(list.files(dir))]
        latest <- paste0(dir, "/", file)}
    else if(file.exists(load)){
        latest <- load}
    else{stop("file does not exist")}
    
    x@meta.data <- read.table(latest, sep = "\t")[colnames(x),]
    return(x)
}

#' calculate_10x_dbr
#'
#' calculate doublet rate for 10x sequencing runs
#' @param ncells number of cells captured
#' @return a vector of doublet rates
#' @export
calculate_10x_dbr <- function(ncells){
    
    cells <- c(500, seq(1000, 10000, 1000))
    percent <- c(0.004, seq(0.008, 0.08, 0.008))
    dbrate <- data.frame(cells, percent)
    model <- lm(percent ~ poly(cells, 3, raw = TRUE), data = dbrate)
    p <- data.frame(cells = ncells)
    predictions <- model %>% predict(p)

    return(predictions)
}

#' intersect_genes
#'
#' intersect common gene names between a list of seurat objects
#' @param x a list of Seurat object
#' @param assay assay name
#' @param min.cells minimum number of cells expressing the gene
#' @export
intersect_genes <- function(x, assay = "RNA", min.cells = 3){
    intersect.list <- list()
    for(i in seq_along(x)){
        keep <- rowSums(x[[i]][[assay]]$counts >= 1) >= min.cells
        intersect.list[[i]] <- rownames(x[[i]][[assay]]$counts)[keep]}
    genes <- Reduce(intersect, intersect.list)
    return(genes)}
