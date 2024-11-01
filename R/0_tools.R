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

#' clean_refmap_prediction
#'
#' clean reference mapping outputs
#' @param x Seurat Object
#' @return Seurat Object cleaned
#' @export
clean_refmap_prediction <- function(x){
    if(any(str_detect(names(x@assays), "^prediction.score"))){
        remove <- names(x@assays)[which(str_detect(names(x@assays), "^prediction.score"))]
        for(i in remove){
            x[[i]] <- NULL}}
    return(x)}

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