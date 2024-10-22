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