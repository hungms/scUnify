#' convert_genedict
#'
#' Convert a list of gene symbols between mouse MGI & human HGNC
#' @param genelist a list of gene symbols
#' @param to organism for the returned gene symbols, either "human" or "mouse"
#' @export
convert_genedict <- function(genelist, to = "human"){
    gs <- list()
    for(l in seq_along(genelist)){
        if(names(genelist)[l] %in% c("Quality","Cycling")){
            gs[[l]] <- genelist[[l]]}
        else{
            if(to == "human"){
                gs[[l]] <- convert_mouse_to_human(genelist[[l]])}
            if(to == "mouse"){
                gs[[l]] <- convert_human_to_mouse(genelist[[l]])}}}
    names(gs) <- names(genelist)
    return(gs)}

#' UCell
#'
#' @export
UCell <- function(x, assay = "RNA"){
    dir = "/camp/home/hungm/scratch/hungm/reference/cell_signatures/LZDZ_signatures/"
    signature_list <- list.files(dir)
    gc_signatures <- list()
    for(x in seq_along(signature_list)){
        gc_signatures[[x]] <- read.csv(paste0(dir, signature_list[x]), sep = ",", header = F)[[1]]}
    names(gc_signatures) <- gsub(".txt", "", signature_list)
    DefaultAssay(query.mapped.hq) <- assay
    x <- AddModuleScore_UCell(x, features=gc_signatures)}