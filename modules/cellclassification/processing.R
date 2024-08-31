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