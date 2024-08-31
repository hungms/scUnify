#' get_biomart
#'
#' Initializes the Ensembl Mart for human and mouse gene symbols.
#' 
#' @param host URL to retrieve archived ensembl biomarts
#' @return set human and mouse biomart objects in local environment
#' @export
get_biomart <- function(host = 'https://dec2021.archive.ensembl.org') {
   assign("human_biomart", biomaRt::useMart("ensembl", host = host, dataset = "hsapiens_gene_ensembl"), envir = .GlobalEnv)
   assign("mouse_biomart", biomaRt::useMart("ensembl", host = host, dataset = "mmusculus_gene_ensembl"), envir = .GlobalEnv)}

#' convert_mouse_to_human
#' 
#' Convert mouse (MGI) gene symbols to human (HGNC) gene symbols; disable one-to-many mapping
#' 
#' @param x A vector of MGI gene symbols.
#' @param unique Whether to disable one-to-many HGNC gene symbols to be returned. Defaults to TRUE to disable one-to-many mapping
#' @return If unique = TRUE, returns a vector of unique HGNC gene symbols that are mapped one-to-one. 
#'         If unique = FALSE, returns a data frame with columns "MGI.symbol" and "HGNC.symbol" that are mapped one-to-many.
#' @export
convert_mouse_to_human <- function(x, unique = T){
   if(!all(c("mouse_biomart", "human_biomart") %in% ls())){
      get_biomart()}
   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse_biomart, attributesL = c("hgnc_symbol"), martL = human_biomart, uniqueRows=T)
   if(unique){ # disable one to many
      humanx <- genesV2[,2][-c(which(duplicated(genesV2[,2])))]}
   else{ # enable one to many
      humanx <- as.data.frame(genesV2[,c(2,1)]) %>% distinct(HGNC.symbol, .keep_all = T)}
   return(humanx)}

#' convert_human_to_mouse
#' 
#' Convert human (HGNC) gene symbols to mouse (MGI) gene symbols; disable one-to-many mapping
#' 
#' @param x A vector of HGNC gene symbols.
#' @param unique Whether to disable one-to-many MGI gene symbols to be returned. Defaults to TRUE to disable one-to-many mapping (TRUE)
#' @return If unique = TRUE, returns a vector of unique MGI gene symbols that are mapped one-to-one. 
#'         If unique = FALSE, returns a data frame with columns "MGI.symbol" and "HGNC.symbol" that are mapped one-to-many.
#' @export
convert_human_to_mouse <- function(x, unique = T){
   if(!all(c("mouse_biomart", "human_biomart") %in% ls(pattern = "_biomart"))){
      get_biomart()}
   genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human_biomart, attributesL = c("mgi_symbol"), martL = mouse_biomart, uniqueRows=T)
   if(unique){ # disable one to many
      mousex <- genesV2[,2][-c(which(duplicated(genesV2[,2])))]}
   else{ # enable one to many
      mousex <- as.data.frame(genesV2[,c(2,1)]) %>% distinct(MGI.symbol, .keep_all = T)}
   return(mousex)}
