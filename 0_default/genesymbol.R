human = useMart("ensembl", host = 'https://dec2021.archive.ensembl.org', dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", host = 'https://dec2021.archive.ensembl.org',dataset = "mmusculus_gene_ensembl")

#' convert_mouse_to_human
#' 
#' @description Convert mouse (MGI) gene symbols to human (HGNC) gene symbols; disable one-to-many mapping
#' @param x A vector of MGI gene symbols.
#' @param unique Whether to disable one-to-many HGNC gene symbols to be returned. Defaults to {TRUE} to disable one-to-many mapping
#' @return If {unique = T}; return a list of dataframe with "MGI.symbol", "HGNC.symbol" as mouse and corresponding human gene symbols. Otherwise return human gene symbols with one-to-one mapping only.
#' @export
#' @examples 
#' convert_mouse_to_human(c("Mki67", "Top2a", "Pcna"), unique = T)
#' 
convert_mouse_to_human <- function(x, unique = T){
   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
   if(unique){ # disable one to many
      humanx <- genesV2[,2][-c(which(duplicated(genesV2[,2])))]}
   else{ # enable one to many
      humanx <- as.data.frame(genesV2[,c(2,1)]) %>% distinct(HGNC.symbol, .keep_all = T)}
   return(humanx)}

#' convert_human_to_mouse
#' 
#' @description Convert human (HGNC) gene symbols to mouse (MGI) gene symbols; disable one-to-many mapping
#' @param x A vector of HGNC gene symbols.
#' @param unique Whether to disable one-to-many MGI gene symbols to be returned. Defaults to {TRUE} to disable one-to-many mapping (TRUE)
#' @return If {unique = T}; return a list of dataframe with "HGNC.symbol", "MGI.symbol" as human and corresponding mouse gene symbols. Otherwise return mouse gene symbols with one-to-one mapping only.
#' @export
#' @examples 
#' convert_mouse_to_human(c("MKI67", "TOP2A", "PCNA"), unique = T)
#' 
convert_human_to_mouse <- function(x, unique = T){
      genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                     values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
   if(unique){ # disable one to many
      mousex <- genesV2[,2][-c(which(duplicated(genesV2[,2])))]}
   else{ # enable one to many
      mousex <- as.data.frame(genesV2[,c(2,1)]) %>% distinct(MGI.symbol, .keep_all = T)}
   return(mousex)}
