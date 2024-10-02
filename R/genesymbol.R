bcr.string <- "^I[Gg][HKLhkl][VDJCAEMGvdjcaemg]|^AC233755"
#' BCR string
#'
#' BCR string
#' @export
"bcr.string"

tcr.string <- "^T[Rr][ABCDGabcdg][VDJCvdjc]"
#' TCR string
#'
#' TCR string
#' @export
"tcr.string"

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
   humanx <- as.data.frame(genesV2[,c(2,1)]) %>% distinct(HGNC.symbol, .keep_all = T)
   if(unique){ # disable one to many
      humanx <- unique(mousex$HGNC.symbol)}
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
   mousex <- as.data.frame(genesV2[,c(2,1)]) %>% distinct(MGI.symbol, .keep_all = T)
   if(unique){ # disable one to many
      mousex <- unique(mousex$MGI.symbol)}
   return(mousex)}


#' convert_human_to_mouse_seurat
#' 
#' Make a Seurat assay converting human (HGNC) gene symbols to mouse (MGI) gene symbols ; disable one-to-many mapping
#' 
#' @param x Seurat object
#' @param orig.assay Original assay name
#' @param new.assay New assay name
#' @param unique Whether to disable one-to-many MGI gene symbols to be returned. Defaults to TRUE to disable one-to-many mapping (TRUE)
#' @return If unique = TRUE, returns a vector of unique MGI gene symbols that are mapped one-to-one. 
#'         If unique = FALSE, returns a data frame with columns "MGI.symbol" and "HGNC.symbol" that are mapped one-to-many.
#' @export
convert_human_to_mouse_seurat <- function(x, orig.assay = "RNA", new.assay = "RNA.MM", unique = T){
    DefaultAssay(x) <- orig.assay
    x[[orig.assay]] <- JoinLayers(x[[orig.assay]])
    genes <- convert_human_to_mouse(rownames(x), unique = F)
    if(unique){
        dup <- genes$MGI.symbol[which(duplicated(genes$MGI.symbol))]
        genes <- genes %>%
            filter(!MGI.symbol %in% dup)}
    counts <- as.data.frame(x[[orig.assay]]$counts)[c(genes$HGNC.symbol),]
    rownames(counts) <- genes$MGI.symbol
    x[[new.assay]] <- CreateAssay5Object(counts, min.feature = 0, min.cell = 0)
    return(x)
}

#' convert_mouse_to_human_seurat
#' 
#' Make a Seurat assay converting mouse (MGI) gene symbols to human (HGNC) gene symbols ; disable one-to-many mapping
#' 
#' @param x Seurat object
#' @param orig.assay Original assay name
#' @param new.assay New assay name
#' @param unique Whether to disable one-to-many HGNC gene symbols to be returned. Defaults to TRUE to disable one-to-many mapping
#' @return If unique = TRUE, returns a vector of unique HGNC gene symbols that are mapped one-to-one. 
#'         If unique = FALSE, returns a data frame with columns "MGI.symbol" and "HGNC.symbol" that are mapped one-to-many.
#' @export
convert_mouse_to_human_seurat <- function(x, orig.assay = "RNA", new.assay = "RNA.HS", unique = T){
    DefaultAssay(x) <- orig.assay
    x[[orig.assay]] <- JoinLayers(x[[orig.assay]])
    genes <- convert_mouse_to_human_full(rownames(x), unique = F)
    if(unique){
       dup <- genes$HGNC.symbol[which(duplicated(genes$HGNC.symbol))]
       genes <- genes %>%
          filter(!HGNC.symbol %in% dup)}
    counts <- as.data.frame(x[[orig.assay]]$counts)[c(genes$MGI.symbol),]
    rownames(counts) <- genes$HGNC.symbol
    x[[new.assay]] <- CreateAssay5Object(counts, min.feature = 0, min.cell = 0)
    return(x)
}
