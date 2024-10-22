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

#' select_top_deg
#'
#' Find top differentially expressed genes for each cluster from Seurat::FindAllMarkers output.
#' @param markers a dataframe of Seurat::FindAllMarkers output
#' @param n no. of genes per cluster
#' @param rank_by variables to rank by, either "avg_log2FC", "p_val_adj", "diff_pct"
#' @param only.pos if TRUE, return only genes with positive avg_log2FC. Defaults to FALSE
#' @export
select_top_deg <- function(markers, n = 5, rank_by = "avg_log2FC", only.pos = F, min.pct = 0.1){

    if(!rank_by %in% c("avg_log2FC", "p_val_adj", "diff_pct", "pxFC")){
        stop('please make sure rank_by is either "avg_log2FC", "p_val_adj", "diff_pct", "pxFC"')}

    genelist <- markers %>%
        filter(!str_detect(gene, paste0("^Mt|^Rp[sl]|", bcr.string))) %>%
        filter(p_val_adj < 0.05) %>%
        mutate(diff_pct = pct.1 - pct.2)

    if(only.pos){
        genelist <- genelist %>%
            filter(avg_log2FC > 0 & pct.1 >= min.pct)}
    else{
        genelist <- genelist %>%
            mutate(
                direction = ifelse(avg_log2FC > 0, "UP", "DOWN"),
                cluster = paste0(cluster, "_", direction),
                remove = case_when(
                    direction == "UP" & pct.1 < min.pct ~ "remove",
                    direction == "DOWN" & pct.1 == 0 ~ "remove",
                    .default = "keep")) %>%
            filter(remove != "remove")}

    if(rank_by == "avg_log2FC"){
        genelist <- genelist %>%
            group_by(cluster) %>%
            slice_max(n=n, order_by = avg_log2FC^2, with_ties = F) %>%
            split(.$cluster) %>%
            lapply(function(df) df$gene)}

    if(rank_by == "p_val_adj"){
        genelist <- genelist %>%
            group_by(cluster) %>%
            slice_min(n=n, order_by = p_val_adj, with_ties = F) %>%
            split(.$cluster) %>%
            lapply(function(df) df$gene)}
    
    if(rank_by == "diff_pct"){
        genelist <- genelist %>%
            group_by(cluster) %>%
            slice_max(n=n, order_by = diff_pct^2, with_ties = F) %>%
            split(.$cluster) %>%
            lapply(function(df) df$gene)}

    if(rank_by == "pxFC"){
        genelist <- genelist %>%
            mutate(pxFC = avg_log2FC*-log10(p_val_adj)) %>%
            group_by(cluster) %>%
            slice_max(n=n, order_by = pxFC^2, with_ties = F) %>%
            split(.$cluster) %>%
            lapply(function(df) df$gene)}

    return(genelist)}

#' convert_deseq2_to_findmarkers
#'
#' 
#' @export
convert_deseq2_to_findmarkers <- function(res, x, assay = "RNA"){
    res <- res %>%
        mutate(gene = rownames(.), p_val = pvalue, p_val_adj = padj, avg_log2FC = log2FoldChange) %>%
        dplyr::select(p_val, p_val_adj, avg_log2FC, gene) %>%
        filter(!is.na(p_val_adj))
    gene <- res$gene
    res$pct.1 <- round(rowSums(x[["RNA"]]$counts[gene,, drop = FALSE] > 0) / ncol(x), digits = 3)
    res$pct.2 <- 1
    return(res)}