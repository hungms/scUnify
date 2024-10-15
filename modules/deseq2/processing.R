convert_deseq2_to_findmarkers <- function(res, x, assay = "RNA"){
    res <- res %>%
        mutate(gene = rownames(.), p_val = pvalue, p_val_adj = padj, avg_log2FC = log2FoldChange) %>%
        dplyr::select(p_val, p_val_adj, avg_log2FC, gene) %>%
        filter(!is.na(p_val_adj))
    gene <- res$gene
    res$pct.1 <- round(rowSums(x[["RNA"]]$counts[gene,, drop = FALSE] > 0) / ncol(x), digits = 3)
    res$pct.2 <- 1
    return(res)}

