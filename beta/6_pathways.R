#library(clusterProfiler)
#library(msigdbr)
#library(fgsea)

for(i in c("hs_gs", "mm_gs")){
    if(i == "hs_gs"){
        org = "Homo sapiens"}
    else if(i == "mm_gs"){
        org = "Mus musculus"}

    gs <- list()

    # hallmark gene set
    HM<-msigdbr(species = org, category ="H")
    gs[[1]] <- HM %>% dplyr::select(., gs_name, gene_symbol) %>% mutate(gs_name = gsub("HALLMARK_", "",gs_name))

    # go gene set
    GO<-msigdbr(species = org, category ="C5")
    gs[[2]] <-GO %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'GOBP')) %>% mutate(gs_name = gsub("GOBP_", "",gs_name))

    # kegg gene set
    C2<-msigdbr(species = org, category ="C2")
    gs[[3]] <-C2 %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'KEGG')) %>% mutate(gs_name = gsub("KEGG_", "",gs_name))

    # reactome gene set
    gs[[4]] <-C2 %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'REACTOME')) %>% mutate(gs_name = gsub("REACTOME_", "",gs_name))

    # reactome gene set
    gs[[5]] <-C2 %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'BIOCARTA')) %>% mutate(gs_name = gsub("BIOCARTA_", "",gs_name))

    # regulatory gene set
    TFT <-msigdbr(species = org, category ="C3")
    gs[[6]] <-TFT %>% dplyr::select(., gs_name, gene_symbol)

    names(gs) <- c("HM", "GO", "KEGG", "REACTOME", "BIOCARTA", "TFT")
    assign(i, gs)}


# seurat 1to1
FindPathways <- function(x, outdir, org = "human", rank = "avg_log2FC"){

    stopifnot(rank %in% c("avg_log2FC", "p_val_adj", "pxFC"))

    if(org == "human"){
        gs <- hs_gs}
    else if(org == "mouse"){
        gs <- mm_gs}

    x <- x %>%
        mutate(pxFC = -log10(p_val_adj)*avg_log2FC)
    deglist <- x %>%
        arrange(desc(.data[[rank]])) %>%
        pull(rank)
    names(deglist) <- x %>%
        arrange(desc(.data[[rank]])) %>%
        rownames(.)

    deglist <- na.omit(deglist)
    deglist <- deglist[which(deglist != 0)]
    deglist = sort(deglist, decreasing = TRUE)
    output <- list()
    for(i in seq_along(gs)){
        output[i] <- GSEA(deglist, TERM2GENE = gs[[i]], pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 1)}
    names(output) <- paste0("comparison_", names(gs))
    qsave(output, file = outdir)}

# plot seurat 1to1
PlotDotPlot <- function(x, collection, group1, group2, top_n = NULL, filter_signif = T, remove_gs_name = F, ...){

    data <- x[which(str_detect(names(x), paste0("_", collection)))][[1]]@result
    data$NES <- as.numeric(data$NES)
    data$qvalue <- as.numeric(data$qvalue)

    if(remove_gs_name == T){
        data$ID <- sub("^[^_]*_", "", data$ID)}

    if(length(top_n) > 0){
        plot_pathways <- data %>%
                    arrange(desc(NES^2)) %>%
            mutate(direction = ifelse(NES > 0, "POS", "NEG")) %>%
            group_by(direction) %>%
            slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
            .$ID

        if(filter_signif){
            plot_pathways <- data %>%
                filter(qvalue < 0.05) %>%
                arrange(desc(NES^2)) %>%
                mutate(direction = ifelse(NES > 0, "POS", "NEG")) %>%
                group_by(direction) %>%
                slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
                .$ID
        }

        data <- data %>%
            filter(ID %in% plot_pathways)}

    comparison <- paste0("Enriched in ", c(group2, group1))
    plot <- data %>%
        mutate(direction = factor(ifelse(NES < 0, comparison[1], comparison[2]), comparison)) %>%
        mutate(psig = case_when(
            qvalue < 0.001 ~ "***",
            qvalue < 0.01 ~ "**",
            qvalue < 0.05 ~ "*",
            .default = "")) %>%
        ggplot(aes(x = fct_reorder(ID, NES), y = NES, color = NES, size = -log10(qvalue))) +
        geom_point(color = "black", aes(size = (-log10(qvalue))*1.1)) +
        geom_point() +
        geom_text(aes(label = psig), color = "black", hjust = 0.5, nudge_x = 0.5) +
        xlab("Pathways") +
        facet_wrap(~direction, scales = "free_x") +
        theme_bw() +
        scale_color_distiller(palette = "RdBu") +
        scale_size(range = c(1, 6)) +
        guides(
            color = guide_colorbar(
                title = "NES",
                title.position = "top",
                direction = "horizontal",
                frame.colour = "black",
                ticks.colour = "black",
                order = 1),
            size = guide_legend(
                title = "-log10 (padj)",
                direction = "horizontal",
                order = 2,
                override.aes = list(fill = "black"),
                title.position = "top")) +
        coord_flip() +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face="bold", size=12),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        scale_x_discrete(expand=c(0.05, 0.05)) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
    return(plot)
}

# seurat 1toall
FindAllPathways <- function(x, outdir, org = "human", rank = "avg_log2FC"){
    stopifnot(rank %in% c("avg_log2FC", "p_val_adj", "diff_pct", "pxFC"))
    if(org == "human"){
        gs <- hs_gs}
    else if(org == "mouse"){
        gs <- mm_gs}

    levels <- sort(unique(x$cluster))
    output <- list()

    if(rank == "pxFC"){
        x <- x %>%
            mutate(
                pxFC = avg_log2FC*p_val_adj,
                diff_pct = pct.1 - pct.2)}

    for(i in seq_along(levels)){
        deglist <- x %>%
            filter(cluster %in% levels[i]) %>%
            arrange(desc(.data[[rank]])) %>%
            .$avg_log2FC
        names(deglist) <- x %>%
            filter(cluster %in% levels[i]) %>%
            arrange(desc(.data[[rank]])) %>%
            .$gene
        deglist <- na.omit(deglist)
        deglist <- deglist[which(deglist != 0)]
        deglist = sort(deglist, decreasing = TRUE)

        for(j in seq_along(gs)){
            output[[length(output) + 1]] <- GSEA(deglist, TERM2GENE = gs[[j]], pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 1)
            names(output)[length(output)] <- paste0(levels[i], "_", names(gs)[j])}}

    qsave(output, file = outdir)}


# plot pathway heatmaps
PlotAllHeatmap <- function(x, collection, top_n = NULL, only_pos = T, plot_qval = F, remove_gs_name = F, fill_na = F, keep_order = F, signif_col = "white",...){

    x <- x[which(str_detect(names(x), paste0("_", collection)))]

    if(keep_order){
        colname_order <- gsub(paste0("_", collection), "", names(x))}

    for(i in seq_along(x)){
        x[[i]]@result$cluster <- gsub(paste0("_", collection), "", names(x)[i])
        x[[i]] <- x[[i]]@result}

    x <- bind_rows(x)

    if(remove_gs_name == T){
        x$ID <- sub("^[^_]*_", "", x$ID)}

    if(length(top_n) > 0){
        if(only_pos){
            plot_pathways <- x %>%
                group_by(cluster) %>%
                filter(NES > 0) %>%
		arrange(desc(NES)) %>%
                slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
                .$ID}
        else{
            plot_pathways <- x %>%
                group_by(cluster) %>%
		arrange(desc(NES^2)) %>%
                slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
                .$ID}

        x <- x %>%
            filter(ID %in% plot_pathways)}

    nes <- x %>%
        dplyr::select(cluster, ID, NES) %>%
        arrange(cluster) %>%
        mutate(cluster = as.factor(cluster), NES = as.numeric(NES)) %>%
        complete(cluster, nesting(ID), fill = list(NES = NA)) %>%
        pivot_wider(names_from = cluster, values_from = NES) %>%
        column_to_rownames("ID")

    qval <- x %>%
        dplyr::select(cluster, ID, qvalue) %>%
        arrange(cluster) %>%
        mutate(cluster = as.factor(cluster), qvalue = as.numeric(qvalue)) %>%
        complete(cluster, nesting(ID), fill = list(qvalue = 1)) %>%
        pivot_wider(names_from = cluster, values_from = qvalue) %>%
        column_to_rownames("ID")

    if(fill_na == T){
        nes[is.na(nes)] <- 0}

    if(plot_qval == F){
        qval[qval != 1] <- 1}

    min <- round(min(nes, na.rm = TRUE))
    max <- round(max(nes, na.rm = TRUE))
    scale = max(c(-min, max))
    col_fun = rev(brewer.pal(12,"RdBu"))

    if(keep_order){
        ht <- Heatmap(
            nes,
            border = TRUE,
            col = col_fun,
            na_col = "black",
            border_gp = gpar(col = "black", lwd = 3),
            rect_gp = gpar(col = "white", lwd = 1),
            heatmap_legend_param = list(
                title = "NES",
                at = c(-scale, scale),
                legend_direction = "vertical",
                labels = c(paste0("-", scale), paste0(" ", scale)),
                width = unit(10, "mm")),
            show_heatmap_legend = T,
            show_column_names = T,
            column_names_side = "top",
            column_title_side = "top",,
            column_dend_reorder = F,
            column_order = colname_order,
            column_dend_side = "bottom",
            show_row_names = T,
            row_title = "Pathways",
            row_dend_reorder = T,
            row_dend_side = "right",
            row_names_side = "left",
            row_names_max_width = unit(30, "cm"),
            column_title_gp = gpar(fontsize = 10,fontface="bold"),
            column_names_gp = gpar(fontsize = 10,fontface="bold"),
            row_title_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 8),
            cell_fun = function(j, i, x, y, w, h, fill) {
                if(qval[i, j] < 0.001) {
                    grid.text("***", x, y, gp = gpar(col = signif_col))
                } else if(qval[i, j] < 0.01) {
                    grid.text("**", x, y, gp = gpar(col = signif_col))
                } else if(qval[i, j] < 0.05) {
                    grid.text("*", x, y, gp = gpar(col = signif_col))
                }},
        ...)}

    else{
        ht <- Heatmap(
            nes,
            border = TRUE,
            col = col_fun,
            na_col = "black",
            border_gp = gpar(col = "black", lwd = 3),
            rect_gp = gpar(col = "white", lwd = 1),
            heatmap_legend_param = list(
                title = "NES",
                at = c(-scale, scale),
                legend_direction = "vertical",
                labels = c(paste0("-", scale), paste0(" ", scale)),
                width = unit(10, "mm")),
            show_heatmap_legend = T,
            show_column_names = T,
            column_names_side = "top",
            column_title_side = "top",,
            column_dend_reorder = T,
            column_dend_side = "bottom",
            show_row_names = T,
            row_title = "Pathways",
            row_dend_reorder = T,
            row_dend_side = "right",
            row_names_side = "left",
            row_names_max_width = unit(30, "cm"),
            column_title_gp = gpar(fontsize = 10,fontface="bold"),
            column_names_gp = gpar(fontsize = 10,fontface="bold"),
            row_title_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 8),
            cell_fun = function(j, i, x, y, w, h, fill) {
                if(qval[i, j] < 0.001) {
                    grid.text("***", x, y, gp = gpar(col = signif_col))
                } else if(qval[i, j] < 0.01) {
                    grid.text("**", x, y, gp = gpar(col = signif_col))
                } else if(qval[i, j] < 0.05) {
                    grid.text("*", x, y, gp = gpar(col = signif_col))
                }},
        ...)}
    return(ht)
}


# pathway to excel
pathways_to_excel <- function(x, diffexp = NULL, outdir){

    sheets <- list()
    collections <- unique(gsub(".*_", "", names(x)))

    if(file.exists(diffexp)){
        diffexp <- read.csv(diffexp, row.names = 1)
        if(!"gene" %in% colnames(diffexp)){
            diffexp$gene <- rownames(diffexp)}
        n <- unique(gsub(paste0("_", paste0(collections, collapse = "|_")), "", names(x)))
        if(!"cluster" %in% colnames(diffexp) & length(n) == 1){
            diffexp$cluster <- gsub(paste0("_", paste0(collections, collapse = "|_")), "", n)}}

    for(c in seq_along(collections)){
        selected.list <- x[which(str_detect(names(x), collections[c]))]
        for(i in seq_along(selected.list)){
            cluster.id <- gsub(paste0("_", paste0(collections, collapse = "|_")), "", names(selected.list)[i])

            diffexp_genes <- diffexp %>%
                filter(cluster == cluster.id & p_val_adj < 0.05) %>%
                .$gene

            edges <- str_split(selected.list[[i]]$core_enrichment, "\\/")

            for(j in seq_along(edges)){
                edges[[j]] <- paste0(intersect(edges[[j]], diffexp_genes), collapse = "/")}

            selected.list[[i]] <- selected.list[[i]]@result %>%
                mutate(
                    cluster = paste0(cluster.id),
                    signif_pathway = ifelse(qvalue < 0.05, "True", "False"),
                    signif_core_enrichment = unlist(edges))}
        sheets[[c]] <- bind_rows(selected.list)}
    names(sheets) <- collections
    writexl::write_xlsx(sheets, outdir)
    }



###
PlotAllHeatmap <- function(gsea, top_n = NULL, only_pos = T, only_signif = F, plot_qval = F, fill_na = F, keep_col_order = F, keep_row_order = F, signif_col = "white",...){
    
    if(is.data.frame(gsea)){
        stopifnot(c("cluster", "ID", "NES", "qvalue") %in% colnames(gsea))}

    colname_order <- NULL
    rowname_order <- NULL
    if(keep_col_order){
        colname_order <- unique(sort(gsea$cluster))}
    if(keep_row_order){
        rowname_order <- unique(sort(gsea$ID))}

    if(length(top_n) > 0){
        if(only_pos){
            plot_pathways <- gsea %>%
                group_by(cluster) %>%
                mutate(NES = as.numeric(NES)) %>%
                filter(NES > 0) %>%
		        arrange(desc(NES)) %>%
                slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
                .$ID}
        else{
            plot_pathways <- gsea %>%
                group_by(cluster) %>%
                mutate(NES = as.numeric(NES)) %>%
		        arrange(desc(NES^2)) %>%
                slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
                .$ID}
        gsea <- gsea %>%
            filter(ID %in% plot_pathways)}

    nes <- gsea %>%
        dplyr::select(cluster, ID, NES) %>%
        arrange(cluster) %>%
        mutate(cluster = as.factor(cluster), NES = NES) %>%
        complete(cluster, nesting(ID), fill = list(NES = NA)) %>%
        pivot_wider(names_from = cluster, values_from = NES) %>%
        column_to_rownames("ID")

    qval <- gsea %>%
        dplyr::select(cluster, ID, qvalue) %>%
        arrange(cluster) %>%
        mutate(cluster = as.factor(cluster), qvalue = qvalue) %>%
        complete(cluster, nesting(ID), fill = list(qvalue = 1)) %>%
        pivot_wider(names_from = cluster, values_from = qvalue) %>%
        column_to_rownames("ID")
    
    if(only_signif){
        keep_signif <- rownames(qval)[which(rowSums(qval < 0.05) > 0)]
        qval <- qval[keep_signif,]
        nes <- nes[keep_signif,]}

    if(fill_na){
        nes[is.na(nes)] <- 0}
    if(!plot_qval){
        qval[qval != 1] <- 1}

    min <- round(min(nes, na.rm = TRUE))
    max <- round(max(nes, na.rm = TRUE))
    scale = max(c(-min, max))
    col_fun = rev(brewer.pal(12,"RdBu"))

    ht <- Heatmap(
        nes,

        # columns
        show_column_names = T,
        cluster_columns = !keep_col_order,
        column_order = colname_order,
        column_names_side = "top",
        column_title_side = "top",
        column_dend_side = "bottom",
        column_title_gp = gpar(fontsize = 10,fontface="bold"),
        column_names_gp = gpar(fontsize = 10,fontface="bold"),

        # rows
        show_row_names = T,
        cluster_rows = !keep_row_order,
        row_order = rowname_order,
        row_names_side = "left",
        row_title_side = "left",
        row_dend_side = "right",
        row_title = "Pathways",
        row_names_max_width = unit(30, "cm"),
        row_title_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 8),

        # rest
        border = TRUE,
        na_col = "black",
        border_gp = gpar(col = "black", lwd = 3),
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(
            title = "NES",
            at = c(-scale, scale),
            legend_direction = "vertical",
            labels = c(paste0("-", scale), paste0(" ", scale)),
            width = unit(10, "mm")),
        col = col_fun,
        show_heatmap_legend = T,
        cell_fun = function(j, i, x, y, w, h, fill) {
            if(qval[i, j] < 0.001) {
                grid.text("***", x, y, gp = gpar(col = signif_col))
            } else if(qval[i, j] < 0.01) {
                grid.text("**", x, y, gp = gpar(col = signif_col))
            } else if(qval[i, j] < 0.05) {
            grid.text("*", x, y, gp = gpar(col = signif_col))
            }},
   ...)
   return(ht)
}


PlotBarPlot <- function(data, group1, group2, top_n = NULL, only.signif = T){

    stopifnot(all(c("NES", "qvalue", "ID") %in% colnames(data)))
    data$NES <- as.numeric(data$NES)
    data$qvalue <- as.numeric(data$qvalue)

    if(length(top_n) > 0){
        plot_pathways <- data %>%
                    arrange(desc(NES^2)) %>%
            mutate(direction = ifelse(NES > 0, "POS", "NEG")) %>%
            group_by(direction) %>%
            slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
            .$ID

        if(only.signif){
            plot_pathways <- data %>%
                filter(qvalue < 0.05) %>%
                arrange(desc(NES^2)) %>%
                mutate(direction = ifelse(NES > 0, "POS", "NEG")) %>%
                group_by(direction) %>%
                slice_min(n = top_n, order_by = qvalue, with_ties = F) %>%
                .$ID
        }

        data <- data %>%
            filter(ID %in% plot_pathways)}

    comparison <- paste0("Enriched in ", c(group2, group1))
    plot <- data %>%
        mutate(direction = factor(ifelse(NES < 0, group2, group1), c(group2, group1))) %>%
        mutate(psig = case_when(
            qvalue < 0.001 ~ "***",
            qvalue < 0.01 ~ "**",
            qvalue < 0.05 ~ "*",
            .default = "")) %>%
        ggplot(aes(x = fct_reorder(ID, NES), y = NES, fill = NES)) +
        geom_col(aes(stroke = psig), width = 0.85, col = "black") +
        geom_text(aes(label = psig, y = NES + 0.3 * sign(NES)), position = position_dodge(width = 0.95)) +
        xlab("Pathways") +
        scale_fill_distiller(palette = "RdBu") +
        guides(
            fill = guide_colorbar(
                title = "NES",
                title.position = "top",
                direction = "vertical",
                frame.colour = "black",
                ticks.colour = "black",
                order = 1)) +
        coord_flip() +
        theme_border() +
        scale_x_discrete(expand=c(0.05, 0.05)) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
    return(plot)
}