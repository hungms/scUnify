#library(alakazam)
#library(scoper)
#library(dplyr)
#library(Seurat)
#library(shazam)
#library(dowser)

#' seurat_add_dandelion
#'
#' add dandelion output to seurat object metadata
#' @param x Seurat object
#' @param vdj dataframe of "all/filtered_contig_dandelion.tsv"
#' @param paired if TRUE, keep cells with a single heavy and a single light chain sequence
#' @export
seurat_add_dandelion <- function(x, vdj, paired = T){
    vdj <- vdj %>%
        filter(filter_contig == "False") %>%
        dplyr::select(!samples)

    if(paired){
        vdj <- vdj %>%
            filter(chain_status == "Single pair" & productive_VDJ == "T")  %>%
            as.data.frame(.)}

    stopifnot(all(rownames(vdj) %in% rownames(x@meta.data)))

    metadata <- x@meta.data %>%
        merge(., vdj, by = 0, all.x = T) %>%
        mutate(vdj_qc = ifelse(filter_contig == "False", "Pass", "Fail")) %>%
        group_by(clone_id) %>%
        mutate(cellxclone = n()) %>%
        ungroup() %>%
        mutate(
            clone_by_count = case_when(
                cellxclone > 30 ~ "Hyperexpanded (30< X)",
                cellxclone > 10 & cellxclone <= 30 ~ "Large (11< X ≤30)",
                cellxclone > 5 & cellxclone <= 10 ~ "Medium (5< X ≤10)",
                cellxclone >= 2 & cellxclone <=5 ~ "Small (1< X ≤5)",
                .default = "Single (0< X ≤1)")) %>%
        mutate(clone_by_count = factor(clone_by_count, c("Single (0< X ≤1)", "Small (1< X ≤5)", "Medium (5< X ≤10)", "Large (11< X ≤30)", "Hyperexpanded (30< X)"))) %>%
        mutate(clone_by_percent = cellxclone/n(),
                clone_by_percent = case_when(
                    clone_by_percent <= 0.0001 ~ "Rare (X≤ 1e-04)",
                    clone_by_percent > 0.0001 & clone_by_percent <= 0.001 ~ "Small (1e-04< X ≤ 0.001)",
                    clone_by_percent > 0.001 & clone_by_percent <= 0.01 ~ "Medium (0.001< X ≤ 0.01)",
                    clone_by_percent > 0.01 & clone_by_percent <= 0.1 ~ "Large (0.01< X ≤ 0.1)",
                    clone_by_percent > 0.1 ~ "Hyperexpanded (0.1< X)")) %>%
        mutate(clone_by_percent = factor(clone_by_percent, c("Rare (X≤ 1e-04)", "Small (1e-04< X ≤ 0.001)", "Medium (0.001< X ≤ 0.01)", "Large (0.01< X ≤ 0.1)", "Hyperexpanded (0.1< X)"))) %>%
        mutate(
            cellxclone = ifelse(clone_id == "No_contig", NA, cellxclone),
            clone_by_percent = ifelse(clone_id == "No_contig", NA, clone_by_percent),
            clone_by_count = ifelse(clone_id == "No_contig", NA, clone_by_count)) %>%
        column_to_rownames("Row.names") %>%
        as.data.frame(.)
    x@meta.data <- metadata[rownames(x@meta.data),]
    return(x)}

#' plot_vdj_qc
#'
#' plot proportion of cells that has VDJ library
#' @param x Seurat object
#' @param group.by column to group cells by
#' @export
plot_vdj_qc <- function(x, group.by = "samples"){
    plot <- x@meta.data %>%
        group_by_at(c(group.by, "vdj_qc")) %>%
        summarize(count = n()) %>%
        group_by_at(group.by) %>%
        mutate(pct = count*100/sum(count)) %>%
        ggplot(aes_string(x = group.by, y = "pct", fill = "vdj_qc")) +
        geom_col(width = 0.85, position = "stack", col = "white") +
        guides(fill = guide_legend(title = "")) +
        theme_line() +
        xlab("") +
        ylab("Proportions (%)")
    return(plot)}

#' plot_vdj
#'
#' plot vdj information for cells
#' @param x Seurat object
#' @param group.by column to group cells by
#' @export
plot_vdj <- function(x, group.by, facet.by = NULL, variable = "isotype"){
    
    if(length(facet.by) > 0){
        group <- c(group.by, facet.by)}
    else{
	group <- group.by}

    if(variable == "isotype"){
        cols <- c("brown", palette_list[["darjeeling_5"]][c(1,5,2,3)])
        names(cols) <- c("IgD", "IgM", "IgA", "IgG", "IgE")
        metadata <- x@meta.data %>%
            filter(str_detect(isotype, "^Ig[DMAGE]$")) %>%
            mutate(isotype = factor(isotype, names(cols)))}
    else if(str_detect(variable, "clone_by")){
        cols <- palette_list[["zissou_5"]]
        names(cols) <- levels(x@meta.data[[variable]])
        metadata <- x@meta.data}

    plot <- metadata %>%
        group_by_at(c(group, variable)) %>%
        summarize(count = n()) %>%
        group_by_at(group) %>%
        mutate(pct = count*100/sum(count)) %>%
        ggplot(aes_string(x = group.by, y = "pct", fill = variable)) +
        geom_col(width = 0.85, position = "stack", col = "white") +
        scale_fill_manual(values = cols) +
        guides(fill = guide_legend(title = "")) +
        theme_line() +
        theme_text() +
        xlab("") +
        ylab("Proportions (%)")

    if(length(facet.by) == 1){
        plot <- plot +
            facet_wrap(as.formula(paste0("~ ", facet.by)), nrow = 1) +
            facet_aes()}

    return(plot)
}

#' plot_shm
#'
#' plot somatic hypermuation
#' @param x Seurat object
#' @param group.by column to group cells by
#' @param facet.by variable to facet by
#' @param cols colors
#' @export
plot_shm <- function(x, group.by, facet.by = NULL, cols = NULL){
    #if(!all(unique(x@meta.data[[group.by]]) %in% names(cols))){
    #    stop('please make sure "cols" is a named vector of colors corresponding to the levels in "group.by"')}

    if(length(facet.by) > 0){
        group <- c(group.by, facet.by)}
    else{
        group <- group.by}

    plot <- x@meta.data %>%
        filter(!is.na(mu_freq)) %>%
        ggplot(aes_string(x = group.by, y = "mu_freq", fill = group.by)) +
        geom_violin(trim = T, adjust = 1.5, drop = T, bw = "nrd0", scale = "width") +
        guides(fill = guide_legend(title = "")) +
        theme_line() +
        xlab("") +
        ylab("Mutation Frequency")

    if(length(cols) > 0){
        plot <- plot +
            scale_fill_manual(values = cols)}

    if(length(facet.by) > 0){
        plot <- plot +
            facet_wrap(as.formula(paste0("~ ", facet.by)), nrow = 1) +
            facet_aes()}

    return(plot)
}



plot_clone_number <- function(x, group.by, cols = NULL){
    plot <- x@meta.data %>%
        group_by_at(group.by) %>%
        summarize(`clones_per_cell` = length(unique(clone_id))*100/n()) %>%
        ggplot(aes_string(x = group.by, y = "clones_per_cell", fill = group.by)) +
        geom_col(width = 0.85, color = "white") +
        theme_line() +
        theme_text() +
        xlab(NULL) +
        ylab("Average Clone Number")
    if(length(cols) > 0){
        plot <- plot + scale_fill_manual(values = cols)}
    return(plot)
}