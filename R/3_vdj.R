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
    x@meta.data[c(colnames(vdj)[which(colnames(vdj) != "samples")], "cellxclone", "clone_by_count", "clone_by_percent")] <- NULL
    
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
        mutate(vdj_qc = ifelse(is.na(filter_contig), "Fail", "Pass")) %>%
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
        mutate(clone_by_percent = cellxclone/n(),
                clone_by_percent = case_when(
                    clone_by_percent <= 0.0001 ~ "Rare (X≤ 1e-04)",
                    clone_by_percent > 0.0001 & clone_by_percent <= 0.001 ~ "Small (1e-04< X ≤ 0.001)",
                    clone_by_percent > 0.001 & clone_by_percent <= 0.01 ~ "Medium (0.001< X ≤ 0.01)",
                    clone_by_percent > 0.01 & clone_by_percent <= 0.1 ~ "Large (0.01< X ≤ 0.1)",
                    clone_by_percent > 0.1 ~ "Hyperexpanded (0.1< X)")) %>%
        mutate(
            cellxclone = ifelse(clone_id == "No_contig", NA, cellxclone),
            clone_by_percent = ifelse(clone_id == "No_contig", NA, clone_by_percent),
            clone_by_count = ifelse(clone_id == "No_contig", NA, clone_by_count)) %>%
        mutate(clone_by_count = factor(clone_by_count, c("Single (0< X ≤1)", "Small (1< X ≤5)", "Medium (5< X ≤10)", "Large (11< X ≤30)", "Hyperexpanded (30< X)"))) %>%
        mutate(clone_by_percent = factor(clone_by_percent, c("Rare (X≤ 1e-04)", "Small (1e-04< X ≤ 0.001)", "Medium (0.001< X ≤ 0.01)", "Large (0.01< X ≤ 0.1)", "Hyperexpanded (0.1< X)"))) %>%
        column_to_rownames("Row.names") %>%
        as.data.frame(.)
    x@meta.data <- metadata[rownames(x@meta.data),]
    return(x)}

#' plot_vdj
#'
#' plot vdj information for cells
#' @param x Seurat object
#' @param variable metadata column to plot : "vdj_qc", "isotype", "clone_by_count", "clone_by_percent"
#' @param ... refer to plot_percent
#' @export
plot_vdj <- function(x, variable = "isotype", ...){

    stopifnot(variable %in% c("vdj_qc", "isotype") | str_detect(variable, "clone_by"))

    if(variable == "isotype"){
        cols <- c("brown", palette_list[["darjeeling_5"]][c(1,5,2,3)])
        names(cols) <- c("IgD", "IgM", "IgA", "IgG", "IgE")
        x@meta.data <- x@meta.data %>%
            filter(str_detect(isotype, "^Ig[DMAGE]$")) %>%
            mutate(isotype = factor(isotype, names(cols)))}

    else if(str_detect(variable, "clone_by")){
        cols <- palette_list[["zissou_5"]]
        names(cols) <- levels(x@meta.data[[variable]])}

    plot <- plot_percent(x, variable = variable, ...)

    return(plot)
}

#' plot_shm
#'
#' plot somatic hypermuation
#' @param x Seurat object
#' @param group.by column to group cells by
#' @param facet.by variable to facet by
#' @param cols colors
#' @param adjust adjust 
#' @export
plot_shm <- function(x, group.by, facet.by = NULL, cols = NULL, adjust = 1){
    #if(!all(unique(x@meta.data[[group.by]]) %in% names(cols))){
    #    stop('please make sure "cols" is a named vector of colors corresponding to the levels in "group.by"')}

    group <- c(group.by, facet.by)

    plot <- x@meta.data %>%
        filter(!is.na(mu_freq)) %>%
        ggplot(aes_string(x = group.by, y = "mu_freq")) +
        geom_violin(aes_string(fill = group.by), adjust = adjust, scale = "width") +
        geom_boxplot(width = 0.25, outliers = F) +
        guides(fill = guide_legend(title = "")) +
        theme_text() +
        xlab("") +
        ylab("Mutation Frequency")

    if(length(cols) > 0){
        plot <- plot +
            scale_fill_manual(values = cols)}

    if(length(facet.by) > 0){
        plot <- plot +
            theme_border() +
            facet_wrap(as.formula(paste0("~ ", facet.by)), nrow = 1) +
            facet_aes()}
    else{
        plot <- plot + theme_line()}

    return(plot)
}

#' plot_clone_number
#'
#' plot clone number
#' @param x Seurat object
#' @param group.by column to group cells by
#' @param cols colors
#' @export
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

#' match_cdr3
#'
#' match CDR3 amino acid sequence
#' @param query query dataframe with "junction_aa_VDJ" as column
#' @param reference reference dataframe with "junction_aa_VDJ" as column
#' @param method method to calculate CDR3 AA distance : "hamming" or "levenshtein"
#' @param match additional metadata columns to match
#' @param ncores no. of cores
#' @export
match_cdr3 <- function(query, reference, method = "hamming", match = c("v_call_VDJ_main", "j_call_VDJ_main"), ncores = 1){

  stopifnot(method %in% c("hamming", "levenshtein"))
  stopifnot(all(match %in% c("v_call_VDJ_main", "j_call_VDJ_main", "junction_aa_VDJ")))
  stopifnot(all(c(match, "junction_aa_VDJ") %in% colnames(query)))
  stopifnot(all(c(match, "junction_aa_VDJ") %in% colnames(reference)))

  query <- query %>%
    dplyr::select(c(match, "junction_aa_VDJ"))
  reference <- reference %>%
    distinct(!!!syms(c(match, "junction_aa_VDJ")))

  query$junction_aa_VDJ_length <- nchar(query$junction_aa_VDJ)
  reference$junction_aa_VDJ_length <- nchar(reference$junction_aa_VDJ)

  doParallel::registerDoParallel(ncores)

  output <- foreach::foreach(i = 1:nrow(reference), .combine=rbind) %dopar% {

    if(method == "hamming") {
      tmp <- query %>% semi_join(reference[i, ], by = c(match, "junction_aa_VDJ_length")) %>% rownames_to_column("cell_barcode")
      if(nrow(tmp) == 0) {
        return(data.frame())}
      cdr3 <- reference[i, ]$junction_aa_VDJ
      tmp$dist <- as.numeric(lapply(tmp$junction_aa_VDJ, function(x){alakazam::seqDist(cdr3, x, alakazam::getAAMatrix()) / nchar(cdr3)}))
      tmp$refseq <- cdr3
      return(tmp)}

    if(method == "levenshtein"){
      if(length(match) > 0) {
        tmp <- query %>% semi_join(reference[i, ], by = match) %>% rownames_to_column("cell_barcode")} 
      else{
        tmp <- query %>% rownames_to_column("cell_barcode")}
      if(nrow(tmp) == 0) {
        return(data.frame())}
      cdr3 <- reference[i, ]$junction_aa_VDJ
      tmp$dist <- as.numeric(lapply(tmp$junction_aa_VDJ, function(x){stringdist::stringdist(cdr3, x, method = 'lv') / max(nchar(cdr3), nchar(x))}))
      tmp$refseq <- cdr3
      return(tmp)}}

  if(nrow(output) > 0){
    output <- output %>% 
      group_by(cell_barcode) %>%
      slice_min(n = 1, order_by = dist, with_ties = F) %>%
      ungroup() %>%
      mutate(sim = ifelse(dist <= 1, 1-dist, 0)) %>%
      dplyr::select(-c("dist")) %>%
      arrange(desc(sim))
    return(output)}
  else{
    stop("No matching sequences")
  }
}

