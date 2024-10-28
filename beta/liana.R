suppressMessages({
library(tidyverse)
library(liana)
library(SingleCellExperiment, quietly = TRUE)
library(reticulate, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(liana, quietly = TRUE)
library(ExperimentHub, quietly = TRUE)
library(decoupleR, quietly = TRUE)
library(SingleCellExperiment)
library(furrr)

})


########## FIND SPECIFIC L-R INTERACTIONS ###########
find_ligand_receptors <- function(genes, resource = "Consensus"){
    if(!any(c("Consensus", "MouseConsensus") %in% resource)){
        stop("make sure to select a valid resource database")}

    db <- select_resource(resource)[[1]][,c(1,2,5,7)]
    genes_db <- db %>%
        filter(target_genesymbol %in% genes | source_genesymbol %in% genes) %>%
        distinct(source_genesymbol, target_genesymbol, .keep_all = T)
    return(genes_db)}

########### TENSOR C2C FACTORIZATION ###########
## convert seurat to sce object ##
convert_seurat_to_sce <- function(x, assay = "RNA"){
    x <- DietSeurat(x, assay = assay)
    x[[assay]] <- as(object = x[[assay]], Class = "Assay")
    sce <- as.SingleCellExperiment(x)
    return(sce)}

## liana context qc ##
liana_context_qc <- function(x, assay = "RNA", samples, celltypes, ncol = 6, min_cells = 3, min_samples = 3, min_prop = 0.05){
    if(class(x) == "Seurat"){
        x <- convert_seurat_to_sce(x, assay = assay)}
    if(class(x) == "SingleCellExperiment"){
        sce <- liana::filter_nonabundant_celltypes(
            x,
            sample_col = samples,
            idents_col = celltype,
            min_cells = min_cells,
            min_samples = min_samples,
            min_prop = min_prop)
        plot <- sce %>%
            get_abundance_summary(
                sample_col = batch,
                idents_col = celltype,
                min_cells = min_cells,
                min_samples = min_samples,
                min_prop = min_prop) %>%
            plot_abundance_summary(ncol = ncol)
        print(plot)
        return(sce)}}

## liana run tensor c2c
liana_context_c2c <- function(sce, resource = "Consensus", samples, celltypes, rank = NULL, outdir = NULL, conda_env = "/camp/home/hungm/.conda/envs/liana"){
    stopifnot(class(sce) == "SingleCellExperiment")
    sce <- liana::filter_nonabundant_celltypes(sce, sample_col = samples, idents_col = celltype)

    if(!length(rank) > 0){
        sce <- liana_bysample(
            sce = sce,
            sample_col = samples,
            idents_col = celltype,
            expr_prop = 0.1, # expression proportion threshold
            inplace=TRUE, # saves inplace to sce
            return_all = FALSE, # whether to return non-expressed interactions
            resource = resource)}

    sce <- liana_tensor_c2c(
        sce = sce,
        score_col = 'LRscore',
        rank = rank,  # set to None to estimate for you data!
        how='outer',  #  defines how the tensor is built
        conda_env = conda_env, # used to pass an existing conda env with cell2cell
        use_available = F, # detect & load cell2cell if available
        lr_sep = "->",
        device = "cpu")

    if(length(outdir) > 0){
        qsave(sce, outdir)}

    return(sce)}

## inspect rank error
plot_rank_error <- function(sce, rank = 10, max = 15){
    plot <- sce@metadata$tensor_res$elbow_metric_raw %>%
        t() %>%
        as.data.frame() %>%
        mutate(rank=row_number()) %>%
        pivot_longer(-rank, names_to = "run_no", values_to = "error") %>%
        group_by(rank) %>%
        summarize(average = mean(error),
                N = n(),
                SE.low = average - (sd(error)/sqrt(N)),
                SE.high = average + (sd(error)/sqrt(N))) %>%
        ggplot(aes(x=rank, y=average), group=1) +
        geom_ribbon(aes(ymin = SE.low, ymax = SE.high), alpha = 0.5) +
        geom_line(col='black', size = 1) +
        geom_point(fill = "white", size = 6, shape = 21) +
        geom_text(aes(label = rank), size = 3) +
        geom_vline(xintercept = rank, colour='red') + # rank of interest
        theme_border() +
        xlim(c(0, max)) +
        labs(y="Error", x="Rank")
    return(plot)}

## set interaction qc
set_interaction_qc <- function(sce, percentile = 0.1){
    keep_lr <- sce@metadata$tensor_res$interactions %>%
        pivot_longer(!lr, names_to = "context", values_to = "loading") %>%
        filter(quantile(loading, percentile) < loading) %>%
        .$lr
    sce@metadata$tensor_res$interactions <- sce@metadata$tensor_res$interactions[which(sce@metadata$tensor_res$interactions$lr %in% keep_lr), ]
    return(sce)}

## format progeny
format_progeny <- function(sce, organism = "human", percentile = 0){
    stopifnot(organism %in% c("human", "mouse"))
    sce <- set_interaction_qc(sce, percentile = percentile)
    progeny <- decoupleR::get_progeny(organism = organism, top=5000)
    progeny_lr <- generate_lr_geneset(sce, resource = progeny, lr_sep = "->")
    return(progeny_lr)}


## find progeny
find_lr_progeny <- function(sce, progeny_lr, sample, p.adj = T, percentile = 0){

  pcol <- ifelse(p.adj, "p_adj", "p_value")
  sce <- set_interaction_qc(sce, percentile = percentile)
  factors <- get_c2c_factors(sce, sample_col = sample)

  mat <- factors$interactions %>%
    column_to_rownames("lr") %>%
    as.matrix()

  cols <- c("black", "white")
  names(cols) <- c("p < 0.05", "ns")

  res <- decoupleR::run_ulm(
    mat = mat,
    network = progeny_lr,
    .source = "set",
    .target = "lr",
    minsize=5) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "fdr"),
      p_plot = -log10(.data[[pcol]]),
      condition = factor(gsub("\\.", " ", condition),  paste0("Factor ", 1:length(unique(condition)))),
      signif = ifelse(.data[[pcol]] <= 0.05, names(cols)[1], "ns"),
      signif = factor(signif, c(names(cols)[1], "ns")))

  psize <- quantile(as.numeric(res$p_plot), probs = c(0.33, 0.66, 1))
  psize <- unname(round(psize,1))

  res %>% # sig/isnig flag
    ggplot(aes(x = source, y = condition)) +
    geom_point(aes(fill = score, size = p_plot, color = signif), shape = 21, stroke = 1.5) +
    scale_fill_gradientn(colors = rev(brewer.pal(6, "RdBu"))) +
    scale_color_manual(values = cols) +
    scale_size_continuous(range = c(0, 6), breaks = psize) +
    theme_border() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x="Pathway", y="Factor", colour="Activity") +
    guides(
      fill = guide_colorbar(
        title = "Activity",
        order = 1,
        title.position = "top",
        direction = "horizontal",
        frame.colour = "black",
        ticks.colour = "black"),
      size = guide_legend(
        title = paste0("-log10(", pcol, ")"),
        order = 2,
        title.position = "top",
        direction = "horizontal",
        override.aes = list(fill = "black"),
        theme = theme(legend.text=element_text(size=8))),
      color = guide_legend(title = "Significance", order = 3)) +
    theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
    )
}


