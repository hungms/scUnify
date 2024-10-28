plot_liana_c2c <- function(sce, sample, celltype, group, group_cols = NULL, celltype_cols = NULL){

    factors <- get_c2c_factors(sce, sample_col = sample, group_col = group)
    colnames(factors$contexts) <- gsub("\\.", " ", colnames(factors$contexts))
    colnames(factors$interactions) <- gsub("\\.", " ", colnames(factors$interactions))
    colnames(factors$senders) <- gsub("\\.", " ", colnames(factors$senders))
    colnames(factors$receivers) <- gsub("\\.", " ", colnames(factors$receivers))

    levels <- levels(sce[[group]])
    factors$contexts[[group]] <- factor(factors$contexts[[group]], levels)
    levels <- levels(sce[[sample]])

    # contexts
    contexts <- factors$contexts %>%
        pivot_longer(cols = -c("context", group), names_to = "factor", values_to = "loadings") %>%
        mutate(context = factor(context, levels)) %>%
        ggplot(aes(x=context, y=loadings, fill=.data[[group]])) +
            geom_col(col = "white", width = 0.85) +
            facet_grid(factor ~ .) +
            theme_bw(base_size = 14) +
            theme(
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                strip.text.y = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
            theme_border() +
            facet_aes() +
            ggtitle('Contexts') +
            ylab(NULL)

    # lr
    lr <- factors$interactions %>%
        pivot_longer(-lr, names_to = "factor", values_to = "loadings") %>%
        ggplot(aes(x=lr, y=loadings)) +
            geom_col(stat="identity", width = 0.85) +
            facet_grid(factor ~ .) +
            theme_bw(base_size = 14) +
            theme(
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                strip.background = element_blank(),
                strip.text.y = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
            ggtitle('Interactions') +
            ylab(NULL) +
            theme(
                panel.grid = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1))

    # Sender cells
    senders <- factors$senders %>%
        pivot_longer(cols = -celltype, names_to = "factor", values_to = "loadings") %>%
        ggplot(aes(x=celltype, y=loadings, fill=celltype)) +
            geom_col(stat="identity", col = "white", width = 0.85) +
            facet_grid(factor ~ .) +
            theme_bw(base_size = 14) +
            theme(
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                strip.background = element_blank(),
                strip.text.y = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
            ylab(NULL) +
            theme_border() +
            facet_aes() +
            ggtitle('Senders')

    # Receiver cells
    receivers <- factors$receivers %>%
        pivot_longer(cols = -celltype, names_to = "factor", values_to = "loadings") %>%
        ggplot(aes(x=celltype, y=loadings, fill=celltype)) +
            geom_col(stat="identity", col = "white", width = 0.85) +
            facet_grid(factor ~ .) +
            theme_bw(base_size = 14) +
            theme(
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                strip.background = element_blank(),
                axis.ticks.x=element_blank(),
                strip.text.y = element_text(size=15, face = "bold"),
                plot.title = element_text(hjust = 0.5)) +
            ylab(NULL) +
            theme_border() +
            facet_aes() +
            ggtitle('Receivers')

    if(length(group_cols) > 0){
        contexts <- contexts +
            scale_fill_manual(values = group_cols)}

    if(length(celltype_cols) > 0){
        senders <- senders +
            scale_fill_manual(values = celltype_cols)
        receivers <- receivers +
            scale_fill_manual(values = celltype_cols)}

    # Assemble overview plot
    overview <- patchwork::wrap_plots(list(contexts, lr, senders, receivers), ncol=4, nrow(1)) + patchwork::plot_layout(guides = "collect")
    grid::grid.draw(overview)}

## boxplot
plot_context_boxplot <- function(sce, sample, group, method="t.test", cols = NULL, nrow = 1, ...){
    factors <- get_c2c_factors(sce, sample_col = sample, group_col = group)
    all <- factors$contexts %>%
        .[[group]] %>%
        unique(.) %>%
        as.character(.) %>%
        combn(., 2, simplify = FALSE)
    levels <- levels(sce[[group]])
    factors$contexts[[group]] <- factor(factors$contexts[[group]], levels)
    plot <- factors$contexts %>%
        dplyr::select(!!all_of(group), starts_with("Factor")) %>%
        pivot_longer(-!!all_of(group), names_to = "fact", values_to = "loadings") %>%
        mutate(
            fact = factor(gsub("\\.", " ", fact), paste0("Factor ", 1:length(unique(fact))))) %>%
        ggplot(aes(x=.data[[group]], y=loadings, fill=.data[[group]])) +
            geom_boxplot(width = 0.5, outlier = F) +
            geom_jitter(aes(fill = .data[[group]]), color = "black", width = 0.1, size  = 3, shape = 21, stroke = 1.1) +
            theme_border() +
            facet_wrap(~fact, nrow = nrow) +
            facet_aes() +
            guides(fill = guide_legend(title = "")) +
            xlab("") +
            ylab("Loadings") +
            ggtitle("") +
            scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
            ggpubr::stat_compare_means(method = method, comparisons = all, ...)

    if(length(cols) > 0){
        plot <- plot + scale_fill_manual(values = cols)}
    return(plot)
}

##
plot_context_heat <- function(sce, sample, group, cols, ...){

    factors <- get_c2c_factors(sce, sample_col = sample, group_col = group)
    colnames(factors$contexts) <- gsub("\\.", " ", colnames(factors$contexts))
    levels <- levels(sce[[group]])
    factors$contexts[[group]] <- factor(factors$contexts[[group]], levels)

    # Samples dictionary
    meta_dict <- factors$contexts %>%
        dplyr::select(context, !!group) %>%
        deframe()

    if(length(cols) > 0){
        names(cols) <- unique(meta_dict)
        col_list <- list(` ` = cols)}
    else{
        col_list <- list()}

    contexts_mat <- factors$contexts %>%
        column_to_rownames("context") %>%
        dplyr::select(starts_with("Factor")) %>%
        t()

    col_fun = rev(brewer.pal(6,"RdBu"))
    contexts_mat %>%
        ComplexHeatmap::Heatmap(
                top_annotation = ComplexHeatmap::HeatmapAnnotation(
                ` ` = meta_dict,
                col = col_list,
                show_annotation_name = FALSE,
                border = T,
                gp = gpar(col = "black", lwd = 1),
                simple_anno_size = grid::unit(0.3, "cm")
                ),
            col = col_fun,
            border = T,
            border_gp = gpar(col = "black", lwd = 3),
	    rect_gp = gpar(col = "white", lwd = 1),
            name = "Context\nLoadings",
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            ...) %>%
            draw(., merge_legend = TRUE)

}

##
plot_interaction_heat <- function(sce, lr_sep=" -> ", n = 5, percentile = 0, cluster_columns=FALSE, ...){

    sce <- set_interaction_qc(sce, percentile = percentile)
    factors <- get_c2c_factors(sce)
    colnames(factors$interactions) <- gsub("\\.", " ", colnames(factors$interactions))

    # Top n Interactions per factor heatmap
    top_lrs <- factors$interactions %>%
        pivot_longer(-lr, names_to = "fact", values_to = "loadings") %>%
        group_by(fact) %>%
        top_n(wt=loadings, n=n) %>%
        pull(lr)

    lrs_mat <- factors$interactions %>%
        filter(lr %in% top_lrs) %>%
        mutate(lr = gsub(as.character(str_glue("\\{lr_sep}")), " -> ", lr)) %>%
        as.data.frame() %>%
        column_to_rownames("lr") %>%
        as.matrix()

    col_fun = rev(brewer.pal(6,"RdBu"))
    lrs_mat %>%
        ComplexHeatmap::Heatmap(
            name = "Interaction\nLoadings",
            col = col_fun,
            border = T,
            rect_gp = gpar(col = "white", lwd = 1),
            border_gp = gpar(col = "black", lwd = 3),
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            cluster_columns=cluster_columns,
            ...)
}

## plot lr density
plot_lr_density <- function(sce, sample, percentile = 0.1){
    factors <- get_c2c_factors(sce, sample_col = sample)
    threshold <- factors$interactions %>%
        pivot_longer(!lr, names_to = "context", values_to = "loading") %>%
        .$loading %>%
        quantile(., percentile)

    factors$interactions %>%
        pivot_longer(!lr, names_to = "context", values_to = "loading") %>%
        ggplot(aes(x = loading)) +
        geom_density(size = 1, fill = "grey") +
        theme_border() +
        geom_vline(xintercept = threshold, color = "red") +
        xlab("Loading") +
        ylab("Density")
        }

## plot progeny weight
plot_progeny_weight <- function(sce, sample, factor, percentile){

  stopifnot(is.numeric(factor))
  sce <- set_interaction_qc(sce, percentile = percentile)
  factors <- get_c2c_factors(sce = sce, sample_col = sample)
  cols <- brewer.pal(6, "RdBu")[c(1,6)]

  # Plot LRs associated with Estrogen
  lrs <-  factors$interactions %>%
    left_join(progeny_lr, by="lr") %>%
    filter(!is.na(set)) %>%
    mutate(loading = .data[[paste0("Factor.", factor)]]) %>%
    dplyr::select(lr, set, mor, loading) %>%
    mutate(lr = gsub(as.character(str_glue("\\^")), " -> ", lr)) %>%
    mutate(
      weight = ifelse(mor >= 0, "Positive", "Negative"),
      weight = factor(weight, c("Positive", "Negative")))
  plot <- lrs %>%
    # only label those that are > x
    mutate(lr = if_else(loading>=0.001 & abs(mor) > 2, lr, "")) %>%
    ggplot(aes(x=mor, y=loading)) +
    stat_smooth(method = "lm", col = "black") +
    geom_point(aes(fill = weight), alpha = 0.8, size = 3, shape = 21) +
    ggrepel::geom_label_repel(aes(label = lr, fill = weight), col = "white", segment.color="black", size = 3) +
    guides(fill = guide_legend(title = "Weight", override.aes = list(label = ""))) +
    scale_fill_manual(values = cols) +
    facet_wrap(~set, ncol = 4, scales = "free") +
    xlab("Mode of Regulation (MoR)") +
    ylab("Loadings") +
    theme_border() +
    facet_aes()
    labs(x="Pathway Weight", y="LR Loading")

  return(plot)
  }

## plot lr dotplot
plot_lr_dotplot  <- function(
    sce, factor, percentile = 0.1,
    x, group.by.x, split.by.x = NULL, diffexp.x = NULL, title.x = "Ligands",
    y, group.by.y, split.by.y = NULL, diffexp.y = NULL, title.y = "Receptors",
    assay = "RNA", slot = "data", scale = T, top_n = 5, plot = T,
    rel_widths = c(1.5, 1.35, 1.2)){

    interaction <- sce@metadata$tensor_res$interactions %>%
        mutate(factor = .data[[paste0("Factor.", factor)]]) %>%
        filter(quantile(factor, percentile) < factor) %>%
        slice_max(n = top_n, order_by = factor, with_ties = F) %>%
        arrange(factor) %>%
        dplyr::select(lr)
    source <- gsub("->.*", "", interaction$lr)
    target <- gsub(".*->", "", interaction$lr)
    interaction <- data.frame(source, target) %>%
        mutate(order = factor(1:nrow(.))) %>%
        separate(col = source, into = c("before", "after"), sep = "_", remove = FALSE) %>%
        dplyr::select(!source) %>%
        pivot_longer(cols = c("before", "after"), names_to = "temp", values_to = "source") %>%
        separate(col = target, into = c("before", "after"), sep = "_", remove = FALSE) %>%
        dplyr::select(!c(target, temp)) %>%
        pivot_longer(cols = c("before", "after"), names_to = "temp", values_to = "target") %>%
        dplyr::select(!temp) %>%
        filter(!is.na(source)) %>%
        filter(!is.na(target)) %>%
        mutate(
            source.diffexp = ifelse(source %in% diffexp.x$gene, T, F),
            target.diffexp = ifelse(target %in% diffexp.y$gene, T, F)) %>%
        filter(source.diffexp == T & target.diffexp == T)

    if(plot == F){
        return(interaction)}

    plot.x <- scDotPlot(
        x,
        features = interaction[["source"]],
        group.by = group.by.x,
        split.by = split.by.x,
        diffexp = diffexp.x,
        assay = assay,
        slot = slot,
        scale = scale,
        palette = "Reds",
        direction = 1,
        coord_flip = T) +
        ggtitle(title.x) +
        theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    plot.x <- plot.x +
        scale_y_discrete(
            position = "left", ,
            breaks = plot.x$data$order,
            labels = plot.x$data$Gene)
    plot.y <- scDotPlot(
        y,
        features = interaction[["target"]],
        group.by = group.by.y,
        split.by = split.by.y,
        diffexp = diffexp.y,
        assay = assay,
        slot = slot,
        scale = scale,
        palette = "Blues",
        direction = 1,
        coord_flip = T) +
        ggtitle(title.y) +
        theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            strip.text.y = element_blank())

    legend.x <- get_legend(plot.x)
    legend.y <- get_legend(plot.y)
    legend <- plot_grid(NULL, legend.x, NULL, legend.y, ncol = 1, row = 4)

    plot.x  <- plot.x + theme(legend.position="none")
    plot.y  <- plot.y + theme(legend.position="none")

    plot_grid(plotlist = list(plot.x, plot.y, legend), ncol = 3, rel_widths = rel_widths, align = "hv", axis = "tb")}



## under development
format_msigdb <- function(sce, organism = "human"){
    msigdb_lr <- get_resource("MSigDB", organism = organism) %>%
        mutate(weight = 1) %>%
        mutate(collection = case_when(
            collection %in% "hallmark" ~ "HM",
            collection %in% c("go_biological_process", "go_cellular_component", "go_molecular_function") ~ "GO",
            collection %in% c("biocarta_pathways") ~ "BIOCARTA",
            collection %in% c("kegg_pathways") ~ "KEGG",
            collection %in% c("reactome_pathways") ~ "REACTOME",
            .default = "remove")) %>%
        filter(collection != "remove") %>%
        group_by(collection) %>%
        mutate(count = n()) %>%
        filter(count == 1) %>%
        ungroup() %>%
        dplyr::select(genesymbol, geneset, weight, collection)
    colnames(msigdb_lr) <- c("lr", "set", "mor", "collection")
    return(msigdb_lr)}