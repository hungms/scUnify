#' theme_border
#'
#' ggplot2 aesthetic option with borders
#' @export
theme_border <- function(){
    theme(
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
        }

#' theme_line
#'
#' ggplot2 aesthetic option with borders
#' @export
theme_line <- function(){
    list(
        scale_y_continuous(expand = c(0, 0)),
        theme(
            axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", size = 1)))
        }

#' facet_aes
#'
#' ggplot2 aesthetic option with facet_wraps
#' @export
facet_aes <- function(){
    theme(
	strip.background = element_blank(),
        strip.text = element_text(face="bold", size=12))}

#' umap_aes
#'
#' ggplot2 aesthetic option with UMAPs
#' @export
umap_aes <- function(){
    list(
	xlab("UMAP1"),
        ylab("UMAP2"),
        theme(
            panel.background = element_rect(fill = "white"),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            aspect.ratio = 1))}

#' scUMAP
#'
#' Plot UMAP from Seurat object
#' @param x Seurat object
#' @param reduction reduction name, defaults to "umap"
#' @param group.by metadata column to group by 
#' @param cols a vector of colors, defaults to NULL
#' @param count If TRUE, count the no. of cells for each group. Defaults to FALSE
#' @param pt.size size of each cell/point on UMAP
#' @param legend.ncol no. of columns for group.by keys in legend
#' @param ... arguments for Seurat::DimPlot()
#' @export
scUMAP <- function(x, reduction = "umap", group.by, cols = NULL, count = F, pt.size = NULL, legend.ncol = 1, ...){
    if(count){
        count.labels <- table(x@meta.data[[group.by]])
        count.labels <- as.factor(paste0(names(count.labels), " (", count.labels, ")"))
        count.labels <- unname(count.labels)}
    else{
        if(!is.factor((x@meta.data[[group.by]]))){
            x@meta.data[[group.by]] <- factor(x@meta.data[[group.by]], sort(unique(x@meta.data[[group.by]])))}
        count.labels <- levels(x@meta.data[[group.by]])}

    plotlist <- DimPlot(x, reduction = reduction, group.by = group.by, cols = cols, pt.size = pt.size, combine = F, ...)
    
    for(i in seq_along(plotlist)){
        plotlist[[i]] <- plotlist[[i]] +
            scale_color_discrete(labels = count.labels) +
            guides(color = guide_legend(title = "", title.position = "top", ncol = legend.ncol, override.aes = list(size=3))) +
	        ggtitle("") +
            umap_aes()}

    if(length(plotlist) == 1){
        plotlist <- plotlist[[1]]
        if(length(cols) > 0){
            plotlist <- plotlist + scale_color_manual(values = cols, labels = count.labels)}}

    return(plotlist)}

#' scDotPlot
#'
#' Plot DotPlot from Seurat object
#' @param x Seurat object
#' @param group.by metadata column to group by 
#' @param split.by metadata column to group by 
#' @param features feature names stored in a vector or in a named list
#' @param assay assay name, defaults to "RNA"
#' @param slot slot name, defaults to "data"
#' @param scale if TRUE, scale gene expression across cells. Defaults to TRUE
#' @param palette a vector of colors, defaults to "Reds" from RColorBrewer
#' @param direction direction of palette, default to 1
#' @param diffexp whether to plot significance, input can be either a dataframe or a csv file path of the Seurat::FindAllMarkers output. Defaults to NULL.
#' @param coord_flip If TRUE, plot features on y-axis. Defaults to TRUE
#' @param unique If TRUE, do not repeat feature names. Defaults to FALSE
#' @export
scDotPlot <- function(
    x, group.by, split.by = NULL, features,
    assay = "RNA", slot = "data", scale = T,
    palette = "Reds", direction = 1, diffexp = NULL,
    coord_flip = T, unique = F){

    # convert genelists to vector, and filter unique genes if necessary
    message("Step 1 : Formatting gene list")
    features_df <- format_features(features = features, unique = unique)
    var <- features_df$Gene

    # get and format gene expression
    message("Step 2 : Summarizing gene expression")
    expdf <- summarize_exp(
        x, features = var, group.by = group.by,
        split.by = split.by, assay = assay, slot = slot, scale = scale)

    # add signif
    message("Step 3 : Annotating expression significance")
    if(length(diffexp) > 0){
        expdf <- annotate_diffexp(expdf = expdf, diffexp = diffexp, split.by = split.by)
        stroke_legend <- TRUE}

    # construct final df
    message("Step 4 : Construct final expression dataframe")
    expdf <- features_df %>%
        merge(., expdf, by = "Gene", all.x = T) %>%
        filter(Gene %in% rownames(x[[assay]])) %>%
        arrange(order)

    # facet and flip if necessary
    message("Step 5 : Plotting gene expression")

    # add palette
    palette <- brewer.pal(9, palette)
    if(direction < 0){
        palette <- rev(palette)}

    if(coord_flip){
        columns <- c("Group", "order", "Split")
        plot <- basicDotPlot(
            df = expdf, x=columns[1], y = columns[2],
            size = "Pct", fill = "Avg",
            stroke_legend = stroke_legend, palette = palette)
        plot <- plot +
            ggh4x::facet_grid2(vars(GeneSet), vars(Split), scales = "free", independent = "none", space = "free", switch = "y") +
            scale_x_discrete(position = "top") +
            scale_y_discrete(
                position = "right",
                breaks = expdf %>% distinct(order, Gene) %>% .$order,
                labels = expdf %>% distinct(order, Gene) %>% .$Gene) +
            theme(
                axis.text.x = element_text(hjust = 0),
                axis.text.x.top = element_text(vjust = 0.5),
                strip.placement = "outside")}
    else{
        columns <- c("order", "Group", "Split")
        plot <- basicDotPlot(
            df = expdf, x=columns[1], y = columns[2],
            size = "Pct", fill = "Avg",
            stroke_legend = stroke_legend, palette = palette)
        plot <- plot +
            ggh4x::facet_grid2(vars(Split), vars(GeneSet), scales = "free", independent = "none", space = "free", switch = "y") +
            scale_x_discrete(
                breaks = expdf %>% distinct(order, Gene) %>% .$order,
                labels = expdf %>% distinct(order, Gene) %>% .$Gene) +
            theme(strip.placement = "outside")}

    return(plot)
}

basicDotPlot <- function(df, x, y, fill, size, stroke_legend = F, palette){
    if("signif" %in% colnames(df)){
        stroke <- "signif"
        stroke_legend <- guide_legend(
            direction = "vertical",
            order = 3,
            override.aes = list(fill = palette[6], size = 6),
            theme = theme(legend.text=element_text(size=10)),
            title = "Significance",
            title.position = "top")}
    else{
	stroke <- NULL
        stroke_legend <- guide_none()}

    # basic dotplot
    plot <- ggplot(df, aes_string(x = x, y = y)) +
        geom_point(aes_string(size = size, fill = fill, stroke = stroke), shape=21) +
        scale_discrete_manual(aesthetics = "stroke", values = c(1.3, 0)) +
        scale_color_manual(values = c("black", "white")) +
        scale_fill_gradientn(colors = palette) +
        scale_size_continuous("% detected", range = c(0,8), breaks = c(0, 33, 66), labels = c("0%", "33%", "66%"))  +
        ylab("") +
        xlab("") +
        guides(
            fill = guide_colorbar(
                title = "Average Expression",
                title.position = "top",
                direction = "horizontal",
                frame.colour = "black",
                ticks.colour = "black",
                order = 1),
            size = guide_legend(
                direction = "horizontal",
                order = 2,
                override.aes = list(fill = "black"),
                theme = theme(legend.text=element_text(size=10)),
                title = "Percent Expressed",
                title.position = "top"),
            stroke = stroke_legend) +
        theme_border() +
        facet_aes() +
        theme(
            axis.text.x = element_text(size=14, color="black"),
            axis.text.y = element_text(size=12, color="black"))
    return(plot)
}

format_features <- function(features, unique = F){
    if(is.character(features)){
        features <- list(` ` = features)}

    if(is.list(features)){
        features_df <- data.frame(
            GeneSet = rep(names(features), sapply(features, length)),
            Gene = unlist(features, use.names = FALSE)) %>%
            mutate(GeneSet = factor(GeneSet, unique(names(features)))) %>%
            group_by(GeneSet) %>%
            mutate(Gene = factor(Gene, rev(unique(Gene)))) %>%
            arrange(Gene) %>%
            ungroup() %>%
            mutate(order = factor(1:nrow(.)))

        if(unique){
            features_df <- features_df %>% distinct(GeneSet, Gene, .keep_all = T)}}
    return(features_df)}

summarize_exp <- function(x, features, group.by, split.by = NULL, assay = "RNA", slot = "data", scale = T){

    # extract expression
    DefaultAssay(x) <- assay
    df <- FetchData(x, vars = c(features), assay = assay, slot = slot)
    features <- features[which(features %in% colnames(df))]
    if(scale){
        df <- as.data.frame(scale(as.matrix(df)))}

    # define groups
    df[["Group"]] <- x@meta.data[[group.by]]
    if(length(split.by) > 0){
        df[["Split"]] <- x@meta.data[[split.by]]}
    else{
        df[["Split"]] <- split.by <- ""}

    # summarize expression
    group_cols <- c("Group", "Gene", "Split")
    expdf <- df %>%
        pivot_longer(cols = unname(features), names_to = "Gene", values_to = "Expression") %>%
        group_by_at(group_cols) %>%
        summarise(
            Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100) %>%
        ungroup() %>%
        arrange(Gene)
    expdf[is.na(expdf)] <- 0
    return(expdf)}

# add differential expression significance
annotate_diffexp <- function(expdf, diffexp, split.by = NULL){

    if(!is.data.frame(diffexp)){
        stopifnot(file.exists(diffexp))
        diffexp <- read.csv(diffexp, row.names = 1)}

    diffexp <- diffexp %>%
        mutate(
            signif = case_when(
                avg_log2FC > 0 & p_val_adj < 0.05 & pct.1 > 0.1 ~ "p < 0.05",
                .default = "ns"),
            signif = factor(signif, c("p < 0.05", "ns")),
            Gene = gene,
            Group = cluster)

    if(length(split.by) > 0){
        group_cols <- c("Group", "Gene", "Split")
        diffexp$Split <- diffexp[[split.by]]}
    else{
        group_cols <- c("Group", "Gene")}

    expdf <- expdf %>%
        #filter(Gene %in% unique(expdf$Gene)) %>%
        #filter(Group %in% unique(expdf$Group)) %>%
        merge(., diffexp, by = group_cols, all.x = T)  %>%
        dplyr::select(!c(gene, cluster))
    expdf[is.na(expdf)] <- "ns"
    
    return(expdf)
    }

#' scFeaturePlot
#'
#' Plot FeaturePlot from Seurat object
#' @param x Seurat object
#' @param features feature names stored in a vector
#' @param ncol no. of columns for the plots
#' @param plot If TRUE, plot with cowplot::plot_grid(). Defaults to FALSE to return a list of ggplot objects
#' @param ... arguments for Seurat::FeaturePlot()
#' @export
scFeaturePlot <- function(x, features, ncol, plot = F, ...){
    suppressMessages({
        plotlist <- FeaturePlot(x, features = features, ncol = ncol, combine = F, ...)
        for(x in seq_along(plotlist)){
            plotlist[[x]] <- plotlist[[x]] + 
                scale_color_viridis() +
                guides(color = guide_colorbar(title = "")) +
		umap_aes()}})
    p <- plotlist
    if(plot){
        p <- plot_grid(plotlist = plotlist, ncol = ncol, align = "hv")}
    return(p)}


#' scDensityUMAP
#'
#' Plot DensityUMAP from Seurat object
#' @param x Seurat object
#' @param reduction reduction name, defaults to "umap"
#' @param split.by metadata column to group by
#' @param adjust adjust from stat_density_2d()
#' @export
scDensityUMAP <- function(x, reduction = "umap", split.by, adjust = 1){
    umap <- as.data.frame(Embeddings(x@reductions[[reduction]]))
    umap2 <- umap
    umap$split <- x@meta.data[[split.by]]
    ncol <- length(unique(umap$split))
    colnames <- colnames(umap)

    plot <- ggplot(umap, aes_string(x = colnames[1], y = colnames[2])) +
    	geom_point(data = umap2, color = "grey80", size = 2) +
        stat_density_2d(aes(fill = ..level..), geom = "polygon", colour = "white", adjust = adjust) +
        scale_fill_distiller(palette = "Reds", direction = 1) +
        guides(fill = guide_none()) +
        facet_wrap(~split, ncol = ncol) +
        scale_x_continuous(expand = c(0.05, 0.05)) +
        scale_y_continuous(expand = c(0.05, 0.05)) +
        umap_aes() +
        facet_aes() +
        theme(
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            strip.text = element_text(face = "bold", size = 14)
            )

    return(plot)
}
