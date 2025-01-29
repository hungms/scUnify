#' theme_text
#' 
#' ggplot2 aesthetic option with texts
#' @export 
theme_text <- function(){
    theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8, size = 10),
        axis.text.y = element_text(size = 8))}

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


#' plot_hvf
#'
#' Plot highly variable features from Seurat object
#' @param x Seurat object
#' @param assay assay name
#' @param nfeatures no. of features to select
#' @param n number of features to label
#' @export
plot_hvf <- function(x, assay = "RNA", nfeatures = 3000, n = 10){
    DefaultAssay(x) <- assay
    x <- FindVariableFeatures(x, nfeatures = nfeatures)
    plot <- VariableFeaturePlot(x) +
        geom_text_repel(data = VariableFeaturePlot(x)$data %>% 
            rownames_to_column("gene") %>%
            slice_max(n = n, order_by = variance.standardized) %>%
            as.data.frame(), aes(label = gene))
    return(plot)}

#' plot_umap
#'
#' Plot UMAP from Seurat object
#' @param x Seurat object
#' @param reduction reduction name, defaults to "umap"
#' @param group.by metadata column to group by 
#' @param split.by metadata column to split by
#' @param cols a vector of colors, defaults to NULL
#' @param count display no. of cells if TRUE
#' @param pt.size geom_point size
#' @param legend.ncol no. of columns for group.by keys in legend
#' @param shuffle shuffle points
#' @export
plot_umap <- function(x, reduction = "umap", group.by, split.by = NULL, cols = NULL, count = T, pt.size = 0.5, legend.ncol = 1, shuffle = F){
    umap <- as.data.frame(Embeddings(x@reductions[[reduction]]))
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap$group <- x@meta.data[[group.by]]
    if(is.factor(x@meta.data[[group.by]])){
        umap$group <- factor(umap$group, levels(x@meta.data[[group.by]]))}
    else{
        umap$group <- factor(umap$group, sort(unique(x@meta.data[[group.by]])))}

    if(count){
        umap <- umap %>%
            group_by(group) %>%
            mutate(groupc = factor(paste0(group, " (", n(), ") "))) %>%
            ungroup() %>%
            arrange(group) %>%
            mutate(group = factor(groupc, unique(groupc)))}


    if(length(split.by) > 0){
        umap$split <- x@meta.data[[split.by]]
        umap$split <- factor(umap$split, levels(x@meta.data[[split.by]]))}

    if(shuffle){
        umap <- umap %>% sample_frac(1)}

    plot <- ggplot(umap, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(color = group), size = pt.size) +
        guides(color = guide_legend(title = "", override.aes = list(size = 3), ncol = legend.ncol)) +
        scale_x_continuous(expand = c(0.05, 0.05)) +
        scale_y_continuous(expand = c(0.05, 0.05)) +
        umap_aes() +
        facet_aes() +
        theme(
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            strip.text = element_text(face = "bold", size = 14))

    # count
    if(count){
        dflabel <- umap %>%
            summarize(count = paste0("n = ", n()))
        dflabel[["UMAP1"]] <- max(umap[["UMAP1"]])
        dflabel[["UMAP2"]] <- min(umap[["UMAP2"]])
        plot <- plot + 
            geom_text(data = dflabel, aes(label = count), color = "black", vjust="inward", hjust="inward")}
    if(length(split.by) > 0){
        plot <- plot + 
            facet_wrap(~split, nrow = 1)}
    if(length(cols) > 0){ # & all(unique(umap$group) %in% names(cols))
        names(cols) <- levels(umap$group)
        plot <- plot + scale_color_manual(values = cols)}
    return(plot)
}


#' plot_dotplot
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
plot_dotplot <- function(
    x, group.by, split.by = NULL, features,
    assay = "RNA", slot = "data", scale = T,
    palette = "RdBu", direction = -1, diffexp = NULL,
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
    palette <- get_palette(palette, n = 9, direction = direction)

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
                strip.text.y = element_text(face = "bold", size = 10),
                strip.background.x = element_blank(),
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
            theme(
                strip.text.x = element_text(face = "bold", size = 10),
                strip.background.y = element_blank(),
                strip.placement = "outside")} 

    plot <- plot  +
        theme(
            legend.text = element_text(size = 8),         # Reduce text size
            legend.title = element_text(size = 10, face = "bold"),      # Reduce title size
            legend.spacing.x =unit(0, "cm"),
            legend.spacing.y =unit(0, "cm")
            #legend.key.size = unit(0.3, "cm")            # Reduce key size
            )

    return(plot)
}

basicDotPlot <- function(df, x, y, fill, size, stroke_legend = F, palette, scale = T){
    if("signif" %in% colnames(df)){
        stroke <- "signif"
        stroke_legend <- guide_legend(
            direction = "vertical",
            order = 3,
            keyheight = unit(0.3, "cm"),  # Scale down key height
            keywidth = unit(0.3, "cm"),   # Scale down key width
            override.aes = list(fill = palette[6], size = 6),
            theme = theme(legend.text=element_text(size=10)),
            title = "FDR of\nupregulated genes",
            title.position = "top")}
    else{
	stroke <- NULL
        stroke_legend <- guide_none()}

    if(scale){
        colorbar_title <- "Scaled Expression"}
    else{
        colorbar_title <- "Average Expression"}

    # basic dotplot
    plot <- ggplot(df, aes_string(x = x, y = y)) +
        geom_point(aes_string(size = size, fill = fill, stroke = stroke), shape=21) +
        scale_discrete_manual(aesthetics = "stroke", values = c(1.3, 0)) +
        scale_color_manual(values = c("black", "white")) +
        scale_fill_gradientn(colors = palette) +
        #scale_size_continuous("% detected", range = c(0,8), breaks = c(0, 33, 66), labels = c("0%", "33%", "66%"))  +
        scale_size(range = c(0, 10), breaks = c(25, 50, 75), labels = c("25%", "50%", "75%")) +
        ylab("") +
        xlab("") +
        guides(
            fill = guide_colorbar(
                title = colorbar_title,
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
        #facet_aes() +
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
                avg_log2FC > 0 & p_val_adj < 0.05 & pct.1 > 0.1 ~ "< 0.05",
                .default = "ns"),
            signif = factor(signif, c("< 0.05", "ns")),
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
    
#' plot_feature
#'
#' Plot FeaturePlot from Seurat object
#' @param x Seurat object
#' @param features feature names stored in a vector
#' @param palette a vector of hexcode colors, or palettes from "RColorBrewer" or "viridis"
#' @param direction direction of the palette vector
#' @param ncol no. of columns for the plots
#' @param plot If TRUE, plot with cowplot::plot_grid(). Defaults to FALSE to return a list of ggplot objects
#' @param ... arguments for Seurat::FeaturePlot()
#' @export
plot_feature <- function(x, assay = "RNA", features, palette = "viridis", direction = 1, ncol, plot = F, ...){
    suppressMessages({
        DefaultAssay(x) <- assay
        plotlist <- FeaturePlot(x, features = features, ncol = ncol, combine = F, ...)
        palette <- get_palette(palette, n = 9, direction = direction)
        for(x in seq_along(plotlist)){
            plotlist[[x]] <- plotlist[[x]] + 
                scale_color_gradientn(colors = palette) +
                guides(color = guide_colorbar(title = "")) +
		umap_aes()}})
    p <- plotlist
    if(plot){
        p <- plot_grid(plotlist = plotlist, ncol = ncol, align = "hv")}
    return(p)}


#' plot_density_umap
#'
#' Plot DensityUMAP from Seurat object
#' @param x Seurat object
#' @param reduction reduction name, defaults to "umap"
#' @param split.by metadata column to group by
#' @param adjust adjust for stat_density_2d() or geom_density_2d
#' @param contour plot bins for geom_density_2d if TRUE
#' @param together plot cells from all facet in background if TRUE
#' @param count display the number of cells if TRUE
#' @param pt.size geom_point size
#' @export
plot_density_umap <- function(x, reduction = "umap", split.by, adjust = 1, contour = T, together = F, count = T, pt.size = 0.5){
    umap <- as.data.frame(Embeddings(x@reductions[[reduction]]))
    umap2 <- umap
    umap$split <- x@meta.data[[split.by]]
    ncol <- length(unique(umap$split))
    colnames <- colnames(umap)

    plot <- ggplot(umap, aes_string(x = colnames[1], y = colnames[2])) +
        geom_point(color = "grey80", size = pt.size) +
        guides(color = guide_colorbar(title = "Density"), fill = guide_none()) +
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

    # count
    if(count){
        dflabel <- umap %>%
            group_by(split) %>%
            summarize(count = paste0("n = ", n())) %>%
            ungroup()
        dflabel[[colnames[1]]] <- max(umap[[colnames[1]]])
        dflabel[[colnames[2]]] <- min(umap[[colnames[2]]])
        plot <- plot + 
            geom_text(data = dflabel, aes_string(label = "count"), vjust="inward",hjust="inward")}

    # together
    if(together){
        plot <- plot + 
            geom_point(data = umap2, color = "grey80", size = pt.size)}

    # contour 
    if(contour){
        plot <- plot + 
            geom_density_2d(aes(color = ..level..), bins = adjust)}
    else{
        plot <- plot + 
            stat_density_2d(aes(fill = ..level..), geom = "polygon", colour = "white", adjust = adjust) +
            scale_fill_distiller(palette = "Reds", direction = 1)}
    return(plot)
}

#' plot_sc_quality
#'
#' quality control plot for single cell
#' @param x Seurat object
#' @param group.by column to group by
#' @param assay assay to use
#' @param cols colors
#' @export
plot_sc_quality <- function(x, group.by, assay = "RNA", cols = NULL){
    stopifnot(all(c(paste0("nFeature_", assay), paste0("nCount_", assay), "pct.mt", "pct.rb") %in% colnames(x@meta.data)))
    p1 <- x@meta.data %>%
        ggplot(aes_string(x = group.by, y = paste0("nFeature_", assay), fill = group.by)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = cols) +
        guides(fill = guide_none()) +
        xlab("") +
        ylab(paste0("nFeature_", assay)) +
        theme_line() +
        theme_text()

    p2 <- x@meta.data %>%
        ggplot(aes_string(x = group.by, y = paste0("nCount_", assay), fill = group.by)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = cols) +
        guides(fill = guide_none()) +
        xlab("") +
        ylab(paste0("nCount_", assay)) +
        theme_line() +
        theme_text()

    p3 <- x@meta.data %>%
        ggplot(aes_string(x = group.by, y = paste0("pct.mt"), fill = group.by)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = cols) +
        guides(fill = guide_none()) +
        xlab("") +
        ylab("Mitochondrial\nFraction (%)") +
        theme_line() +
        theme_text()

    p4 <- x@meta.data %>%
        ggplot(aes_string(x = group.by, y = paste0("pct.rb"), fill = group.by)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = cols) +
        guides(fill = guide_none()) +
        xlab("") +
        ylab("Ribosomal\nFraction (%)") +
        theme_line() +
        theme_text()
    
    plot <- plot_grid(p1, p2, p3, p4, ncol = 2)
    return(plot)
}


#' plot_similarity_heatmap
#'
#' plot similarity heatmap between seurat clusters
#' @param correlation.matrix a matrix of pearson correlation score
#' @param annotations a vector of cluster annotations for ComplexHeatmap::HeatmapAnnotation(). Defaults to colnames(correlation.matrix)
#' @param color colors for cluster annotations
#' @export
plot_similarity_heatmap <- function(correlation.matrix, annotations = colnames(correlation.matrix), color = NULL, cluster = T){

    col_fun = rev(brewer.pal(12,"RdBu"))
    if(all(annotations == colnames(correlation.matrix))){
        ha <- NULL}
    else{
        ha <- HeatmapAnnotation(
            Cluster = annotations,
            col = list(Cluster = color),
            show_annotation_name = F)}
    ht <- Heatmap(correlation.matrix,
        name = "Pearson\nCorrelation",
        top_annotation = ha,
        border = T,
        col = col_fun,
        na_col = "black",
        column_title_gp = gpar(fontsize = 12),
        border_gp = gpar(col = "black", lwd = 3),
        rect_gp = gpar(col = "white", lwd = 1),
        row_title_gp = gpar(fontsize = 12),
        cluster_rows = cluster,
        cluster_columns = cluster,
        heatmap_legend_param = list(
            title = "Pearson\nCorrelation",
            legend_direction = "vertical")
        )
    return(ht)}

#' plot_percent
#'
#' plot percentage for cells
#' @param x Seurat object
#' @param group.by column to group by
#' @param variable variable to color by
#' @param facet.by column to facet by
#' @param cols vector of colors
#' @param legend.ncol no. of legend columns
#' @export
plot_percent <- function(x, group.by, variable, facet.by = NULL, cols = NULL, legend.ncol = 1){
    
    stopifnot(c(group.by, variable, facet.by) %in% colnames(x@meta.data))

    group <- c(group.by, facet.by)
    metadata <- x@meta.data

    plot <- metadata %>%
        filter(!is.na(!!sym(variable))) %>%
        group_by_at(c(group, variable)) %>%
        summarize(count = n()) %>%
        group_by_at(group) %>%
        mutate(pct = count*100/sum(count)) %>%
        ggplot(aes_string(x = group.by, y = "pct", fill = variable)) +
        geom_col(width = 0.85, position = "stack", col = "white") +
        guides(fill = guide_legend(title = "", ncol = legend.ncol)) +
        theme_line() +
        theme_text() +
        xlab("") +
        ylab("Frequency (%)")
    
    if(length(cols) > 0){
        plot <- plot +
            scale_fill_manual(values = cols)}

    if(length(facet.by) == 1){
        plot <- plot +
            facet_wrap(as.formula(paste0("~ ", facet.by)), nrow = 1) +
            facet_aes()}

    return(plot)
}