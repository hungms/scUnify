process_seurat <- function(x, assay = "RNA", nfeatures = 2000, dims = 1:10, res = 0.4, samples = NULL, vars.to.regress, ...){
    DefaultAssay(x) <- assay
    if(length(samples) == 0){
        x <- SCTransform(x, vst.flavor = "v2", 
            variable.features.n = nfeatures, 
            assay = assay, 
            vars.to.regress = vars.to.regress)}
    else{
        list <- SplitObject(x, split.by = samples)
        for(i in seq_along(list)){
            list[[i]] <- SCTransform(
                list[[i]], 
                vst.flavor = "v2", 
                variable.features.n = nfeatures, 
                assay = assay, 
                vars.to.regress = vars.to.regress, 
                return.only.var.genes = FALSE)}
        VariableFeatures(x) <- SelectIntegrationFeatures(list)}
    x <- RunPCA(x, npcs = 30)
    print(ElbowPlot(x, reduction = "pca"))
    x <- RunUMAP(x, dims = dims, reduction = "pca", verbose = F)
    x <- FindNeighbors(x, dims = dims, reduction = "pca", verbose = F)
    x <- FindClusters(x, resolution = res, verbose = F)
    print(scUMAP(x, reduction = "umap", group.by = "seurat_clusters"))
    if(!length(samples) == 0){
        print(scUMAP(x, reduction = "umap", group.by = "samples"))}
    return(x)}



plot_umap <- function(x, reduction = "umap", group.by, split.by = NULL, cols = NULL, count = T, pt.size = 0.5, legend.ncol = 1, shuffle = F){
    umap <- as.data.frame(Embeddings(x@reductions[[reduction]]))
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap$group <- x@meta.data[[group.by]]
    umap$group <- factor(umap$group, levels(x@meta.data[[group.by]]))

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