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