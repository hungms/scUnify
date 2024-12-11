#' calculate_fractions
#'
#' Calculate fractions of mitochondrial, ribosomal, haemoglobin, TCR, BCR and MHC genes
#' @param x Seurat object
#' @export
calculate_fractions <- function(x){
    x <- PercentageFeatureSet(x, pattern = "^[Mm][Tt]-", col.name = "pct.mt")
    x <- PercentageFeatureSet(x, pattern = "^R[Pp][SsLl]", col.name = "pct.rb")
    x <- PercentageFeatureSet(x, pattern = "^H[B][ABDEGMPQZ]?\\d*$|^H[b][abdegmpqz]?\\d*", col.name = "pct.hb")
    x <- PercentageFeatureSet(x, pattern = "^T[Rr][ABCDGabcdg][VDJCvdjc]", col.name = "pct.tcr")
    x <- PercentageFeatureSet(x, pattern = "^I[Gg][HKLhkl][VDJCAEMGvdjcaemg]", col.name = "pct.bcr")
    x <- PercentageFeatureSet(x, pattern =  "^HLA-|^H2-", col.name = "pct.mhc")
    return(x)}

#' run_htodemux
#'
#' HTODemux algorithm to demultiplex cell hashtags
#' @param x Seurat object
#' @param process if TRUE, run cluster algorithms on HTO counts
#' @param add.one +1 pseudocount for each hashtag oligo in each cell
#' @export
run_htodemux <- function(x, process = F, add.one = F){
    if("HTO" %in% names(x@assays)){
        DefaultAssay(x) <- "HTO"
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) + 1}
        x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR", margin = 2)
        x <- HTODemux(x, assay = "HTO", positive.quantile = 0.99)
        if(process){
            x <- process_hto(x)}
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) - 1}
        DefaultAssay(x) <- "RNA"}
    return(x)}

#' run_multiseqdemux
#'
#' MULTIseqDemux algorithm to demultiplex cell hashtags
#' @param x Seurat object
#' @param process if TRUE, run cluster algorithms on HTO counts
#' @param add.one +1 pseudocount for each hashtag oligo in each cell
#' @export
run_multiseqdemux <- function(x, process = F, add.one = F){
    if("HTO" %in% names(x@assays)){
        DefaultAssay(x) <- "HTO"
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) + 1}
        x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR", margin = 2)
        x <- MULTIseqDemux(x, assay = "HTO")
        x@meta.data$MULTI.global <- ifelse(x@meta.data$MULTI_ID == "Doublet", "Doublet", ifelse(x@meta.data$MULTI_ID == "Negative", "Negative", "Singlet"))
        if(process){
            x <- process_hto(x)}
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) - 1}
        DefaultAssay(x) <- "RNA"}
    return(x)}

#' process_hto
#'
#' Run clustering algorithms on HTO counts
#' @param x Seurat object
#' @param assay name of HTO assay, defaults to "HTO"
#' @export
process_hto <- function(x, assay = "HTO"){
    assayname <- tolower(assay)
    if(!is.list(x)){
        x <- list(x)}

    for(i in seq_along(x)){
        nfeatures <- nrow(x[[i]][[assay]])
        if(nfeatures < 5){
            next}
        features <- rownames(x[[i]][[assay]])
        DefaultAssay(x[[i]]) <- assay
        x[[i]] <- ScaleData(x[[i]], features = features, verbose = F)
        x[[i]] <- RunPCA(x[[i]], features = features, approx = FALSE, reduction.name = paste0("pca_", assayname), verbose = F)
        x[[i]] <- RunUMAP(x[[i]], reduction = paste0("pca_", assayname), dims = 1:nfeatures, reduction.name = paste0("umap_", assayname), verbose = F)}

    if(length(x) == 1){
        x <- x[[1]]}

    return(x)
}

#' run_dsb
#' 
#' Normalize CITEseq expressions by DSB normalization
#' 
#' @param x a seurat object
#' @param dir original directory
#' @param denoise.counts denoise counts
#' @param use.isotype.control use isotype controls
#' @param isotype.control.name.vec isotype control name vec
#' @param min.count minimum count required for each protein stained
#' @return a list of seurat objects with DSB normalised "data" in ADT assay
#' @export
run_dsb <- function(x, dir, denoise.counts = T, use.isotype.control = F, isotype.control.name.vec = NULL, min.count = NULL){
    stopifnot(length(x) == length(dir))
    raw.dir <- gsub("/outs/.*", "/outs/multi/count/raw_feature_bc_matrix", dir)

    adt.cell <- Read10X(dir)[["Antibody Capture"]][rownames(x[["ADT"]]),]
    adt.raw <- Read10X(raw.dir)[["Antibody Capture"]][rownames(x[["ADT"]]),]
    
    if(length(min.count) > 0){
        if(is.numeric(min.count)){
            keep1 <- rownames(adt.cell)[which(apply(adt.cell, 1, max) > min.count)]
            keep2 <- rownames(adt.raw)[which(apply(adt.raw, 1, max) > min.count)]
            keep <- intersect(keep1, keep2)
            adt.cell <- adt.cell[keep,]
            adt.raw <- adt.raw[keep, -c(which(colnames(adt.raw) %in% colnames(adt.cell)))]}}

    adt.raw <- adt.raw[, -c(which(colnames(adt.raw) %in% colnames(adt.cell)))]
    out = dsb::DSBNormalizeProtein(
            cell_protein_matrix = adt.cell, 
            empty_drop_matrix = adt.raw, 
            denoise.counts = denoise.counts, 
            use.isotype.control = use.isotype.control,
            isotype.control.name.vec = isotype.control.name.vec)
    x[["ADT"]]$data <- as(out, "sparseMatrix")
    
    return(x)
    }

#' calculate_mad
#'
#' Calculate median absolute deviation (MAD) and identify poor quality cells
#' @param x Seurat object
#' @param columns a vector of metadata columns to calculate MAD
#' @param samples a metadata column of sample names
#' @param stdev label poor quality cells as stdev*MAD away from median for each column. Defaults to 5
#' @export
calculate_mad <- function(x, columns, samples, stdev = 5){
    mad_output <- list()
    stopifnot(is.vector(columns))

    for(i in unique(x@meta.data[[samples]])){
        data <- x@meta.data[which(x@meta.data[[samples]] == i),]

        qclist <- list()
        for(j in seq_along(columns)){
            selected <- data[[columns[j]]]
            median <- median(selected)
            mad <- median(abs(selected - median))
            qclist[[j]] <- ifelse(selected >= (median - stdev*mad) & selected <= (median + stdev*mad), "Pass", "Fail")
            names(qclist)[j] <- columns[j]}

        qclist <- as.data.frame(qclist)
        qclist$softqc <- apply(qclist, 1, function(row) ifelse(all(row == "Pass"), "Pass", "Fail"))
        rownames(qclist) <- rownames(data)
        mad_output[[i]] <- qclist}
    
    mad_output <- bind_rows(mad_output)[colnames(x),]
    x@misc$softqc <- mad_output
    x@meta.data$softqc <- mad_output$softqc
	
    qc_report(x, column = "softqc", samples = samples)
    return(x)
    }

#' qc_report
#'
#' Report cells remaining from each sample after quality control
#' @param x Seurat object
#' @param column a metadata column with values of "Pass" or "Fail"
#' @param samples a metadata column of sample names
#' @export
qc_report <- function(x, column, samples){

    for(i in unique(x@meta.data[[samples]])){
        data <- x@meta.data[which(x@meta.data[[samples]] == i),]
        n = length(which(data[[column]] == "Pass"))
        pct = n*100/nrow(data)
        pct = round(pct, 1)
        message(paste0(pct, "% (", n, ") of cells remains - ", i))}
        }

#' run_doubletfinder
#'
#' predict doublets with DoubletFinder
#' @param x Seurat object
#' @param assay assay name, defaults to "RNA"
#' @param dims no. of PCs, defaults to 1:10
#' @param truth whether cell demultiplexing (ground-truth) result is available. Defaults to NULL, else set to metadata column with values of "Singlet" and "Doublet".
#' @param clusters metadata column containing cluster labels
#' @param dbr doublet rate
#' @param ncores cores used for parallel processing, defaults to 1
#' @export
run_doubletfinder <- function(x, assay = "RNA", dims = 1:10, truth = NULL, clusters, dbr, ncores = 1){

    x@meta.data <- x@meta.data %>%
        mutate(cell_barcode = rownames(.))
    obj <- x

    if(!"data" %in% Layers(obj, assay = assay)){
        obj <- NormalizeData(obj, verbose = F)}

    if(str_detect(assay, "SCT")){
        sct <- TRUE}
    else{
        sct <- FALSE}

    if(length(truth) == 0){
        message("Step 1 : No ground truth")
        suppressMessages({
            sweep_res <- paramSweep(obj, PCs = dims, sct = sct, num.cores = ncores)
            sweep_stat <- summarizeSweep(sweep_res, GT = FALSE)})}

    else if (truth %in% colnames(obj@meta.data)){
        message("Step 1 : Filter singlet and doublet from ground truth")
        cells <- colnames(obj)[which(obj@meta.data[[truth]] %in% c("Singlet", "Doublet"))]
        obj <- subset(obj, cells = cells)

        suppressMessages({
            sweep_res <- paramSweep(obj, PCs = dims, sct = sct, num.cores = ncores)
            gt.calls <- obj@meta.data[rownames(sweep_res[[1]]), truth]
            sweep_stat <- summarizeSweep(sweep_res, GT = TRUE, GT.calls = gt.calls)})}

    bcmvn <- find.pK(sweep_stat)
    top_pK <- bcmvn %>%
        mutate(pK = as.numeric(as.character(pK))) %>%
        filter(pK <= 0.1) %>%
        top_n(1, BCmetric) %>%
        .$pK
    message(paste0("Value of top pK is ", top_pK))

    message("Step 2 : Estimate Doublets")
    prop <- modelHomotypic(x@meta.data[[clusters]])
    nExp_poi <- round(dbr*nrow(x@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-prop))
    obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = top_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
    obj@meta.data %>% dplyr::select(starts_with("pANN_0.25")) %>% colnames() -> pann
    obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = top_pK, nExp = nExp_poi.adj, reuse.pANN = pann, sct = sct)
    class <- obj@meta.data %>%
        dplyr::select(c("cell_barcode", tail(colnames(obj@meta.data),1)))
    colnames(class)[2] <- "DoubletFinder"

    x@meta.data <- x@meta.data %>%
	as.data.frame(.) %>%
        merge(., class, by = "cell_barcode", all.x = T) %>%
        column_to_rownames("cell_barcode") %>%
        as.data.frame(.)
    return(x)
}

#' run_scdblfinder
#'
#' predict doublets with scDblFinder
#' @param x Seurat object
#' @param assay assay name, defaults to "RNA"
#' @param truth whether cell demultiplexing (ground-truth) result is available. Defaults to NULL, else set to metadata column with values of "Singlet" and "Doublet".
#' @param samples metadata column containing sample labels
#' @param clusters metadata column containing cluster labels
#' @param dbr doublet rate
#' @param ncores cores used for parallel processing, defaults to 1
#' @export
run_scdblfinder <- function(x, assay = "RNA", samples, truth = NULL, clusters, dbr, ncores = 1){
    DefaultAssay(x) <- assay
    sce <- convert_seurat_to_sce(x)

    if(length(truth) > 0){
        sce <- sce[,colData(sce)[[truth]] %in% c("Singlet", "Doublet")]}

    sce <- scDblFinder(sce, clusters = clusters, samples = samples, dbr=dbr, BPPARAM=BiocParallel::MulticoreParam(ncores))
    sce$scDblFinder <- ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet")
    sce$cell_barcode <- colnames(sce)
    sce.meta <- colData(sce)[,c("cell_barcode", "scDblFinder")]

    x@meta.data <- x@meta.data %>%
        as.data.frame(.) %>%
        rownames_to_column("cell_barcode") %>%
        merge(., as.data.frame(sce.meta), by = "cell_barcode", all.x = T) %>%
        column_to_rownames("cell_barcode") %>%
        as.data.frame(.)
    return(x)
}

#' calculate_cellcycle
#'
#' Calculate cell cycle score and predict cell cycle phase by Seurat
#' @param x Seurat object
#' @param org organism, either "human" or "mouse"
#' @param remove_genes If TRUE, remove genes from current assay after calculating cell cycle phase. Defaults to FALSE 
#' @param assay assay name, defaults to "RNA"
#' @export
calculate_cellcycle <- function(x, org = "human", remove_genes = F, assay = "RNA"){

    if(!org %in% c("human", "mouse")){
        stop('organism must be either "human" or "mouse"')}
    sgenes <- cc.genes$s.genes
    g2mgenes <- cc.genes$g2m.genes
    if(org == "mouse"){
        sgenes <- convert_human_to_mouse(cc.genes$s.genes, unique = T)
        g2mgenes <- convert_human_to_mouse(cc.genes$g2m.genes, unique = T)}
    
    x <- CellCycleScoring(x, s.features = sgenes, g2m.features = g2mgenes, set.ident = F, assay = assay)
    
    if(remove_genes){
        x <- remove_genes(x, features = c(sgenes, g2mgenes), from.assay = assay, to.assay = "CC")}
    return(x)    
    }


#' process_seurat
#'
#' Wrapper function to cluster cells by gene expression
#' @param x Seurat object
#' @param assay assay name of gene expression. Defaults to "RNA"
#' @param nfeatures number of variable features to use. Defaults to 2000
#' @param dims no. of PCs to use. Defaults to 1:10
#' @param res clustering resolution. Defaults to 0.4
#' @param samples sample column
#' @param vars.to.regress variables to regress
#' @export
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

#' integrate_v4
#'
#' Wrapper function to integrate cells by samples
#' @param x Seurat object
#' @param split.by variable to split object by
#' @param assay assay name of gene expression. Defaults to "RNA"
#' @param nfeatures no. of features to select for integration
#' @param npcs no. of PCs to computec
#' @param vars.to.regress variables to regress
#' @param method integration method, see Seurat::indIntegrationAnchors()
#' @param k.weight k.weight, see Seurat::IntegrateData(), reduce when no. of cells in samples are low.
#' @export
integrate_v4 <- function(x, split.by, assay = "RNA", nfeatures = 3000, npcs = 50, vars.to.regress = NULL, method = "rpca", k.weight = 100){
    options(Seurat.object.assay.version = "v4")
    list <- SplitObject(x, split.by = split.by)
    for(i in seq_along(list)){
        list[[i]] <- SCTransform(list[[i]], vst.flavor = "v2", variable.features.n = nfeatures, assay = assay, vars.to.regress = vars.to.regress, return.only.var.genes = FALSE)
        list[[i]] <- RunPCA(list[[i]], verbose = FALSE, assay = "SCT", npcs = npcs)}
    features <- SelectIntegrationFeatures(list, nfeatures = nfeatures)
    list <- PrepSCTIntegration(object.list = list, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features, reduction = method)
    integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = rownames(x[[assay]]), k.weight = k.weight)
    VariableFeatures(integrated[["SCT"]]) <- features
    VariableFeatures(integrated[["integrated"]]) <- features
    options(Seurat.object.assay.version = "v5")
    return(integrated)}


#' qc_before_recluster
#'
#' Subset cells under < 200 features after QC
#' @param x Seurat object
#' @param assay assay name of gene expression. Defaults to "RNA"
#' @param min.features minimum genes expressed by each cell. Defaults to 200.
#' @param min.cells minimum cells expressing each gene. Defaults to 3.
#' @export
qc_before_recluster <- function(x, assay = "RNA", min.features = 200, min.cells = 3){
    x[[assay]] <- CreateAssay5Object(counts = x[[assay]]$counts, data = x[[assay]]$data, min.features = min.features, min.cells = min.cells)
    x <- subset(x, cells = colnames(x[[assay]]$counts))
    return(x)}

#' calculate_correlation
#'
#' calculate pearson correlation of gene expression between columns of gene expression matrix
#' @param matrix gene expression matrix, clusters (columns) by features (rows)
#' @export
calculate_correlation <- function(matrix){
    id <- colnames(matrix)
    ncol <- ncol(matrix)
    correlation.matrix <- matrix(NA, nrow = ncol, ncol = ncol)
    rownames(correlation.matrix) <- id
    colnames(correlation.matrix) <- id
    for (i in 1:ncol) {
    for (j in 1:ncol) {
        if (i < j) {
        pearson_cor <- cor(x = matrix[,i], y = matrix[,j])
        correlation.matrix[i, j] <- pearson_cor
        correlation.matrix[j, i] <- pearson_cor}}}
    correlation.matrix <- replace(correlation.matrix, is.na(correlation.matrix), 1)
    return(correlation.matrix)}

#' calculate_cluster_similarity
#'
#' calculate pearson correlation similarity between seurat clusters
#' @param x Seurat object
#' @param cluster cluster metadata column containing cluster labels
#' @param assay assay name for gene expression. Defaults to "RNA"
#' @param slot slot name for gene expression. Defaults to "data"
#' @param variable.features if TRUE, calculate correlation of variable features only
#' @export
calculate_cluster_similarity <- function(x, cluster, assay = "RNA", slot = "data", variable.features = F){
    if(variable.features){
        selected.features <- VariableFeatures(x)}
    else{
	selected.features <- rownames(x[[assay]])}

    pseudobulk <- AverageExpression(x, assays = assay, features = selected.features, group.by = cluster, slot = slot)
    pseudobulk <- pseudobulk[[assay]]
    correlation.matrix <- calculate_correlation(pseudobulk)

    return(correlation.matrix)}