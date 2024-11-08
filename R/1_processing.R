#' create_seurat_object
#'
#' Create Seurat object from cellranger-multi outputs
#' @param dir a vector of directories
#' @param samples a vector of sample names for each directory
#' @param hto_str prefix to identify HTO tag names if HTO library is present. Defaults to NULL
#' @export
create_seurat_object <- function(dir, samples, hto_str = NULL, adt_normalize = T){
    stopifnot(length(samples) == length(dir))

    x <- list()
    for(i in seq_along(samples)){

       	# read 10x output
        message(paste0(samples[i], " --- Loading Sample ", i))
        message("Step 1 : Adding RNA counts")
        mtx <- Read10X(dir[i])

        # add rna counts
        any.adt <- "Antibody Capture" %in% names(mtx)
        if(any.adt){rna <- mtx$`Gene Expression`}
        else{rna <- mtx}
        seu_obj <- suppressWarnings(CreateSeuratObject(as.matrix(rna), assay="RNA", min.cells = 3, min.features = 200)) # import gene expression/adt to seurat, remove cells with less than 200 features and features expressed in less than 3 cells

        # check for modalities
        message("Step 2 : Check for Modalities")
        if(any.adt){
            any.hto <- any(str_detect(rownames(mtx$`Antibody Capture`), hto_str))
            any.adt <- any(!(str_detect(rownames(mtx$`Antibody Capture`), hto_str)))}
        else{
            any.hto <- F}
        adt.skip <- ifelse(any.adt, "", "--- SKIPPED")
        hto.skip <- ifelse(any.hto, "", " --- SKIPPED")

       	# add hto counts
        message(paste0("Step 3 : Adding HTO counts", hto.skip))
        if(any.hto){
            hash <- mtx$`Antibody Capture`[c(which(str_detect(rownames(mtx$`Antibody Capture`), hto_str))),]
            seu_obj[["HTO"]] <- CreateAssay5Object(as.matrix(hash[,colnames(seu_obj[["RNA"]]$counts)]), min.cells = 0, min.features = 0)
            }

       	# add adt counts
        message(paste0("Step 4 : Adding ADT counts", adt.skip))
        if(any.adt){
            prot <- mtx$`Antibody Capture`
            if(any.hto){
                prot <- prot[-c(which(str_detect(rownames(prot), hto_str))),]}
            seu_obj[["ADT"]] <- CreateAssay5Object(as.matrix(prot[,colnames(seu_obj[["RNA"]]$counts)]), min.cells = 0, min.features = 0)
            if(adt_normalize){
                seu_obj <- dsb_normalize(seu_obj, dir = dir[i], denoise.counts = T, use.isotype.control = F, isotype.control.name.vec = NULL)}

            }

       	# add sample id in metadata
        message("Step 5 : Creating Seurat Object")
        seu_obj <- set_assay_keys(seu_obj)
        seu_obj@meta.data$samples <- rep(samples[i], ncol(seu_obj))
        seu_obj <- RenameCells(seu_obj, new.names = paste0(samples[i], "_", colnames(seu_obj)))
        x[[i]] <- seu_obj}

    names(x) <- paste0(samples)
    return(x)}

#' calculate_fractions
#'
#' Calculate fractions of mitochondrial, ribosomal, haemoglobin, TCR, BCR and MHC genes
#' @param x Seurat object
#' @export
calculate_fractions <- function(x){
    x <- PercentageFeatureSet(x, pattern = "^[Mm][Tt]-", col.name = "pct.mt")
    x <- PercentageFeatureSet(x, pattern = "^R[Pp][SsLl]", col.name = "pct.rb")
    x <- PercentageFeatureSet(x, pattern = "^H[B][ABDEGPQ]?\\d*$|^H[b][abdegpq]?\\d*", col.name = "pct.hb")
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

#' dsb_normalize
#' 
#' Normalize CITEseq expressions by DSB normalization
#' 
#' @param x a seurat object
#' @param dir original directory
#' @param denoise.counts denoise counts
#' @param use.isotype.control use isotype controls
#' @param isotype.control.name.vec isotype control name vec
#' @return a list of seurat objects with DSB normalised "data" in ADT assay
#' @export
dsb_normalize <- function(x, dir, denoise.counts = T, use.isotype.control = F, isotype.control.name.vec = NULL){
    stopifnot(length(x) == length(dir))
    raw.dir <- gsub("/outs/.*", "/outs/multi/count/raw_feature_bc_matrix", dir)

    adt.cell <- Read10X(dir)[["Antibody Capture"]][rownames(x[["ADT"]]),]
    adt.raw <- Read10X(raw.dir)[["Antibody Capture"]][rownames(x[["ADT"]]),]
    
    keep1 <- rownames(adt.cell)[which(apply(adt.cell, 1, max) > 10)]
    keep2 <- rownames(adt.raw)[which(apply(adt.raw, 1, max) > 10)]
    keep <- intersect(keep1, keep2)
    
    adt.cell <- adt.cell[keep,]
    adt.raw <- adt.raw[keep, -c(which(colnames(adt.raw) %in% colnames(adt.cell)))]
    x[["ADT"]]$data = dsb::DSBNormalizeProtein(
            cell_protein_matrix = adt.cell, 
            empty_drop_matrix = adt.raw, 
            denoise.counts = denoise.counts, 
            use.isotype.control = use.isotype.control,
            isotype.control.name.vec = isotype.control.name.vec)
    
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
        x <- remove_genes(x, features = c(sgenes, g2mgenes), orig.assay = assay, new.assay = "CC")}
    return(x)    
    }

#' remove_genes
#'
#' Remove genes from current assay to new assay
#' @param x Seurat object
#' @param features a vector genes/features to remove
#' @param orig.assay current assay name, defaults to "RNA"
#' @param new.assay new assay name to store removed genes
#' @export
remove_genes <- function(x, features = NULL, orig.assay = "RNA", new.assay){

    if(orig.assay == "RNA"){
        stopifnot("data" %in% Layers(x, assay = "RNA"))}

    DefaultAssay(x) <- orig.assay
    keep <- which(rownames(x) %in% c(features))


    feature.counts <- x[[orig.assay]]$counts[keep, colnames(x)]
    new.counts <- x[[orig.assay]]$counts[-c(keep), colnames(x)]

    if("data" %in% Layers(x, assay = orig.assay)){
        feature.data <- x[[orig.assay]]$data[keep, colnames(x)]
        new.data <- x[[orig.assay]]$data[-c(keep), colnames(x)]
        x[[new.assay]] <- CreateAssay5Object(counts = as(feature.counts, "sparseMatrix"), data = as(feature.data, "sparseMatrix"), min.features = 0, min.cells = 0)
        x[[orig.assay]] <- CreateAssay5Object(counts = as(new.counts, "sparseMatrix"), data = as(new.data, "sparseMatrix"), min.features = 0, min.cells = 0)}
    else{
        x[[new.assay]] <- CreateAssay5Object(counts = as(feature.counts, "sparseMatrix"), min.features = 0, min.cells = 0)
        x[[orig.assay]] <- CreateAssay5Object(counts = as(new.counts, "sparseMatrix"), min.features = 0, min.cells = 0)}

    return(x)}

#' return_genes
#'
#' Return genes from an assay to a specific assay
#' @param x Seurat object
#' @param from.assay assay name from which genes/features are taken
#' @param to.assay assay name to which genes/features are stored, defaults to "RNA"
#' @export
return_genes <- function(x, from.assay, to.assay = "RNA"){
    counts <- rbind(x[[from.assay]]$counts, x[[to.assay]]$counts)
    data <- rbind(x[[from.assay]]$data, x[[to.assay]]$data)
    x[[to.assay]] <- CreateAssay5Object(counts = counts, data = data, min.cells = 0, min.features = 0)
    x[[from.assay]] <- NULL
    return(x)}

#' remove_vdj_genes
#'
#' Remove VDJ genes from an assay
#' @param x Seurat object
#' @param bcr if TRUE, remove BCR-VDJ genes. Defaults to FALSE
#' @param tcr if TRUE, remove TCR-VDJ genes. Defaults to FALSE
#' @param orig.assay current assay name, defaults to "RNA"
#' @export
remove_vdj_genes <- function(x, bcr = T, tcr = T, orig.assay = "RNA"){
    if(bcr){
        bcr.genes <- rownames(x[[orig.assay]])[which(str_detect(rownames(x[[orig.assay]]), bcr.string))]
        x <- remove_genes(x, features = bcr.genes, orig.assay = orig.assay, new.assay = "BCR")}
    if(tcr){
        tcr.genes <- rownames(x[[orig.assay]])[which(str_detect(rownames(x[[orig.assay]]), tcr.string))]
        x <- remove_genes(x, features = tcr.genes, orig.assay = orig.assay, new.assay = "TCR")}
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
#' @param vars.to.regress variables to regress
#' @param method integration method, see Seurat::indIntegrationAnchors()
#' @param k.weight k.weight, see Seurat::IntegrateData(), reduce when no. of cells in samples are low.
#' @export
integrate_v4 <- function(x, split.by, assay = "RNA", nfeatures = 3000, vars.to.regress = NULL, method = "rpca", k.weight = 100){
    list <- SplitObject(x, split.by = split.by)
    for(i in seq_along(list)){
        list[[i]] <- SCTransform(list[[i]], vst.flavor = "v2", variable.features.n = nfeatures, assay = assay, vars.to.regress = vars.to.regress, return.only.var.genes = FALSE)
        list[[i]] <- RunPCA(list[[i]], npcs = 50, verbose = FALSE, assay = "SCT")}
    features <- SelectIntegrationFeatures(list, nfeatures = nfeatures)
    list <- PrepSCTIntegration(object.list = list, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features, reduction = method)
    integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = rownames(x[[assay]]), k.weight = k.weight)
    VariableFeatures(integrated[["SCT"]]) <- features
    VariableFeatures(integrated[["integrated"]]) <- features
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