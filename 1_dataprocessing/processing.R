########## SETUP ###########
# load R libraries
.libPaths("/nemo/lab/caladod/working/Matthew/.conda/envs/seurat5/lib/R/library")
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e12)
options(dplyr.summarise.inform = FALSE)
set.seed(123)

bcr.string <- "^I[Gg][HKLhkl][VDJCAEMGvdjcaemg]|^AC233755"
tcr.string <- "^T[Rr][ABCDGabcdg][VDJCvdjc]"

########### LOAD FUNCTIONS #############
create_seurat_object <- function(dir, samples, adt = F, hto = F, hto_str = NULL){
    stopifnot(length(samples) == length(dir))

    x <- list()
    for(i in seq_along(samples)){

       	# read 10x output
        message(paste0(samples[i], " --- Loading Sample ", i))
        message("Step 1 : Adding RNA counts")
        mtx <- Read10X(dir[i])
        if(adt == T | hto == T){
             rna <- mtx$`Gene Expression`}
        else{
             rna <- mtx}

       	# import gene expression/adt to seurat, remove cells with less than 200 features and features expressed in less than 3 cells
        seu_obj <- suppressWarnings(CreateSeuratObject(as.matrix(rna), assay="RNA", min.cells = 3, min.features = 200))

       	# add hto counts
        if(hto){
            message("Step 2 : Adding HTO counts")
            if(length(hto_str) == 0){stop("please provide HTO string pattern")}
            if(any(str_detect(rownames(mtx$`Antibody Capture`), hto_str))){
                hash <- mtx$`Antibody Capture`[c(which(str_detect(rownames(mtx$`Antibody Capture`), hto_str))),]
                seu_obj[["HTO"]] <- CreateAssay5Object(as.matrix(hash[,colnames(seu_obj[["RNA"]]$counts)]), min.cells = 0, min.features = 0)}}

       	# add adt counts
        if(adt){
            message("Step 3 : Adding ADT counts")
            prot <- mtx$`Antibody Capture`
            if(hto){
                prot <- prot[-c(which(str_detect(rownames(prot), hto_str))),]}
            seu_obj[["ADT"]] <- CreateAssay5Object(as.matrix(prot[,colnames(seu_obj[["RNA"]]$counts)]), min.cells = 0, min.features = 0)}

       	# add sample id in metadata
        message("Final Step : Creating Seurat Object")
        seu_obj <- set_assay_keys(seu_obj)
        seu_obj@meta.data$samples <- rep(samples[i], ncol(seu_obj))
        seu_obj <- RenameCells(seu_obj, new.names = paste0(samples[i], "_", colnames(seu_obj)))
        x[[i]] <- seu_obj}

    names(x) <- paste0(samples)
    return(x)}

########### QC FUNCTIONS ###########
## calculate fractions of mitochondrial, ribosomal, haemoglobin, TCR, BCR and MHC-II genes
calculate_fractions <- function(x){
    x <- PercentageFeatureSet(x, pattern = "^[Mm][Tt]-", col.name = "pct.mt")
    x <- PercentageFeatureSet(x, pattern = "^R[Pp][SsLl]", col.name = "pct.rb")
    x <- PercentageFeatureSet(x, pattern = "^H[Bb].*-", col.name = "pct.hb")
    x <- PercentageFeatureSet(x, pattern = "^T[Rr][ABCDGabcdg][VDJCvdjc]", col.name = "pct.tcr")
    x <- PercentageFeatureSet(x, pattern = "^I[Gg][HKLhkl][VDJCAEMGvdjcaemg]", col.name = "pct.bcr")
    x <- PercentageFeatureSet(x, pattern =  "^HLA-", col.name = "pct.mhc")
    return(x)}

## demultiplex cell hashtags
htodemux <- function(x, cluster = F, add.one = F){
    if("HTO" %in% names(x@assays)){
        DefaultAssay(x) <- "HTO"
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) + 1}
        x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR", margin = 2)
        x <- HTODemux(x, assay = "HTO", positive.quantile = 0.99)
        if(cluster){
            x <- htocluster(x)}
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) - 1}
        DefaultAssay(x) <- "RNA"}
    return(x)}

multiseqdemux <- function(x, cluster = F, add.one = F){
    if("HTO" %in% names(x@assays)){
        DefaultAssay(x) <- "HTO"
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) + 1}
        x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR", margin = 2)
        x <- MULTIseqDemux(x, assay = "HTO")
        x@meta.data$MULTI.global <- ifelse(x@meta.data$MULTI_ID == "Doublet", "Doublet", ifelse(x@meta.data$MULTI_ID == "Negative", "Negative", "Singlet"))
        if(cluster){
            x <- htocluster(x)}
        if(add.one){
            x@assays$HTO$counts <- as.matrix(x@assays$HTO$counts) - 1}
        DefaultAssay(x) <- "RNA"}
    return(x)}

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

## calculate mad scores
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

qc_report <- function(x, column, samples){

    for(i in unique(x@meta.data[[samples]])){
        data <- x@meta.data[which(x@meta.data[[samples]] == i),]
        n = length(which(data[[column]] == "Pass"))
        pct = n*100/nrow(data)
        pct = round(pct, 1)
        message(paste0(pct, "% (", n, ") of cells remains - ", i))}
        }

## run doubletfinder
run_doubletfinder <- function(x, assay = "RNA", dims = 1:10, truth = NULL, cluster, dbr, ncores = 1){

    x@meta.data <- x@meta.data %>%
        mutate(cell_barcode = rownames(.))
    obj <- x

    if(!"data" %in% Layers(obj, assay = assay)){
        obj <- NormalizeData(obj, verbose = F)}

    if(length(truth) == 0){
        message("Step 1 : No ground truth")
        suppressMessages({
            sweep_res <- paramSweep(obj, PCs = dims, sct = FALSE, num.cores = ncores)
            sweep_stat <- summarizeSweep(sweep_res, GT = FALSE)})}

    else if (truth %in% colnames(obj@meta.data)){
        message("Step 1 : Filter singlet and doublet from ground truth")
        cells <- colnames(obj)[which(obj@meta.data[[truth]] %in% c("Singlet", "Doublet"))]
        obj <- subset(obj, cells = cells)

        suppressMessages({
            sweep_res <- paramSweep(obj, PCs = dims, sct = FALSE, num.cores = ncores)
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
    prop <- modelHomotypic(x@meta.data[[cluster]])
    nExp_poi <- round(dbr*nrow(x@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-prop))
    obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = top_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    obj@meta.data %>% select(starts_with("pANN_0.25")) %>% colnames() -> pann
    obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = top_pK, nExp = nExp_poi.adj, reuse.pANN = pann, sct = FALSE)
    class <- obj@meta.data %>%
        select(c("cell_barcode", tail(colnames(obj@meta.data),1)))
    colnames(class)[2] <- "DoubletFinder"

    x@meta.data <- x@meta.data %>%
	as.data.frame(.) %>%
        merge(., class, by = "cell_barcode", all.x = T) %>%
        column_to_rownames("cell_barcode") %>%
        as.data.frame(.)
    return(x)
}

## run scdblfinder
run_scdblfinder <- function(x, clusters, samples, dbr, truth = NULL, ncores = 1){
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

## calculate cell cycle
calculate_cellcycle <- function(x, org = "human", remove_genes = T, orig.assay = "RNA"){

    if(!org %in% c("human", "mouse")){
        stop('organism must be either "human" or "mouse"')}
    sgenes <- cc.genes$s.genes
    g2mgenes <- cc.genes$g2m.genes
    if(org == "mouse"){
        sgenes <- convert_human_to_mouse(cc.genes$s.genes)
        g2mgenes <- convert_human_to_mouse(cc.genes$g2m.genes)}
    
    x <- CellCycleScoring(x, s.features = sgenes, g2m.features = g2mgenes, set.ident = F, assay = orig.assay)
    
    if(remove_genes){
        x <- remove_genes(x, features = c(sgenes, g2mgenes), orig.assay = orig.assay, new.assay = "CC")}
    return(x)    
    }

## remove genes from assay
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

## return genes from assay
return_genes <- function(x, from_assay, to_assay = "RNA"){
    counts <- rbind(x[[from_assay]]$counts, x[[to_assay]]$counts)
    data <- rbind(x[[from_assay]]$data, x[[to_assay]]$data)
    x[[to_assay]] <- CreateAssay5Object(counts = counts, data = data, min.cells = 0, min.features = 0)
    x[[from_assay]] <- NULL
    return(x)}

## remove vdj genes from assay
remove_vdj_genes <- function(x, bcr = T, tcr = T, orig.assay = "RNA"){
    if(bcr){
        bcr.genes <- rownames(x[[orig.assay]])[which(str_detect(rownames(x[[orig.assay]]), bcr.string))]
        x <- remove_genes(x, features = bcr.genes, orig.assay = orig.assay, new.assay = "BCR")}
    if(tcr){
        tcr.genes <- rownames(x[[orig.assay]])[which(str_detect(rownames(x[[orig.assay]]), tcr.string))]
        x <- remove_genes(x, features = tcr.genes, orig.assay = orig.assay, new.assay = "TCR")}
    return(x)
        }

## Clustering
process_seurat <- function(x, assay = "RNA", dims = 1:10, reduction = NULL, harmony = F, group.by.vars = "samples"){
    DefaultAssay(x) <- assay
    if(!"data" %in% Layers(x, assay = assay)){
        x <- NormalizeData(x, verbose = F)}
    if(length(reduction) == 0){
        x <- FindVariableFeatures(x)
        x <- ScaleData(x, verbose = F)
        x <- RunPCA(x, verbose = F)
        reduction <- "pca"}
    if(harmony){
        stopifnot(length(group.by.vars) == 1)
        x <- RunHarmony(x, reduction.use = reduction, group.by.vars = group.by.vars, verbose = F)
        reduction <- "harmony"}
    print(ElbowPlot(x, reduction = reduction))
    x <- RunUMAP(x, dims = dims, reduction = reduction, verbose = F)
    x <- FindNeighbors(x, dims = dims, reduction = reduction, verbose = F)
    x <- FindClusters(x, resolution = 0.4, verbose = F)
    print(scUMAP(x, reduction = "umap", group.by = "seurat_clusters", cols = kelly))
    if(harmony){
        print(scUMAP(x, reduction = "umap", group.by = group.by.vars, cols = kelly))}
    return(x)}

qc_before_recluster <- function(x, assay = "RNA"){
    x[[assay]] <- CreateAssay5Object(counts = x[[assay]]$counts, data = x[[assay]]$data, min.features = 200, min.cells = 3)
    x <- subset(x, cells = colnames(x[[assay]]$counts))
    return(x)}

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

calculate_cluster_similarity <- function(x, cluster, assay = "RNA", slot = "data", variable.features = F){
    if(variable.features){
        selected.features <- VariableFeatures(x)}
    else{
	selected.features <- rownames(x[[assay]])}

    pseudobulk <- AverageExpression(x, assay = assay, features = selected.features, group.by = cluster, slot = slot)
    pseudobulk <- pseudobulk[[assay]]
    correlation.matrix <- calculate_correlation(pseudobulk)

    return(correlation.matrix)}


plot_similarity_heatmap <- function(correlation.matrix, annotations = colnames(correlation.matrix), color = NULL){


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
        cluster_rows = T,
        cluster_columns = T,
        heatmap_legend_param = list(
            title = "Pearson\nCorrelation",
            legend_direction = "vertical")
        )
    return(ht)}

## Differential expression
select_top_deg <- function(df, n = 5, rank_by = "avg_log2FC", only.pos = F){

    if(!rank_by %in% c("avg_log2FC", "p_val_adj")){
        stop('please make sure rank_by is either "avg_log2FC" or "p_val_adj"')}

    genelist <- markers %>%
        filter(!str_detect(gene, paste0("^Mt|^Rp[sl]|", bcr.string))) %>%
        distinct(gene, .keep_all = T) %>%
        filter(p_val_adj < 0.05)

    if(only.pos){
        genelist <- genelist %>%
            filter(avg_log2FC > 0 & pct.1 > 0.1)}

    if(rank_by == "avg_log2FC"){
        genelist <- genelist %>%
            group_by(cluster) %>%
            slice_max(n=n, order_by = avg_log2FC, with_ties = F) %>%
            split(.$cluster) %>%
            lapply(function(df) df$gene)}

    if(rank_by == "p_val_adj"){
        genelist <- genelist %>%
            group_by(cluster) %>%
            slice_min(n=n, order_by = p_val_adj, with_ties = F) %>%
            split(.$cluster) %>%
            lapply(function(df) df$gene)}

    return(genelist)}

## Troubleshoot
set_assay_keys <- function(x){
    for(assay in names(x@assays)){
        x[[assay]]@key <- paste0(tolower(assay), "_")}
    return(x)}


## seurat to adata
convert_seurat_to_anndata <- function(x, h5ad, columns = NULL, pca = "pca", umap = "umap", assay = "RNA", return_genes = T){

    # downgrade to v4
    options(Seurat.object.assay.version = "v3")
    if(!is.object(x)){
        if(file.exists(x)){
            x <- qread(paste0(x))}}

    if(!length(columns) > 0){
        message(paste0("selecting all columns from metadata"))
        columns <- colnames(x@meta.data)}
    else{
	message(paste0("selecting columns from metadata: ", paste0(columns, collapse = ", ")))
        x@meta.data <- x@meta.data[,columns]}
    for(i in seq_along(colnames(x@meta.data))){
        x@meta.data[[i]] <- as.character(x@meta.data[[i]])}

    if(return_genes & !str_detect(assay, "SCT")){
        return_assays <- intersect(c("CC", "BCR", "TCR", "MHC"), names(x@assays))
        if(length(return_assays) > 0){
            message(paste0("returning assays (", paste0(return_assays, collapse = ", "), ") to ", assay, " assay"))
            for(i in return_assays){
                x <- return_genes(x, from_assay = i, to_assay = assay)}}}

    message("storing PCA, UMAP to the appropriate reduction slot")
    DefaultAssay(x) <- assay
    x@reductions[["pca"]] <- x@reductions[[pca]]
    x@reductions[["umap"]] <- x@reductions[[umap]]
    x@reductions <- x@reductions[c("pca", "umap")]

    if(any(!str_detect(names(x@assays), "SCT"))){
        v5assays <- names(x@assays)[which(!str_detect(names(x@assays), "SCT"))]
        message(paste0("converting V5 assays (", paste0(v5assays, collapse = ", "), ") to V3 assays"))
        for(i in v5assays){
            x[[i]] <- as(object = x[[i]], Class = "Assay")}}

    h5seurat <- gsub("\\.h5ad", ".h5Seurat", h5ad)
    SaveH5Seurat(x, filename = h5seurat)
    Convert(h5seurat, dest=h5ad)

    options(Seurat.object.assay.version = "v5")
    }

## convert seurat to singlecellexperiment
convert_seurat_to_sce <- function(x, assay = "RNA"){
    x <- DietSeurat(x, assay = assay)
    x[[assay]] <- as(object = x[[assay]], Class = "Assay")
    sce <- as.SingleCellExperiment(x)
    return(sce)}
