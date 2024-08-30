source("/camp/home/hungm/working/Matthew/library/R/functions/default/genesymbol.R")
bcr.string <- "^I[Gg][HKLhkl][VDJCAEMGvdjcaemg]|^AC233755"
tcr.string <- "^T[Rr][ABCDGabcdg][VDJCvdjc]"

## Convert functions
convert_genedict <- function(genelist, to = "human"){
    gs <- list()
    for(l in seq_along(genelist)){
        if(names(genelist)[l] %in% c("Quality","Cycling")){
            gs[[l]] <- genelist[[l]]}
        else{
            if(to == "human"){
                gs[[l]] <- convert_mouse_to_human(genelist[[l]])}
            if(to == "mouse"){
                gs[[l]] <- convert_human_to_mouse(genelist[[l]])}}}
    names(gs) <- names(genelist)
    return(gs)}

## All Immune Cells
immune.list.mm <- list(
    Quality = c("nFeature_RNA", "nCount_RNA", "pct.mt", "pct.rb", "pct.bcr", "pct.tcr"),
    Immune = c("Ptprc"),
    `T-cell` =  c("Cd3d", "Cd4", "Cd8a", "Cd8b"),
    NK = c("Nkg7", "Gnly"),
    `B-cell` = c("Cd19", "Cd79a", "Cd79b", "Ms4a1"),
    Plasma = c("Jchain", "Sdc1", "Irf4", "Prdm1"),
    Dendritic = c("Cst3", "Itgax"),
    Monocytes = c("Lyz2", "Itgam", "Cd14"),
    Cycling = c("S.Score", "G2M.Score"))


## B Cells
b.list.mm <- list(
    "B-cells" = c("Cd19", "Cd79a", "Ms4a1"),
    "Activated" = c("Fos", "Cd69", "Dusp2"),
    "Naive" = c("Ccr7", "Sell", "Txnip", "Bank1"),
    "Memory" = c("Klf2", "Cd38", "Hhex", "Gpr183", "Ccr6"),
    "GC-LZ" = c("Fcer2a", "Cd40", "Il4i1", "Myc"),
    "GC-DZ" = c("Bcl6", "Aicda", "Mef2b"),
    "Cycling" = c("S.Score", "G2M.Score"),
    "Plasma-cells" = c("Jchain", "Irf4", "Prdm1", "Sdc1"))

b.deep.list.mm <- list(
    "B-cells" = c("Cd19", "Cd79a", "Ms4a1", "Cd79b"),
    "Activated" = c("Fos", "Cd69", "Dusp2", "Junb"), 
    "Unswitched" = c("Ighd", "Ighm", "Ighg1", "Igha"),
    "Naive" = c("Ccr7", "Sell", "Txnip", "Bank1", "Fcmr", "Bank1",  "Cd22", "Cd24", "Cd27"),
    "Memory" = c("Klf2", "Cd38", "Hhex", "Tnfrsf13b", "Cd38", "Celf2", "Crip1", "Nt5e", "Itgam", "Fcrl5", "Fcrl4", "Txnip", "Bank1", "Gpr183", "Ccr6"), 
    "Act-Mem" = c("Lrrk2", "Blvrb", "Cd72", "Satb1"),
    "GC-Mem" = c("Cd80", "Pdcd1lg2", "Myadm", "Vim", "Pdlim1"),
    "GC-LZ" = c("Fcer2a", "Cd40", "Il4i1", "Myc", "Bcl2a1a", "Cd83"),
    "GC-DZ" = c("Bcl6", "Aicda", "Mef2b", "Ccnd3", "Cxcr4", "Foxo1"), 
    "Cycling" = c("S.Score", "G2M.Score"),
    "Plasma-cells" = c("Jchain", "Irf4", "Prdm1", "Sdc1"))

pc.list.hs <- list(
    PC_Identity = c("CD19", "PTPRC", "MS4A1", "CD38", "CXCR4", "CXCL12"),
    BCR = c("IGHA1", "IGHA2", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "JCHAIN", "pct.bcr"), 
    Differentiation = c("IRF4", "PRDM1", "SDC1", "STAT3", "LYN", "IKZF3"),
    Survival = c("TNFSF13", "TNFRSF13B", "TNFRSF13C", "TNFRSF17"),
    AntiApoptosis = c("MCL1", "BCL2", "BCL2L1", "BCL2A1", "CD27", "ZNF667", "BIRC3", "CFLAR", "RAF1"),
    Antibody = c("SLC3A2", "SLC1A5", "ELL2", "EAF2"),
    AX_SLPC = c("EPCAM", "TIGIT", "LAG3", "CXCR3"),
    AX_LLPC = c("ZBTB20", "FBXW7", "FOXO1", "DOCK2", "CD93", "IL6ST", "IKZF1", "TCF12", "CREB3L2", "IGF1R"),
    BLIMP1 = c("CD28", "VAV1", "GRB2", "PLCG1", "CDKN2C", "CCND2", "CCNE", "CDKN1B"),
    Autophagy = c("ATG2A", "ATG5", "ATG10", "GABARAPL1", "DEPTOR"),
    Kurosaki = c("ITGB7", "KLF2", "S1PR1"),
    Metabolism = c("ATP12A", "TFBS", "SRF", "ENPP1"),
    ERstress = c("XBP1", "ATF4", "ATF6", "ERO1LB", "MANF", "ERN1", "NOS2"),
    Interaction = c("LEPR", "VCAM1", "ICAM1", "TNFRSF13B", "SLC3A2", "CD44", "CD37", "ITGAL", "ITGB2", "ITGB1", "ITGA4"),
    IgA = c("TGFB1", "CCR9", "CCR10"))

## Convert all to human
for(i in c("immune.list.mm", "b.list.mm", "b.deep.list.mm")){
    assign(gsub("\\.mm", "\\.hs", i), convert_genedict(get(i)))}

## Convert all to mouse
for(i in c("pc.list.hs")){
    assign(gsub("\\.hs", "\\.mm", i), convert_genedict(get(i), to = "mouse"))}
