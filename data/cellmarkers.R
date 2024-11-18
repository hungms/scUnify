## All Immune Cells
immune.mm <- list(
    Immune = c("Ptprc"),
    `T-cell` =  c("Cd3d", "Cd4", "Cd8a", "Cd8b"),
    NK = c("Nkg7", "Gnly"),
    `B-cell` = c("Cd19", "Cd79a", "Cd79b", "Ms4a1"),
    Plasma = c("Jchain", "Sdc1", "Irf4", "Prdm1"),
    Dendritic = c("Cst3", "Itgax"),
    Monocytes = c("Lyz2", "Itgam", "Cd14"),
    Cycling = c("Mki67", "Pcna", "Top2a"))

immune.hs <- list(
    Immune = c("PTPRC"),
    `T-cell` =  c("CD3D", "CD4", "CD8A", "CD8B"),
    NK = c("NKG7", "GNLY"),
    `B-cell` = c("CD19", "CD79A", "CD79B", "MS4A1"),
    Plasma = c("JCHAIN", "SDC1", "IRF4", "PRDM1"),
    Dendritic = c("CST3", "ITGAX"),
    Monocytes = c("LYZ2", "ITGAM", "CD14"),
    Cycling = c("MKI67", "PCNA", "TOP2A"))

## B Cells
b.mm <- list(
    "B-cells" = c("Cd19", "Cd79a", "Ms4a1"),
    "Activated" = c("Fos", "Cd69", "Dusp2"),
    "Naive" = c("Ccr7", "Sell", "Txnip", "Bank1"),
    "Memory" = c("Klf2", "Cd38", "Hhex", "Gpr183", "Ccr6"),
    "GC-LZ" = c("Fcer2a", "Cd40", "Il4i1", "Myc"),
    "GC-DZ" = c("Bcl6", "Aicda", "Mef2b"),
    "Cycling" = c("Mki67", "Pcna", "Top2a"),
    "Plasma-cells" = c("Jchain", "Irf4", "Prdm1", "Sdc1"))

b.hs <- list(
    "B-cells" = c("CD19", "CD79A", "MS4A1"),
    "Activated" = c("FOS", "CD69", "DUSP2"),
    "Naive" = c("CCR7", "SELL", "TXNIP", "BANK1"),
    "Memory" = c("KLF2", "CD38", "HHEX", "GPR183", "CCR6"),
    "GC-LZ" = c("FCER2A", "CD40", "IL4I1", "MYC"),
    "GC-DZ" = c("BCL6", "AICDA", "MEF2B"),
    "Cycling" = c("MKI67", "PCNA", "TOP2A"),
    "Plasma-cells" = c("JCHAIN", "SDC1", "IRF4", "PRDM1"))

b.complete.mm <- list(
    "B-cells" = c("Cd19", "Cd79a", "Ms4a1", "Cd79b"),
    "Activated" = c("Fos", "Cd69", "Dusp2", "Junb"), 
    "Unswitched" = c("Ighd", "Ighm", "Ighg1", "Igha"),
    "Naive" = c("Ccr7", "Sell", "Txnip", "Bank1", "Fcmr", "Bank1",  "Cd22", "Cd24", "Cd27"),
    "Memory" = c("Klf2", "Cd38", "Hhex", "Tnfrsf13b", "Celf2", "Crip1", "Nt5e", "Itgam", "Fcrl5", "Fcrl4", "Txnip", "Bank1", "Gpr183", "Ccr6"), 
    "Act-Mem" = c("Lrrk2", "Blvrb", "Cd72", "Satb1"),
    "GC-Mem" = c("Cd80", "Pdcd1lg2", "Myadm", "Vim", "Pdlim1"),
    "GC-LZ" = c("Fcer2a", "Cd40", "Il4i1", "Myc", "Bcl2a1a", "Cd83"),
    "GC-DZ" = c("Bcl6", "Aicda", "Mef2b", "Ccnd3", "Cxcr4", "Foxo1"), 
    "Cycling" = c("Mki67", "Pcna", "Top2a"),
    "Plasma-cells" = c("Jchain", "Irf4", "Prdm1", "Sdc1"))

b.complete.hs <- list(
    "B-cells" = c("CD19", "CD79A", "MS4A1", "CD79B"),
    "Activated" = c("FOS", "CD69", "DUSP2", "JUNB"), 
    "Unswitched" = c("IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1"),
    "Naive" = c("CCR7", "SELL", "TXNIP", "BANK1", "FCMR", "BANK1", "CD22", "CD24", "CD27"),
    "Memory" = c("KLF2", "CD38", "HHEX", "TNFRSF13B", "CD38", "CELF2", "CRIP1", "NT5E", "ITGAM", "FCRL5", "FCRL4", "TXNIP", "BANK1", "GPR183", "CCR6"), 
    "Act-Mem" = c("LRRK2", "BLVRB", "CD72", "SATB1"),
    "GC-Mem" = c("CD80", "PDCD1LG2", "MYADM", "VIM", "PDLIM1"),
    "GC-LZ" = c("FCER2", "CD40", "IL4I1", "MYC", "BCL2A1", "CD83"),
    "GC-DZ" = c("BCL6", "AICDA", "MEF2B", "CCND3", "CXCR4", "FOXO1"), 
    "Cycling" = c("MKI67", "PCNA", "TOP2A"),
    "Plasma-cells" = c("JCHAIN", "IRF4", "PRDM1", "SDC1"))

## Plasma Cells
pc.complete.mm <- list(
    PC_Identity = c("Cd19", "Ptprc", "Ms4a1", "Cd38", "Jchain", "Cxcr4", "Cxcl12"),
    Differentiation = c("Irf4", "Prdm1", "Sdc1", "Stat3", "Lyn", "Ikzf3"),
    Survival = c("Tnfsf13", "Tnfrsf13b", "Tnfrsf13c", "Tnfrsf17"),
    AntiApoptosis = c("Mcl1", "Bcl2", "Bcl2l1", "Bcl2a1d", "Cd27", "Zfp667", "Birc3", "Cflar", "Raf1"),
    Antibody = c("Slc3a2", "Slc1a5", "Ell2", "Eaf2"), 
    `Short-Lived` = c("Epcam", "Tigit", "Lag3", "Cxcr3"),
    AX_LLPC = c("Zbtb20", "Fbxw7", "Foxo1", "Dock2", "Cd93", "Il6st", "Ikzf1", "Tcf12", "Creb3l2"), 
    BLIMP1 = c("Cd28", "Vav1", "Grb2", "Plcg1", "Cdkn2c", "Ccnd2", "Cdkn1b"),
    Autophagy = c("Atg2a", "Atg5", "Atg10", "Gabarap", "Gabarapl1", "Gabarapl2", "Deptor", "Rps6", "Map1lc3a", "Sqstm1"),
    Kurosaki = c("Itgb7", "Klf2", "S1pr1"),
    Metabolism = c("Atp12a", "Srf", "Enpp1"),
    ERstress = c("Xbp1", "Atf4", "Atf6", "Manf", "Ern1", "Nos2"),
    Interaction = c("Lepr", "Vcam1", "Icam1", "Tnfrsf13b", "Slc3a2", "Cd44", "Cd37", "Itgal", "Itgb2l", "Itgb1", "Itga4"),
    IgA = c("Tgfb1", "Ccr9", "Ccr10"))

pc.complete.hs <- list(
    PC_Identity = c("CD19", "PTPRC", "MS4A1", "CD38", "JCHAIN", "CXCR4", "CXCL12"),
    Differentiation = c("IRF4", "PRDM1", "SDC1", "STAT3", "LYN", "IKZF3"),
    Survival = c("TNFSF13", "TNFRSF13B", "TNFRSF13C", "TNFRSF17"),
    AntiApoptosis = c("MCL1", "BCL2", "BCL2L1", "BCL2A1", "CD27", "ZNF667", "BIRC3", "CFLAR", "RAF1"),
    Antibody = c("SLC3A2", "SLC1A5", "ELL2", "EAF2"), 
    `Short-Lived` = c("EPCAM", "TIGIT", "LAG3", "CXCR3"),
    AX_LLPC = c("ZBTB20", "FBXW7", "FOXO1", "DOCK2", "CD93", "IL6ST", "IKZF1", "TCF12", "CREB3L2"), 
    BLIMP1 = c("CD28", "VAV1", "GRB2", "PLCG1", "CDKN2C", "CCND2", "CCNE", "CDKN1B"),
    Autophagy = c("ATG2A", "ATG5", "ATG10", "GABARAP", "GABARAPL1", "GABARAPL2", "DEPTOR", "RPS6", "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "SQSTM1"),
    Kurosaki = c("ITGB7", "KLF2", "S1PR1"),
    Metabolism = c("ATP12A", "TFBS", "SRF", "ENPP1"),
    ERstress = c("XBP1", "ATF4", "ATF6", "ERO1LB", "MANF", "ERN1", "NOS2"),
    Interaction = c("LEPR", "VCAM1", "ICAM1", "TNFRSF13B", "SLC3A2", "CD44", "CD37", "ITGAL", "ITGB2", "ITGB1", "ITGA4"),
    IgA = c("TGFB1", "CCR9", "CCR10"))

