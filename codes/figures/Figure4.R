source("utils.R")
dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)


# Figure 4

## SG
gex_s <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_SG_merged_GEX.rds")))
gex_s <- subset(gex_s, cells = doubletcells, invert = TRUE)
Idents(gex_s) = dplyr::left_join(data.frame(cluster = gex_s@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "SG",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = paste0("SG-",cluster,":",clu_name),
                cluster = factor(cluster, levels = paste0("SG-",cluster_df[cluster_df$cell=="SG",]$cluster,":",cluster_df[cluster_df$cell=="SG",]$clu_name))) %>%
  .$cluster

png(paste0(dir,"/output/figure/UMAP_SG.png"), width = 6, height = 6, units = "in", res = 300)
CellDimPlot(
  srt = gex_s, group.by = c("seurat_clusters"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5, pt.alpha = 0.8, 
) +
  scale_color_manual(labels = paste0(cluster_df[cluster_df$cell=="SG",]$clu_name),
                     values = manual_colors,
                     name = "") + 
  guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                               title = "",
                               ncol = 1)
  )
dev.off()

markers = c(
  "PDPN", "THY1", "ICAM1","VCAM1",
  "EPCAM", 
  "LYZ","PIP","MUC7",
  "MUC5B","CEACAM6","BPIFB2", 
  "TFCP2L1","KRT19", 
  "KRT5","KRT14","TP63","CAV1", 
  "ACTA2", 
  "CDH5", 
  "MCAM", 
  "CCL21","CCL19","TNC", "CCL8", "CCL2", 
  "CD82","CXCL9","CCL19","CD34","TNFSF13B","PDPN", "THY1", "ICAM1","VCAM1","PRG4","FN1" 
)

pdf(paste0(dir,"/output/figure/AveExpression_SG.pdf"), width = 12, height = 5)
plot_aveExp(gex_s, 
            markers,
            fontsize = 18,
            fontsize_row = 18,
            cellwidth = 18, cellheight = 18)
dev.off()

for(gene in unique(markers)){
  p = Plot_Density_Custom(gex_s, pt.size = 1, 
                          features = gene) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    guides(colour = "none")
  png(paste0(dir,"/output/figure/Expression_",gene,"_SG.png"), width = 2.3, height = 2.5, units = "in", res = 300)
  plot(p)
  dev.off() 
}

pdf(paste0(dir,"/output/figure/DotPlot_SG.pdf"), width = 6, height = 6.3)
DotPlot_scCustom(gex_s, features = unique(markers), flip_axes = T, x_lab_rotate = TRUE)
dev.off()

## Other leukocytes
gex_o <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_OtherLym_merged_GEX.rds")))
gex_o <- subset(gex_o, cells = doubletcells, invert = TRUE)
Idents(gex_o) = dplyr::left_join(data.frame(cluster = gex_o@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "other",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = clu_name) %>%
  .$cluster

png(paste0(dir,"/output/figure/UMAP_OL.png"), width = 6, height = 6, units = "in", res = 300)
CellDimPlot(
  srt = gex_o, group.by = c("seurat_clusters"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 2, pt.alpha = 1, 
) +
  scale_color_manual(labels = paste0(cluster_df[cluster_df$cell=="other",]$clu_name),
                     values = manual_colors,
                     name = "") + 
  guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                               title = "",
                               ncol = 1)
  )
dev.off()

markers = c("CD14", "FCGR3A", 
            "CD163", 
            "NCAM1","FCER1G","KLRF1", 
            "NKG7","TYROBP","PRF1",
            "GZMB","GZMH","XCL1","GZMK",
            "CLEC7A",
            "CD1C"
)

pdf(paste0(dir,"/output/figure/AveExpression_OL.pdf"), width = 12, height = 5)
plot_aveExp(gex_o, 
            markers,
            fontsize = 18,
            fontsize_row = 18,
            cellwidth = 18, cellheight = 18)
dev.off()

for(gene in unique(markers)){
  p = Plot_Density_Custom(gex_o, pt.size = 1, 
                          features = gene) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    guides(colour = "none")
  png(paste0(dir,"/output/figure/Expression_",gene,"_OL.png"), width = 2.3, height = 2.5, units = "in", res = 300)
  plot(p)
  dev.off() 
}

pdf(paste0(dir,"/output/figure/DotPlot_OL.pdf"), width = 4.5, height = 4)
DotPlot_scCustom(gex_o, features = unique(markers), flip_axes = T, x_lab_rotate = TRUE)
dev.off()

## Other leukocytes
gex_o <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_OtherLym_merged_GEX.rds")))
gex_o <- subset(gex_o, cells = doubletcells, invert = TRUE)
Idents(gex_o) = dplyr::left_join(data.frame(cluster = gex_o@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "other",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = clu_name) %>%
  .$cluster



dbs = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome")


gex_s <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_SG_merged_GEX.rds")))
gex_s <- subset(gex_s, cells = doubletcells, invert = TRUE)
Idents(gex_s) = dplyr::left_join(data.frame(cluster = gex_s@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "SG",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = paste0("SG-", cluster,":",clu_name)) %>%
  .$cluster
gex_s@meta.data$cluster = Idents(gex_s)
unique(gex_s@meta.data$cluster)
for(clu in unique(gex_s@meta.data$cluster)){
  gex_sub <- subset(gex_s, cluster == clu)
  gex_sub <- subset(gex_sub, Disease %in% c("pSS(SSA/cent)","centSS","pSS(SSA)","nonSS"))

  # DE test
  for(test.use in c("wilcox")){
    tryCatch({
      gex_sub <- RunDEtest(srt = gex_sub, group_by = "Disease", 
                           test.use = test.use,
                           fc.threshold = 1, only.pos = FALSE)
    }, error = function(e) {
      print(paste("Error test.use =", test.use))
    })
  }
  genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(gex_sub@assays$RNA$scale.data), value = TRUE)
  genes_use <- setdiff(row.names(gex_sub@assays$RNA$scale.data), genes_exclude)
  gex_sub@tools$DEtest_Disease$AllMarkers_wilcox = gex_sub@tools$DEtest_Disease$AllMarkers_wilcox %>%
    dplyr::filter(gene %in% genes_use)
  gex_sub <- AnnotateFeatures(gex_sub, species = "Homo_sapiens", db = c("TF", "CSPA"))
  ## Enrichment test
  gex_sub <- RunEnrichment(
    srt = gex_sub, group_by = "Disease", 
    db = dbs,
    species = "Homo_sapiens",
    DE_threshold = "avg_log2FC > log2(1.5) & p_val < 0.05"
  )
  # GSEA test
  gex_sub <- RunGSEA(
    srt = gex_sub, group_by = "Disease", 
    db = dbs, 
    species = "Homo_sapiens",
    DE_threshold = "p_val < 0.05"
  )
  
  saveRDS(gex_sub, file=paste0(dir,"/output/DEGsEnrichmentPathways_SG_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  
}



dbs = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome")

for(clu in c(unique(gex_s@meta.data$cluster),
             unique(gex_o@meta.data$cluster))){
  gex_sub <- subset(gex_s, cluster == clu)
  gex_sub <- subset(gex_sub, Disease %in% c("pSS(SSA/cent)","centSS","pSS(SSA)","nonSS"))
  
  # DE test
  for(test.use in c("wilcox")){
    tryCatch({
      gex_sub <- RunDEtest(srt = gex_sub, group_by = "Disease", 
                           test.use = test.use,
                           fc.threshold = 1, only.pos = FALSE)
    }, error = function(e) {
      print(paste("Error test.use =", test.use))
    })
  }
  genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(gex_sub@assays$RNA$scale.data), value = TRUE)
  genes_use <- setdiff(row.names(gex_sub@assays$RNA$scale.data), genes_exclude)
  gex_sub@tools$DEtest_Disease$AllMarkers_wilcox = gex_sub@tools$DEtest_Disease$AllMarkers_wilcox %>%
    dplyr::filter(gene %in% genes_use)
  gex_sub <- AnnotateFeatures(gex_sub, species = "Homo_sapiens", db = c("TF", "CSPA"))
  ## Enrichment test
  gex_sub <- RunEnrichment(
    srt = gex_sub, group_by = "Disease", 
    db = dbs,
    species = "Homo_sapiens",
    DE_threshold = "avg_log2FC > log2(1.5) & p_val < 0.05"
  )
  # GSEA test
  gex_sub <- RunGSEA(
    srt = gex_sub, group_by = "Disease", 
    db = dbs, 
    species = "Homo_sapiens",
    DE_threshold = "p_val < 0.05"
  )
  
  saveRDS(gex_sub, file=paste0(dir,"/output/DEGsEnrichmentPathways_SG_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  
}


gsea_df = data.frame()
for(clu in c(unique(gex_s@meta.data$cluster),
             unique(gex_o@meta.data$cluster))){
  message(clu)
  gex_sub = readRDS(paste0(dir,"/output/DEGsEnrichmentPathways_SG_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  gsea_df = rbind(gex_sub@tools$GSEA_Disease_wilcox$enrichment %>%
                    dplyr::mutate(cluster = clu),
                  gsea_df)
}


pdf(paste0(dir,"/output/figure/gsea_selectedPathways_byDisease_SGOL.pdf"), width = 16, height = 3.5)
plot(
  ggplot(gsea_df %>%
           dplyr::filter(p.adjust<0.05) %>%
           dplyr::filter(grepl("SG-",cluster) | grepl("OL-",cluster)) %>%
           dplyr::filter(NES>0) %>%
           dplyr::filter(Description %in% c(
             unique(grep("Mitochondrial", gsea_df$Description,value = TRUE)),
             
             unique(grep("defense response to bacterium", gsea_df$Description,value = TRUE)),
             unique(grep("^adaptive immune response$", gsea_df$Description,value = TRUE)),
             unique(grep("Interferon", gsea_df$Description,value = TRUE)),
             unique(grep("Toll Like Receptor", gsea_df$Description,value = TRUE)),
             unique(grep("RIG−I−like", gsea_df$Description,value = TRUE)),
             unique(grep("NOD−like", gsea_df$Description,value = TRUE)),
             unique(grep("NF−kappa B", gsea_df$Description,value = TRUE)),
             unique(grep("TNF signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Extrafollicular and follicular B cell activation", gsea_df$Description,value = TRUE)),
             unique(grep("^cytokine-mediated signaling pathway$", gsea_df$Description,value = TRUE))
           )) %>%
           dplyr::mutate(Disc = factor(Description, levels = c(
             unique(grep("Mitochondrial", gsea_df$Description,value = TRUE)),
             
             unique(grep("defense response to bacterium", gsea_df$Description,value = TRUE)),
             unique(grep("^adaptive immune response$", gsea_df$Description,value = TRUE)),
             unique(grep("Interferon", gsea_df$Description,value = TRUE)),
             unique(grep("Toll Like Receptor", gsea_df$Description,value = TRUE)),
             unique(grep("RIG−I−like", gsea_df$Description,value = TRUE)),
             unique(grep("NOD−like", gsea_df$Description,value = TRUE)),
             unique(grep("NF−kappa B", gsea_df$Description,value = TRUE)),
             unique(grep("TNF signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Extrafollicular and follicular B cell activation", gsea_df$Description,value = TRUE)),
             unique(grep("^cytokine-mediated signaling pathway$", gsea_df$Description,value = TRUE))
           ))) %>%
           dplyr::mutate(Groups = factor(Groups, levels = c("nonSS","centSS","pSS(SSA)","pSS(SSA/cent)"))), 
         aes(x = Groups, y = Description, color = NES, size = setSize)) +
    geom_point(aes(shape = as.logical(p.adjust < 0.05)), stroke = 2) +  # Adjusted p-value for significance
    scale_color_gradient2(low = "#0072B5FF", mid = "white", high = "#BC3C29FF", midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x labels for better visibility
          legend.position = "right") +
    labs(color = "NES", size = "setSize", shape = "Significant (p.adjust < 0.05)") +
    facet_grid(.~cluster, scales = "free_y")
)
dev.off()

