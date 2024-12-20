source("utils.R")

dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)

# Figure 1
gex <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_merged_GEX.rds")))
gex <- subset(gex, cells = doubletcells, invert = TRUE)
Idents(gex) = dplyr::left_join(data.frame(cluster = gex@meta.data$seurat_clusters),
                               cluster_df[cluster_df$cell == "ALL",] %>%
                                 dplyr::mutate(cluster = factor(cluster)),
                               by = "cluster") %>%
  dplyr::mutate(cluster = clu_name) %>%
  .$cluster

pdf(paste0(dir,"/output/figure/UMAP_all.pdf"), width = 6, height = 6)
CellDimPlot(
  srt = gex, group.by = c("seurat_clusters"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 2
) +
  scale_color_manual(labels = paste0(cluster_df[cluster_df$cell=="ALL",]$clu_name),
                     values = manual_colors,
                     name = "")
dev.off()

cluster_COLORS <- manual_colors[1:length(unique(Idents(gex)))]
names(cluster_COLORS) = unique(Idents(gex))
pdf(paste0(dir,"/output/figure/Violin_all.pdf"), width = 4, height = 6)
Stacked_VlnPlot(seurat_object = gex, 
                features = c("PTPRC",
                             "CD19","MS4A1", "CD79A", "SDC1",
                             "CD3E", 
                             "FCGR3A","CD68", 
                              "THY1", "VCAM1", 
                             "EPCAM",
                             "MUC7", 
                             "MUC5B",
                             "TFCP2L1", 
                             "KRT5","CAV1" 
                             ), 
                vln_linewidth = 0.1,
                             x_lab_rotate = TRUE,
                colors_use = cluster_COLORS) 
dev.off()

pdf(paste0(dir,"/output/figure/Dotplot_all.pdf"), width = 5, height = 3)
DotPlot_scCustom(gex, features = c("CD79A", 
                                   "CD3E", 
                                   "JCHAIN",
                                   "CD68", 
                                   "MUC7",
                                   "LYZ"), 
                 flip_axes = T, x_lab_rotate = TRUE)
dev.off()
for (gene in c("CD79A", 
               "CD3E", 
               "JCHAIN",
               "CD68", 
               "MUC7", 
               "LYZ")){
  p1 <- Plot_Density_Custom(gex, pt.size = 0.1, features = gene)
  png(paste0(dir,"/output/figure/Expression_",gene,"_all.png"), width = 3, height = 2.5, units = "in", res = 300)
  plot(p1)
  dev.off()
}

# Extended Figure 1
gex_cellranger = readRDS(paste0(dir,"/output/SeuratObj_merged_GEX_cellranger.rds"))
Idents(gex_cellranger) = dplyr::left_join(data.frame(cluster = gex_cellranger@meta.data$seurat_clusters),
                               cluster_df[cluster_df$cell == "ALL_cellranger",] %>%
                                 dplyr::mutate(cluster = factor(cluster)),
                               by = "cluster") %>%
  dplyr::mutate(cluster = clu_name) %>%
  .$cluster
pdf(paste0(dir,"/output/figure/Dotplot_all_woCellBender.pdf"), width = 5, height = 3)
DotPlot_scCustom(gex_cellranger, features = c("CD79A", 
                                   "CD3E", 
                                   "JCHAIN",
                                   "CD68", 
                                   "MUC7",
                                   "LYZ"), 
                 flip_axes = T, x_lab_rotate = TRUE)
dev.off()
for (gene in c("CD79A", 
               "CD3E", 
               "JCHAIN",
               "CD68", 
               "MUC7",
               "LYZ")){
  p1 <- Plot_Density_Custom(gex_cellranger, pt.size = 0.1, features = gene)
  png(paste0(dir,"/output/figure/Expression_",gene,"_all_woCellBender.png"), width = 3, height = 2.5, units = "in", res = 300)
  plot(p1)
  dev.off()
}


LISI_df = rbind(
  readRDS(file=paste0(dir,"/output/LISI_harmony_ALL.rds")) %>%
    dplyr::mutate(group = "after"),
  readRDS(file=paste0(dir,"/output/LISI_PCA_ALL.rds")) %>%
    dplyr::mutate(group = "before")
) 

LISI_df %>%
  dplyr::group_by(key,group) %>%
  dplyr::summarize(median = median(val))

pdf(paste0(dir,"/output/figure/LISI_all.pdf"), width = 10, height = 3)
plot(
  ggplot(LISI_df) +
    geom_density(
      aes(x  = val, color = group) 
    ) +
    labs(
      title = "effect of batch correction",
      x = "LISI score",
      y = "Density"
    ) +
    facet_wrap( ~ key, scales = "free", ncol = 6) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.text.x=element_text(size=15, color="black"),
          strip.text.y=element_text(size=15, color="black"),
          legend.position = "bottom",
          plot.title = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size =15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          legend.text =  element_text(size = 15),
          legend.key.size = grid::unit(0.5, "lines"),
          legend.title = element_text(size = 0.8, hjust = 0)) 
)
dev.off()
