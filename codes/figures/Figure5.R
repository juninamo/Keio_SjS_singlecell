source("utils.R")
dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)


# Figure 5

ste = readRDS(paste0(dir,"/output/Visium_merged.rds"))

cluster_COLORS <- manual_colors[1:length(unique(Idents(ste)))]
names(cluster_COLORS) = sort(unique(Idents(ste)))

png(paste0(dir,"/output/figure/UMAP_Visium.png"), width = 6, height = 6, units = "in", res = 300)
CellDimPlot(
  srt = ste, group.by = c("seurat_clusters"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5, pt.alpha = 0.8
) +
  
  scale_color_manual(values = cluster_COLORS,
                     name = "") + 
  guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                               title = "",
                               ncol = 1)
  )
dev.off()


png(paste0(dir,"/output/figure/Visium_merged.png"), width = 15, height = 5, units = "in", res = 300)
SpatialDimPlot(ste, 
               label = FALSE, 
               ncol = 5,
               cols = cluster_COLORS )
dev.off()

pdf(paste0(dir,"/output/figure/Visium_merged_vlnplot.pdf"), width = 5, height = 5)
Stacked_VlnPlot(seurat_object = ste, 
                features = c("CD19","CD3E","CD4","CD8A","FCGR3A","EPCAM","THY1"), x_lab_rotate = TRUE,
                colors_use = cluster_COLORS,
                vln_linewidth = 0.01) 
dev.off()



dbs = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome")

clu_marker = data.frame()
for(clu in unique(ste@meta.data$seurat_clusters)){
  message(clu)
  gex_sub <- subset(ste, seurat_clusters == clu)
  gex_sub <- subset(gex_sub, Disease %in% c("pSS(SSA/cent)","centSS","pSS(SSA)","nonSS"))
  clu_marker <- rbind(presto::wilcoxauc(ste@assays$SCT$data, ste@meta.data$Disease) %>% dplyr::mutate(seurat_clusters = clu),
                      clu_marker)
  
}

library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

fgseaRes <- data.frame()
for(term in c("H","C2","C5","C7")){
  msigdbr_df <- msigdbr(species = "human", category = term)
  pathways = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
  for(clu in unique(ste@meta.data$seurat_clusters)){
    message(clu)
    for(disease_group in unique(clu_marker$group)){
      clu_marker_sub = clu_marker %>%
        dplyr::filter(seurat_clusters == clu & group == disease_group) %>%
        dplyr::filter(pval<0.05)
      gene_list <- clu_marker_sub$logFC
      names(gene_list) <- clu_marker_sub$feature
      gene_list_sorted <- sort(gene_list, decreasing = TRUE)
      entrez_ids <- bitr(names(gene_list_sorted), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
      gene_list_entrez <- gene_list_sorted[entrez_ids$SYMBOL]
      names(gene_list_entrez) <- entrez_ids$ENTREZID
      fgseaRes <- rbind(fgsea(pathways = pathways, 
                              stats    = gene_list_entrez,
                              minSize  = 5,
                              maxSize  = 500) %>%
                          dplyr::mutate(Disease = disease_group,
                                        cluster = clu),
                        fgseaRes)
    }
  }
}


df = fgseaRes %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::select(Disease, pathway, NES, cluster) %>%
  tidyr::pivot_wider(names_from = Disease, values_from = NES)  %>%
  na.omit()
df$Disease <- factor(df$cluster)

Heatmap(df[, -c(1,2)], 
        name = "NES", 
        split = df$cluster, 
        row_title = "Pathway", row_names_side = "left",
        row_title_side = "left", row_title_gp = gpar(fontsize = 10),
        column_title = "Disease Specific Clusters", column_title_side = "top",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_row_dend = FALSE)

tmp = fgseaRes %>%
  dplyr::filter(padj<0.05) %>%
  #dplyr::filter(cluster == "4") %>%
  dplyr::filter(NES>0) 
pdf(paste0(dir,"/output/figure/gsea_selectedPathways_byDisease_Visium.pdf"), width = 6, height = 3)
plot(
  ggplot(fgseaRes %>%
           dplyr::filter(padj<0.05) %>%
           dplyr::filter(NES>0) %>%
           dplyr::filter(cluster == "4") %>%
           dplyr::filter(grepl("HALLMARK_",pathway)) %>%
           dplyr::mutate(Disease = factor(Disease, levels = c("nonSS","centSS","pSS(SSA)","pSS(SSA/cent)")),
                         pathway = gsub("HALLMARK_","",pathway,.)) %>%
           dplyr::filter(pathway %in% c("INTERFERON_ALPHA_RESPONSE",
                                        "INTERFERON_GAMMA_RESPONSE",
                                        
                                        "TNFA_SIGNALING_VIA_NFKB",
                                        "IL6_JAK_STAT3_SIGNALING",
                                        "TGF_BETA_SIGNALING",
                                        "EPITHELIAL_MESENCHYMAL_TRANSITION",
                                        "APOPTOSIS",
                                        
                                        "OXIDATIVE_PHOSPHORYLATION"
           )) %>%
           dplyr::mutate(pathway = factor(pathway, levels = c("INTERFERON_ALPHA_RESPONSE",
                                                              "INTERFERON_GAMMA_RESPONSE",
                                                              
                                                              "TNFA_SIGNALING_VIA_NFKB",
                                                              "IL6_JAK_STAT3_SIGNALING",
                                                              "TGF_BETA_SIGNALING",
                                                              "EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                              "APOPTOSIS",
                                                              
                                                              "OXIDATIVE_PHOSPHORYLATION"
           ))),  
         aes(x = Disease, y = pathway, color = NES, size = size)) +
    geom_point(aes(shape = as.logical(padj < 0.05)), stroke = 2) +  # Adjusted p-value for significance
    scale_color_gradient2(low = "#0072B5FF", mid = "white", high = "#BC3C29FF", midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x labels for better visibility
          legend.position = "right") +
    labs(color = "NES", size = "setSize", shape = "Significant (p.adjust < 0.05)") +
    facet_grid(.~cluster, scales = "free_y")
)
dev.off()
