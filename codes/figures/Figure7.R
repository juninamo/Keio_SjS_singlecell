source("utils.R")
dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)


# Figure 7

ste = readRDS(paste0(dir,"/output/Visium_merged.rds"))

cluster_COLORS <- manual_colors[1:length(unique(Idents(ste)))]
names(cluster_COLORS) = sort(unique(Idents(ste)))


Stacked_VlnPlot(seurat_object = ste, 
                features = c("CD3E","CD8A","GZMK","GZMB"), x_lab_rotate = TRUE,
                colors_use = cluster_COLORS,
                vln_linewidth = 0.01) 

pdf(paste0(dir,"/output/figure/Visium_GZMKB.pdf"), width =  8, height = 5)
SpatialFeaturePlot(subset(ste, CD3E > 0.5 & CD8A > 0.5), features = c("GZMK","GZMB"))
dev.off()

gzmb_df = data.frame(gzm = ste@assays$SCT$data["GZMB",(ste@assays$SCT$data["CD3E",] > 0.5) & (ste@assays$SCT$data["CD8A",] > 0.5)],
                     t(ste@assays$SCT$data[,(ste@assays$SCT$data["CD3E",] > 0.5) & (ste@assays$SCT$data["CD8A",] > 0.5)])
)
gzmk_df = data.frame(gzm = ste@assays$SCT$data["GZMK",(ste@assays$SCT$data["CD3E",] > 0.5) & (ste@assays$SCT$data["CD8A",] > 0.5)],
                     t(ste@assays$SCT$data[,(ste@assays$SCT$data["CD3E",] > 0.5) & (ste@assays$SCT$data["CD8A",] > 0.5)])
)

calculate_correlations <- function(data) {
  gzm_col <- data$gzm
  other_columns <- setdiff(names(data), "gzm")
  correlations <- sapply(other_columns, function(col) cor(gzm_col, data[[col]], use = "complete.obs", method = "spearman"))
  return(correlations)
}

gzmb_correlations <- sort(na.omit(calculate_correlations(gzmb_df)), decreasing = TRUE)
gzmb_correlations

gzmk_correlations <- sort(na.omit(calculate_correlations(gzmk_df)), decreasing = TRUE)
gzmk_correlations

genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", names(gzmb_correlations), value = TRUE)
genes_use <- setdiff(names(gzmb_correlations), genes_exclude)
gzmb_correlations = gzmb_correlations[genes_use]
genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", names(gzmk_correlations), value = TRUE)
genes_use <- setdiff(names(gzmk_correlations), genes_exclude)
gzmk_correlations = gzmk_correlations[genes_use]

gzmb_correlations_df = data.frame(cor = gzmb_correlations,
                                  gene = names(gzmb_correlations)) %>%
  dplyr::filter(gene != "GZMB")
limit_b <- max(abs(gzmb_correlations_df$cor)) * c(-1, 1)
gzmk_correlations_df = data.frame(cor = gzmk_correlations,
                                  gene = names(gzmk_correlations)) %>%
  dplyr::filter(gene != "GZMK")
limit_k <- max(abs(gzmk_correlations_df$cor)) * c(-1, 1)
g_b = ggplot(rbind(head(gzmb_correlations_df,10),
                   tail(gzmb_correlations_df,10)) %>%
               dplyr::arrange(cor), 
             aes(reorder(gene, cor), cor, fill = cor)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_distiller(type = "div", limit = limit_b)
g_k = ggplot(rbind(head(gzmk_correlations_df,10),
                   tail(gzmk_correlations_df,10)) %>%
               dplyr::arrange(cor), 
             aes(reorder(gene, cor), cor, fill = cor)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_distiller(type = "div", limit = limit_k)

pdf(paste0(dir,"/output/figure/Visium_corgenes.pdf"), width = 7, height = 3)
plot(g_b + g_k)
dev.off()

library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

fgseaRes <- data.frame()
for(clu in c("GZMK","GZMB")){
  if(clu == "GZMB"){
    gene_list = gzmb_correlations
  } else if(clu == "GZMK"){
    gene_list = gzmk_correlations
  }
  gene_list_sorted <- sort(gene_list, decreasing = TRUE)
  entrez_ids <- bitr(names(gene_list_sorted), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  gene_list_entrez <- gene_list_sorted[entrez_ids$SYMBOL]
  names(gene_list_entrez) <- entrez_ids$ENTREZID
  for(term in c("H","C2","C5","C7")){
    msigdbr_df <- msigdbr(species = "human", category = term)
    pathways = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
    fgseaRes <- rbind(fgsea(pathways = pathways, 
                            stats    = gene_list_entrez,
                            minSize  = 5,
                            maxSize  = 500) %>%
                        dplyr::mutate(cluster = clu,
                                      term = term),
                      fgseaRes)
  }
}

tmp_k = fgseaRes %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::filter(NES>0)  %>%
  dplyr::filter(cluster == "GZMK")
tmp_b = fgseaRes %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::filter(NES>0)  %>%
  dplyr::filter(cluster == "GZMB")
g = ggvenn::ggvenn(list(`GZMK pathways` = tmp_k$pathway,
                        `GZMB pathways` = tmp_b$pathway), 
                   show_elements = FALSE, show_percentage = TRUE, 
                   digits = 1, fill_color = c("yellow", "green"), 
                   fill_alpha = 0.5, stroke_color = "white", stroke_alpha = 0.5, 
                   stroke_size = 0, stroke_linetype = "solid", set_name_color = "black", 
                   set_name_size = 4, text_color = "black", text_size = 3.5, label_sep = ",") 
g
pdf(paste0(dir,"/output/figure/Visium_Venn_GZMKB.pdf"), width = 3, height = 2)
plot(g)
dev.off()

pathways = c("GOBP_FIBROBLAST_APOPTOTIC_PROCESS",
             "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
             
             "GOBP_RESPONSE_TO_INTERLEUKIN_1",                                                         
             "GOBP_RESPONSE_TO_INTERLEUKIN_2",  
             "GOBP_RESPONSE_TO_INTERLEUKIN_4",
             "GOBP_RESPONSE_TO_INTERLEUKIN_6",
             "GOBP_RESPONSE_TO_INTERLEUKIN_15",
             
             "GOBP_INTERFERON_BETA_PRODUCTION",
             "NATSUME_RESPONSE_TO_INTERFERON_BETA_UP",
             "WP_TYPE_III_INTERFERON_SIGNALING" ,
             
             "GOBP_ACUTE_INFLAMMATORY_RESPONSE",                                                      
             "GOBP_CHRONIC_INFLAMMATORY_RESPONSE",                                              
             
             "GOBP_T_CELL_MEDIATED_CYTOTOXICITY"        )

pdf(paste0(dir,"/output/figure/Visium_GSEA_GZMKB.pdf"), width = 6.3, height = 3.5)
plot(
  ggplot(fgseaRes %>%
           dplyr::filter(padj<0.05) %>%
           dplyr::filter(NES>0) %>%
           dplyr::filter(pathway %in% pathways) %>%
           dplyr::mutate(pathway = factor(pathway, levels = pathways)),  
         aes(x = cluster, y = pathway, color = NES, size = size)) +
    geom_point(aes(shape = as.logical(padj < 0.05)), stroke = 2) + 
    scale_color_gradient2(low = "#0072B5FF", mid = "white", high = "#BC3C29FF", midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  
          legend.position = "right") +
    labs(color = "NES", size = "setSize", shape = "Significant (p.adjust < 0.05)")
)
dev.off()



