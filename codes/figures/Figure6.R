source("utils.R")
dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)


# Figure 6

ste = readRDS(paste0(dir,"/output/Visium_merged.rds"))
cluster_COLORS <- manual_colors[1:length(unique(Idents(ste)))]
names(cluster_COLORS) = sort(unique(Idents(ste)))

n_factors = 4
outfile = paste0(dir,"/output/Visium_ST_model_merged_gpu_GroupSlide_factor",n_factors,".hdf5")
MOFAobject <- load_model(outfile)
plot_top_weights(MOFAobject, factors = "all", nfeatures = 20, sign = "positive", abs=FALSE) + 
  facet_wrap(~factor, nrow = 2, scales = "free")
for(i in 1:n_factors){
  pdf(paste0(dir,"/output/figure/Visium_merged_weights_Factor",i,".pdf"), width = 2.5, height = 3.5)
  plot(plot_top_weights(MOFAobject, factors = i, nfeatures = 20, sign = "positive", abs=FALSE) + 
         facet_wrap(~factor, nrow = 2, scales = "free") +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

Z <- get_factors(MOFAobject, factors = "all", as.data.frame = TRUE)
Z <- Z[complete.cases(Z), ]
df <- get_covariates(MOFAobject, as.data.frame = TRUE) %>% 
  merge(Z, by = "sample", suffixes = c(".covariate", ".factor"))
covariates_dt <- df %>%
  tidyr::pivot_wider(names_from="covariate", values_from="value.covariate") 
covariates_dt <- mutate(covariates_dt, color_by = value.factor) # for compatibility with .add_legend
covariates.names <- c(colnames(covariates_dt)[ncol(covariates_dt)-1], colnames(covariates_dt)[ncol(covariates_dt)])
covariates_dt %<>%
  dplyr::mutate(slide = stringr::str_split(sample,pattern="_",simplify=TRUE) %>% .[,1]) %>%
  dplyr::mutate(
    imagerow = dplyr::case_when(
      slide == "slide1" ~ imagerow,
      slide == "slide2" ~ imagerow - 20000,
      slide == "slide3" ~ imagerow - 40000,
      slide == "slide4" ~ imagerow - 60000,
      slide == "slide5" ~ imagerow - 80000
    ),
    imagecol = dplyr::case_when(
      slide == "slide1" ~ imagecol,
      slide == "slide2" ~ imagecol - 20000,
      slide == "slide3" ~ imagecol - 40000,
      slide == "slide4" ~ imagecol - 60000,
      slide == "slide5" ~ imagecol - 80000
    )) %>%
  dplyr::left_join(.,
                   ste@meta.data %>%
                     tibble::rownames_to_column("sample") %>%
                     dplyr::select(-slide),
                   by="sample")

ggplot(covariates_dt, aes(x=.data$imagecol,
                          y=.data$imagerow * -1,
                          col = .data$value.factor)) +
  geom_point(size = 1, shape=21, alpha = 1) +
  scale_fill_gradient2(low = "#0072B5FF", mid="gray90", high = "#BC3C29FF") + 
  scale_color_gradient2(low = "#0072B5FF", mid="gray90", high = "#BC3C29FF") + 
  facet_grid(factor ~ slide) + coord_fixed() +
  theme_bw() +
  theme(
    axis.text = element_text(size = rel(0.9), color = "black"),
    axis.title = element_text(size = rel(1.0), color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  ) + guides(col = guide_colorbar(title = "Factor value"))


covariates_dt <- df %>%
  tidyr::pivot_wider(names_from="covariate", values_from="value.covariate") 
covariates_dt <- mutate(covariates_dt, color_by = value.factor) # for compatibility with .add_legend
covariates.names <- c(colnames(covariates_dt)[ncol(covariates_dt)-1], colnames(covariates_dt)[ncol(covariates_dt)])
covariates_dt %<>%
  dplyr::left_join(.,
                   ste@meta.data %>%
                     tibble::rownames_to_column("sample") %>%
                     dplyr::select(-slide),
                   by="sample") %>%
  dplyr::mutate(
    Disease = factor(Disease, levels = c("nonSS","pSS(SSA)","centSS","pSS(SSA/cent)")))

for(disease_group in c("nonSS",
                       "pSS(SSA)",
                       "pSS(SSA/cent)",
                       "centSS")){
  disease_label = gsub(" ","_",gsub("\\/","_",gsub("\\)","_",gsub("\\(","_",disease_group))))
  g_list <- vector("list", length(unique(covariates_dt[covariates_dt$Disease == disease_group, ]$sample_id)))
  names(g_list) = unique(covariates_dt[covariates_dt$Disease == disease_group, ]$sample_id)
  size = 1
  for (clu in unique(covariates_dt[covariates_dt$Disease == disease_group, ]$sample_id)){
    g_list[[clu]] = ggplot(covariates_dt[covariates_dt$Disease == disease_group, ] %>%
                             dplyr::filter(sample_id == clu), 
                           aes(x=.data$imagecol,
                               y=.data$imagerow * -1,
                               col = .data$value.factor)) +
      geom_point(size = size, shape=16, alpha = 0.8) +
      scale_fill_gradient2(low = "#0072B5FF", mid="gray90", high = "#BC3C29FF") + 
      scale_color_gradient2(low = "#0072B5FF", mid="gray90", high = "#BC3C29FF") + 
      facet_grid(factor ~ Disease) + 
      coord_fixed() +
      theme_void() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
      ) + 
      guides(col = guide_colorbar(title = "Factor value"))
  }
  g0 = patchwork::wrap_plots(g_list, nrow = 1)
  pdf(paste0(dir,"/output/figure/factors_Visium_",disease_label,".pdf"), width = 15, height = 10)
  plot(g0)
  dev.off()
}


factors.dt <- Reduce(rbind,get_factors(MOFAobject, factors="all")) %>% 
  data.table::as.data.table(keep.rownames = T) %>% 
  data.table::setnames("rn","sample")
factors.dt = dplyr::left_join(ste@reductions$umap@cell.embeddings %>%
                                as.data.frame() %>%
                                tibble::rownames_to_column("sample"),
                              factors.dt,
                              by="sample")
all(factors.dt$sample == rownames(ste@meta.data))
data_to_plot <- dplyr::left_join(ste@meta.data %>%
                                   tibble::rownames_to_column("sample"),
                                 factors.dt,
                                 by="sample") %>%
  dplyr::select(seurat_clusters, sample, starts_with("Factor")) %>%
  tidyr::pivot_longer(cols = -c(seurat_clusters,sample), names_to = "factor", values_to = "factor_value") %>%
  dplyr::mutate(factor = gsub("factor_value.","",factor)) 

png(paste0(dir,"/output/figure/Visium_merged_factors_byClusters.png"), width = 7, height = 5, units = "in", res = 300)
ggplot(data_to_plot %>%
         dplyr::filter(factor == "Factor1"), 
       aes(x = factor, y = factor_value, color = seurat_clusters)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 0.5, alpha = 0.3) +
  geom_boxplot(position = position_dodge(width = 0.8),
               alpha = 0.5, outlier.shape = NA,
               width=0.6) +
  scale_color_manual(values = cluster_COLORS) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(color="black", size=20)
  ) +
  labs(
    x = "Factor",
    y = "Factor value",
    color = "Cluster",
    title = "Factor values per cell cluster"
  ) +
  guides(color = guide_legend(ncol = 5, override.aes = list(size = 3, alpha = 1)))
dev.off()

data_to_plot <- dplyr::left_join(ste@meta.data %>%
                                   tibble::rownames_to_column("sample"),
                                 factors.dt,
                                 by="sample") %>%
  dplyr::select(Disease, sample, starts_with("Factor")) %>%
  tidyr::pivot_longer(cols = -c(Disease,sample), names_to = "factor", values_to = "factor_value") %>%
  dplyr::mutate(factor = gsub("factor_value.","",factor)) 
ggplot(data_to_plot, aes(x = factor, y = factor_value, color = Disease)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
             size = 0.5, alpha = 0.3) +
  geom_boxplot(position = position_dodge(width = 0.5),
               alpha = 0.5, outlier.shape = NA,
               width=0.4) +
  scale_color_manual(values = c(
    "nonSS" = "#E41A1C",
    "pSS(SSA)" = "#377EB8",
    "centSS" = "#4DAF4A",
    "pSS(SSA/cent)" = "#984EA3"
  )) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(color="black", size=20)
  ) +
  labs(
    x = "Factor",
    y = "Factor value",
    color = "Cluster",
    title = "Factor values per cluster"
  ) +
  guides(color = guide_legend(ncol = 5, override.aes = list(size = 3, alpha = 1)))


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
data_to_plot <- dplyr::left_join(ste@meta.data %>%
                                   tibble::rownames_to_column("sample"),
                                 factors.dt,
                                 by="sample") %>%
  dplyr::mutate(lymph_score = T.score + B.plasma.score ) %>%
  dplyr::select(Disease, sample, lymph_score, starts_with("Factor")) %>%
  tidyr::pivot_longer(cols = -c(Disease,sample,lymph_score), names_to = "factor", values_to = "factor_value") %>%
  dplyr::mutate(factor = gsub("factor_value.","",factor))  %>%
  dplyr::mutate(Disease = factor(Disease,levels = c("nonSS", 
                                                    "pSS(SSA)",
                                                    "centSS",
                                                    "pSS(SSA/cent)"
  ))) %>%
  dplyr::group_by(factor) %>%
  dplyr::mutate(density = get_density(lymph_score, factor_value)) %>%
  dplyr::ungroup() 
pdf(paste0(dir,"/output/figure/Factor_values_per_lymphocyte_score.pdf"), width = 3.5, height = 3)
ggplot(data_to_plot, aes(x = lymph_score, y = factor_value, color = density)) +
  geom_point_rast(size = 0.01, 
                  raster.dpi=300) +
  scale_color_viridis() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~factor, scales = "free") +
  theme_bw(base_size = 8) +
  labs(
    x = "Lymphocyte score",
    y = "Factor value",
    color = "Cluster",
    title = "Factor values per lymphocyte score"
  ) 
dev.off()


dplyr::left_join(ste@meta.data %>%
                   tibble::rownames_to_column("sample"),
                 factors.dt,
                 by="sample") %>%
  dplyr::select(Disease, sample_id, GreenSpan, starts_with("Factor")) %>%
  tidyr::pivot_longer(cols = -c(Disease,sample_id, GreenSpan), names_to = "factor", values_to = "factor_value") %>%
  dplyr::mutate(factor = gsub("factor_value.","",factor))  %>%
  dplyr::mutate(Disease = factor(Disease,levels = c("nonSS", 
                                                    "pSS(SSA)",
                                                    "centSS",
                                                    "pSS(SSA/cent)"
  ))) %>%
  dplyr::mutate(factor_positive = ifelse(factor_value > 0, "positive", "negative")) %>%
  dplyr::group_by(factor_positive,sample_id,Disease,factor,GreenSpan) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::group_by(sample_id,factor) %>%
  dplyr::mutate(total = sum(count)) %>%
  dplyr::arrange(sample_id) %>%
  dplyr::mutate(proportion = 100*count/total) %>%
  dplyr::filter(factor_positive == "positive") %>%
  ggscatter(data = ., x = "GreenSpan", y = "proportion", 
            size = 1.5, alpha = 0.6,
            add = "reg.line") +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2) +
  theme_minimal() +
  facet_wrap(~factor, scales = "free") 

W <- get_weights(MOFAobject)[["single_view"]]
n_features = 20
fac="Factor1"
top_weights = rownames(W)[order(W[,fac], decreasing = TRUE)][1:n_features]
print(W[top_weights,fac] )

gex_s <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_SG_merged_GEX.rds")))
gex_s <- subset(gex_s, cells = doubletcells, invert = TRUE)
Idents(gex_s) = dplyr::left_join(data.frame(cluster = gex_s@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "SG",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = paste0("SG-",cluster,":",clu_name)) %>%
  .$cluster
Idents(gex_s) = factor(Idents(gex_s),
                       levels = cluster_df[cluster_df$cell == "SG",] %>%
                         dplyr::mutate(cluster = factor(cluster)) %>%
                         dplyr::mutate(cluster = paste0("SG-",cluster,":",clu_name)) %>%
                         .$cluster
)
cluster_COLORS <- manual_colors[1:length(unique(Idents(gex_s)))]
names(cluster_COLORS) = sort(unique(Idents(gex_s)))
pdf(paste0(dir,"/output/figure/SG_Factor1.pdf"), width = 4, height = 4.3)
plot(
  Stacked_VlnPlot(seurat_object = gex_s, 
                  features = top_weights[1:10], x_lab_rotate = TRUE,
                  colors_use = cluster_COLORS,
                  vln_linewidth = 0.1)
)
dev.off()


pathway_df_all=data.frame()
for(slide_id in 1:5){
  sce <- zellkonverter::readH5AD(paste0(dir,"/output/COMMOT_slide",slide_id,"_withoutObsp.h5ad"))
  pathway_df = cbind(reducedDims(sce)[["commot.user_database.sum.sender"]],
                     reducedDims(sce)[["commot.user_database.sum.receiver"]]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_barcodes") %>%
    dplyr::mutate(cell_barcodes = paste0("slide",slide_id,"_",cell_barcodes)) %>%
    dplyr::left_join(.,ste@meta.data,by="cell_barcodes") %>%
    na.omit()
  pathway_df_all = rbind(pathway_df_all,pathway_df)
}
head(pathway_df_all[,1:5])

Z <- get_factors(MOFAobject, factors = "all", as.data.frame = TRUE)
Z <- Z[complete.cases(Z), ]
Z %<>%
  dplyr::select(-group) %>%
  tidyr::pivot_wider(names_from="factor", values_from="value")

pathway_df_all = dplyr::left_join(pathway_df_all,Z,by=c("cell_barcodes"="sample"))
head(pathway_df_all[,1:5])
head(pathway_df_all[,(ncol(pathway_df_all)-4):ncol(pathway_df_all)])

pathways = sort(gsub("s-","",colnames(reducedDims(sce)[["commot.user_database.sum.sender"]])))
message(paste0("Total number of pathways: ",length(pathways)))


cor_res = data.frame()
for(pathwayName in pathways){
  df = data.frame(sender = pathway_df_all[,paste0("s-",pathwayName)],
                  receiver = pathway_df_all[,paste0("r-",pathwayName)],
                  Factor1 = pathway_df_all[,paste0("Factor1")]) 
  cor_res = rbind(cor_res,
                  data.frame(pathway = pathwayName,
                             cor_sender = cor(df$sender, df$receiver, method = "spearman"),
                             pval_sender = cor.test(df$sender, df$receiver, method = "spearman")$p.value,
                             cor_receiver = cor(df$receiver, df$Factor1, method = "spearman"),
                             pval_receiver = cor.test(df$receiver, df$Factor1, method = "spearman")$p.value
                  ))
}
cor_res %<>%
  dplyr::mutate(padj_sender = p.adjust(pval_sender, method = "BH"),
                padj_receiver = p.adjust(pval_receiver, method = "BH")) 
message(paste0("Number of significant pathways: ",sum(cor_res$padj_sender < 0.05 & 
                                                        abs(cor_res$cor_sender)>0.3 &
                                                        cor_res$padj_receiver < 0.05 & 
                                                        abs(cor_res$cor_receiver)>0.3, na.rm = TRUE)))
cor_res %>%
  dplyr::filter(padj_sender < 0.05 & abs(cor_sender)>0.3 &
                  padj_receiver < 0.05 & abs(cor_receiver)>0.3) %>%
  write.csv(paste0(dir,"/output/COMMOT_Pathway_Correlation.csv"), row.names = FALSE)

cor_res %>%
  filter(padj_sender < 0.05, padj_receiver < 0.05) %>%
  as.data.frame() %>%
  ggplot() +
  geom_step(aes(x = sort(cor_sender), y = cumsum(sort(cor_sender, decreasing = TRUE))), color = "blue") +
  geom_step(aes(x = sort(cor_receiver), y = cumsum(sort(cor_receiver, decreasing = TRUE))), color = "red") +
  labs(x = "Correlation Coefficient", y = "Cumulative Correlation", title = "Cumulative Correlation Curves")


pathwayName="CXCL14-CXCR4"

Sender = stringr::str_split(pathwayName, pattern = "-", simplify = TRUE)[1]
Receiver = stringr::str_split(pathwayName, pattern = "-", simplify = TRUE)[1]

data.frame(sender = pathway_df_all[,paste0("s-",pathwayName)],
           receiver = pathway_df_all[,paste0("r-",pathwayName)],
           Factor1 = pathway_df_all[,paste0("Factor1")]) %>%
  ggplot(aes(x = sender, y = Factor1)) +
  geom_point_rast(size = 0.01, 
                  raster.dpi=300) +
  scale_color_viridis() +
  stat_cor(method = "spearman", cor.coef.name = "rho", p.method = "spearman", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw(base_size = 8) +
  labs(
    x = paste0("Sender signal [",Sender,"]"),
    y = "Factor1",
    color = "Density",
    title = "Sender vs Factor1"
  )
data.frame(sender = pathway_df_all[,paste0("s-",pathwayName)],
           receiver = pathway_df_all[,paste0("r-",pathwayName)],
           Factor1 = pathway_df_all[,paste0("Factor1")]) %>%
  ggplot(aes(x = receiver, y = Factor1)) +
  geom_point_rast(size = 0.01, 
                  raster.dpi=300) +
  scale_color_viridis() +
  stat_cor(method = "spearman", cor.coef.name = "rho", p.method = "spearman", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw(base_size = 8) +
  labs(
    x = paste0("Receiver signal [",Receiver,"]"),
    y = "Factor1",
    color = "Density",
    title = "Receiver vs Factor1"
  )




