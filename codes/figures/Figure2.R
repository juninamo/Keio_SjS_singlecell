source("utils.R")
dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)

# Figure 2
gex_t <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_T_merged_GEX.rds")))
gex_t <- subset(gex_t, cells = doubletcells, invert = TRUE)
Idents(gex_t) = dplyr::left_join(data.frame(cluster = gex_t@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "T",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = paste0("T-",cluster,":",clu_name),
                cluster = factor(cluster, levels = paste0("T-",cluster_df[cluster_df$cell=="T",]$cluster,":",cluster_df[cluster_df$cell=="T",]$clu_name))) %>%
  .$cluster
gex_t@meta.data$cluster = Idents(gex_t)

cluster_COLORS <- manual_colors[1:length(unique(Idents(gex_t)))]
names(cluster_COLORS) = sort(unique(Idents(gex_t)))

png(paste0(dir,"/output/figure/UMAP_T.png"), width = 6, height = 6, units = "in", res = 300)
CellDimPlot(
  srt = gex_t, group.by = c("cluster"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5, pt.alpha = 0.8
) +
  
  scale_color_manual(values = cluster_COLORS,
                     name = "") + 
  guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                               title = "",
                               ncol = 1)
  )
dev.off()

markers = c(
  "CCR7",
  "CD4", "CD8A", 
  "CXCR3","CXCR4","CXCR5","CXCR6","CX3CR1","CCR2","CCR4","CCR6", 
  "CCL5","CCL4","CXCL13",
  "FOXP3","IL2RA",
  "PDCD1","CTLA4","ICOS","CD40LG", "TIGIT",
  "MAF",
  "RORC", 
  "GATA3", 
  "TBX21", 
  "HLA-DRA","HLA-DRB1", 
  "B3GAT1",
  "CD27", 
  "CD38", 
  "KLRB1", "SLC4A10", 
  "GZMB","GZMK","GZMA","GZMH","IFNG",
  "GNLY","PRF1","NKG7", 
  "MKI67", 
  "ZEB2",
  "MT-CO1"
)

pdf(paste0(dir,"/output/figure/AveExpression_T.pdf"), width = 20, height = 5)
plot_aveExp(gex_t, 
            markers,
            fontsize = 18,
            fontsize_row = 18,
            cellwidth = 18, cellheight = 18)
dev.off()


Stacked_VlnPlot(seurat_object = gex_t, 
                features = markers, 
                vln_linewidth = 0.1,
                x_lab_rotate = TRUE,
                colors_use = cluster_COLORS) 

for(gene in markers){
  p = Plot_Density_Custom(gex_t, pt.size = 0.5, 
                          features = gene) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
          ) +
    guides(colour = "none")
  png(paste0(dir,"/output/figure/Expression_",gene,"_T.png"), width = 2.3, height = 2.5, units = "in", res = 300)
  plot(p)
  dev.off() 
}

pdf(paste0(dir,"/output/figure/DotPlot_T.pdf"), width = 6, height = 8)
DotPlot_scCustom(gex_t, features = unique(markers)[!grepl("^MT-",unique(markers))], flip_axes = T, x_lab_rotate = TRUE)
dev.off()

data_type = "TCR"
contig_list <- list()
for (sample_id in sample_ids) {
  tryCatch({
    S <- read.csv(paste0(dir,"/data/",sample_id,"_",data_type,"/outs/filtered_contig_annotations.csv"))
    contig_list[[sample_id]] <- S
  }, error = function(e) {
    print(paste("Error reading CSV for sample_id =", sample_id))
  })
}
combined <- combineTCR(contig_list, 
                       samples = sample_ids, 
                       ID = sample_ids)
for(feature in colnames(meta)){
  combined <- addVariable(combined, 
                          name = feature, 
                          variables = meta[,feature])
}

seurat_t <- combineExpression(combined, 
                              gex_t, 
                              cloneCall = "strict", 
                              chain = "both",
                              group.by = "sample", 
                              proportion = FALSE, 
                              cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
pdf(paste0(dir,"/output/figure/ClonalSpaceHomeostasis_UMAP_T.pdf"), width = 4, height = 2)
plot(ggplot() +
       geom_point(
         data = cbind(seurat_t[["umap"]]@cell.embeddings,
                      seurat_t@meta.data) %>% 
           dplyr::mutate(alpha = ifelse(is.na(cloneType),0.2,1)) %>%
           dplyr::sample_frac(1L),
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "cloneType", alpha = "alpha"),
         size = 1, stroke = 0, shape = 16
       ) +
       labs(
         title = paste0("UMAP colored by cloneType groups"),
         x = "UMAP1",
         y = "UMAP2"
       ) + 
       scale_color_manual(name="cloneType", values=colorblind_vector(5)) +
       theme_classic(base_size = 10) +
       theme(
         legend.position = "right",
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         panel.grid = element_blank(),
         plot.title = element_text(color="black", size=10, face = "italic")
       ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)),
                  fill = guide_legend(override.aes = list(alpha = 1, size = 2)),
                  alpha="none")
)
dev.off()

t_clone = read.table(paste0(dir,"/output/TCR_contig_Frequency_Sequence_merged.txt"),sep="\t",header=TRUE)

top_clone_n = 20

for (disease_group in c("pSS(SSA)","centSS","pSS(SSA/cent)")){
  disease_label = gsub(" ","_",gsub("\\/","_",gsub("\\)","_",gsub("\\(","_",disease_group))))
  top_clone = t_clone %>%
    dplyr::filter(Disease %in% disease_group) %>%
    dplyr::arrange(-Frequency.of.clonotypes) %>%
    dplyr::distinct(combination.of.the.nucleotide.and.gene.sequence.CTstrict., .keep_all=TRUE) %>%
    head(top_clone_n) %>%
    .$combination.of.the.nucleotide.and.gene.sequence.CTstrict.
  df = t_clone %>%  
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    dplyr::filter(combination.of.the.nucleotide.and.gene.sequence.CTstrict. %in% top_clone) %>%
    merge(.,cluster_df[cluster_df$cell=="T",],by.x="clusters_T",by.y="cluster") %>%
    dplyr::filter(grepl("CD8$",clu_name)) %>%
    dplyr::group_by(combination.of.the.nucleotide.and.gene.sequence.CTstrict.) %>%
    dplyr::mutate(Frequency = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-Frequency) %>%
    dplyr::mutate(value = 1)
  pdf(paste0(dir,"/output/figure/ExpandedClones_CD8T_",disease_label,".pdf"), width = 3, height = 2)
  plot(
    ggplot(df, aes(x = reorder(combination.of.the.nucleotide.and.gene.sequence.CTstrict.,-Frequency), 
                   y = value, 
                   fill = factor(clusters_T))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name ="cluster", 
                        labels = paste0("T-",cluster_df[cluster_df$cell=="T" & grepl("CD8$",cluster_df$clu_name) & cluster_df$cluster %in% unique(df$clusters_T),]$cluster,":",cluster_df[cluster_df$cell=="T" & cluster_df$cluster %in% unique(df$clusters_T) & grepl("CD8$",cluster_df$clu_name),]$clu_name),
                        values=manual_colors) +
      ggtitle(paste0("CD8+ cell clusters in expanded clones [", disease_group,"]")) +
      labs(x = paste0("Top ",top_clone_n, " expanded clone"), 
           y = "Cells in clone") +
      theme_minimal() +
      theme(strip.text.x=element_text(size=10, color="black"),
            strip.text.y=element_text(size=10, color="black"),
            legend.position = "none",
            plot.title = element_text(size=10),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.text =  element_text(size = 10),
            legend.key.size = grid::unit(0.5, "lines"),
            legend.title = element_text(size = 0.8, hjust = 0)) + 
      guides(fill = guide_legend(override.aes = list(size=1,alpha=1),
                                 title = "cluster",
                                 ncol = 1))
  )
  dev.off()
}


for (disease_group in c("pSS(SSA)","centSS","pSS(SSA/cent)")){
  disease_label = gsub(" ","_",gsub("\\/","_",gsub("\\)","_",gsub("\\(","_",disease_group))))
  top_clone = t_clone %>%
    dplyr::filter(Disease %in% disease_group) %>%
    dplyr::arrange(-Frequency.of.clonotypes) %>%
    dplyr::distinct(combination.of.the.nucleotide.and.gene.sequence.CTstrict., .keep_all=TRUE) %>%
    head(top_clone_n) %>%
    .$combination.of.the.nucleotide.and.gene.sequence.CTstrict.
  
  df = t_clone %>%  
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    dplyr::filter(combination.of.the.nucleotide.and.gene.sequence.CTstrict. %in% top_clone) %>%
    merge(.,cluster_df[cluster_df$cell=="T",],by.x="clusters_T",by.y="cluster") %>%
    dplyr::filter(!grepl("CD8$",clu_name) & !grepl("DN$",clu_name) & !grepl("MT-high",clu_name) & !grepl("Proliferating",clu_name)) %>%
    dplyr::group_by(combination.of.the.nucleotide.and.gene.sequence.CTstrict.) %>%
    dplyr::mutate(Frequency = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-Frequency) %>%
    dplyr::mutate(value = 1)
  pdf(paste0(dir,"/output/figure/ExpandedClones_CD4T_",disease_label,".pdf"), width = 3, height = 2)
  plot(
    ggplot(df, aes(x = reorder(combination.of.the.nucleotide.and.gene.sequence.CTstrict.,-Frequency), 
                   y = value, 
                   fill = factor(clusters_T))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name ="cluster", 
                        labels = paste0("T-",cluster_df[cluster_df$cell=="T" & !grepl("CD8$",cluster_df$clu_name) & !grepl("DN$",cluster_df$clu_name) & !grepl("MT-high",cluster_df$clu_name) & !grepl("Proliferating",cluster_df$clu_name) & cluster_df$cluster %in% unique(df$clusters_T),]$cluster,":",cluster_df[cluster_df$cell=="T" & !grepl("CD8$",cluster_df$clu_name) & !grepl("DN$",cluster_df$clu_name) & !grepl("MT-high",cluster_df$clu_name) & !grepl("Proliferating",cluster_df$clu_name) & cluster_df$cluster %in% unique(df$clusters_T),]$clu_name),
                        values=manual_colors) +
      ggtitle(paste0("CD4+ cell clusters in expanded clones [", disease_group,"]")) +
      labs(x = paste0("Top ",top_clone_n, " expanded clone"), 
           y = "Cells in clone") +
      theme_minimal() +
      theme(strip.text.x=element_text(size=10, color="black"),
            strip.text.y=element_text(size=10, color="black"),
            legend.position = "none",
            plot.title = element_text(size=10),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.text.x = element_blank(),
            legend.text =  element_text(size = 10),
            legend.key.size = grid::unit(0.5, "lines"),
            legend.title = element_text(size = 0.8, hjust = 0)) + 
      guides(fill = guide_legend(override.aes = list(size=1,alpha=1),
                                 title = "cluster",
                                 ncol = 1))
  )
  dev.off()
}


seurat_t@meta.data$cluster = dplyr::left_join(seurat_t@meta.data,
                                              cluster_df[cluster_df$cell=="T",] %>%
                                                dplyr::mutate(seurat_clusters = cluster),
                                              by="seurat_clusters") %>%
  dplyr::mutate(cluster = paste0(cell,"-",cluster,":",clu_name)) %>%
  .$cluster
cluster_colors = manual_colors[1:length(unique(seurat_t@meta.data$cluster))]
names(cluster_colors) = paste0("T-",names(cluster_colors),":",cluster_df[cluster_df$cell=="T",]$clu_name)

Idents(seurat_t) = dplyr::left_join(data.frame(cluster = seurat_t@meta.data$seurat_clusters),
                                    cluster_df[cluster_df$cell == "T",] %>%
                                      dplyr::mutate(cluster = factor(cluster)),
                                    by = "cluster") %>%
  dplyr::mutate(cluster = paste0("T-",cluster,":",clu_name)) %>%
  .$cluster
Idents(seurat_t) = factor(Idents(seurat_t),
                          levels = c("T-0:GZMB+GNLY+ CD8","T-4:GZMB+GNLY+ CD8","T-11:GZMA+ CD8","T-1:GZMK/B+ CD8","T-9:GZMK+GZMA+GZMH+ CD8", "T-10:GZMK+GZMA+ CD8",
                                     "T-2:Th1(CXCR3+IFNG+)", "T-3:Tph/Tfh" ,"T-5:CMCD4(CCR7+)", "T-6:Th17(CCR6+RORC+)","T-7:Treg")
)
seurat_t@meta.data$cluster = Idents(seurat_t)
cd8_clu = paste0("T-",cluster_df[cluster_df$cell=="T" & grepl("CD8$",cluster_df$clu_name),]$cluster,":",cluster_df[cluster_df$cell=="T" & grepl("CD8$",cluster_df$clu_name),]$clu_name)

circles <- getCirclize(subset(seurat_t, idents = cd8_clu),
                       group.by = "cluster")
#Graphing the chord diagram
pdf(paste0(dir,"/output/figure/CD8T_Circlize.pdf"), width = 4, height = 4)
chordDiagram(circles, self.link = 1, grid.col = cluster_COLORS)
dev.off()

cd4_clu = c("T-2:Th1(CXCR3+IFNG+)", "T-3:Tph/Tfh" ,"T-5:CMCD4(CCR7+)", "T-6:Th17(CCR6+RORC+)","T-7:Treg")
circles <- getCirclize(subset(seurat_t, idents = cd4_clu),
                       group.by = "cluster")
#Graphing the chord diagram
pdf(paste0(dir,"/output/figure/CD4T_Circlize.pdf"), width = 4, height = 4)
chordDiagram(circles, self.link = 1, grid.col = cluster_COLORS)
dev.off()



db_airr <- do.call("rbind", lapply(sample_ids, function(sample_id) {
  filepath <- paste0(dir, "/annotation/imm/", sample_id, "_", data_type, "/heavy_parse-select_productive_parse-select_clone-pass_germ-pass.tsv")
  if (!file.exists(filepath)) {
    print(paste("File doesn't exist:", filepath))
    return(NULL)
  }
  data <- airr::read_rearrangement(filepath)
  data$sample_id <- sample_id
  return(data)
}))
db_airr = dplyr::left_join(db_airr,
                           meta %>%
                             dplyr::rename(sample_id = sample),
                           by="sample_id") %>% 
  dplyr::filter(productive) %>%
  dplyr::mutate(cell_id = stringr::str_split(sequence_id, pattern="_",simplify = TRUE)[,3])

cellranger_airr <- do.call("rbind", lapply(sample_ids, function(sample_id) {
  filepath <- paste0(dir, "/data/", sample_id, "_", data_type, "/outs/airr_rearrangement.tsv")
  data <- airr::read_rearrangement(filepath)
  data$sample_id <- sample_id
  return(data)
}))
cellranger_airr = cellranger_airr %>%
  dplyr::distinct(sequence, .keep_all = TRUE)
db_airr = dplyr::left_join(db_airr,
                           cellranger_airr[,c("sequence",setdiff(colnames(cellranger_airr),colnames(db_airr)))],
                           by="sequence")
meta_cell = cbind(
  gex_t@meta.data[,c("seurat_clusters",colnames(gex_t@meta.data)[grepl("score",colnames(gex_t@meta.data))])],
  as.data.frame(t(gex_t@assays$RNA@scale.data[all_markers,]))) %>%
  dplyr::mutate(cell_id = stringr::str_split(rownames(.),pattern="_",simplify=TRUE)[,2])
db_airr = dplyr::left_join(db_airr,meta_cell,by = "cell_id")
# remove cells with multiple heavy chain
multi_heavy <- table(filter(db_airr, locus == "TRB" | locus == "TRD")$cell_id)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]

db_airr <- filter(db_airr, !cell_id %in% multi_heavy_cells)
cat(paste("There are", nrow(db_airr), "rows in the data after filtering out cells with multiple heavy chains.\n"))

pdf(paste0(dir,"/output/figure/Diversity_T.pdf"), width = 2.5, height = 1.7)
alphaDiversity_custom("Disease", "T")
dev.off()

pdf(paste0(dir,"/output/figure/Diversity_CD4T.pdf"), width = 2.5, height = 1.7)
alphaDiversity_custom_cd4("Disease", "T")
dev.off()

pdf(paste0(dir,"/output/figure/Diversity_CD8T.pdf"), width = 2.5, height = 1.7)
alphaDiversity_custom_cd8("Disease", "T")
dev.off()


g = alphaDiversity_custom_fixed_list("cluster", "T", remove_clusters = c("T-8:DN","T-12:MT-high (low quality)","T-13:Proliferating"))
pdf(paste0(dir,"/output/figure/Diversity_fixed_T.pdf"), width = 13, height = 1.7)
plot(g[[1]] + g[[3]] + g[[5]])
dev.off()

set.seed(1234)
sample_curve <- alphaDiversity(db_airr[db_airr$Disease %in% c("pSS(SSA)","centSS","pSS(SSA/cent)"),] %>%
                                 dplyr::left_join(.,
                                                  cluster_df[cluster_df$cell=="T",] %>%
                                                    dplyr::mutate(seurat_clusters = cluster),
                                                  by="seurat_clusters") %>%
                                 dplyr::mutate(cluster = paste0(clu_name)) %>%
                                 dplyr::mutate(Disease_cluster = paste0(Disease,"_",cluster)), 
                               group=c("Disease_cluster"), 
                               clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)
for(Q in c(0,1,2,3,4)){
  pdf(paste0(dir,"/output/figure/Diversity_fixed_DiseaseCluster_T_",Q,".pdf"), width = 8, height = 8)
  plotDiversityTest_nested(sample_curve, as.integer(Q), 
                           legend_title="Disease",
                           main_title = paste0("Diversity [q=",Q,")"),
                           colors = disease_colors)
  dev.off()
}


dbs = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome")
for(clu in c("T-0:GZMB+GNLY+ CD8","T-4:GZMB+GNLY+ CD8","T-11:GZMA+ CD8","T-1:GZMK/B+ CD8","T-9:GZMK+GZMA+GZMH+ CD8", "T-10:GZMK+GZMA+ CD8",
             "T-2:Th1(CXCR3+IFNG+)", "T-3:Tph/Tfh" ,"T-5:CMCD4(CCR7+)", "T-6:Th17(CCR6+RORC+)","T-7:Treg")){
  gex_sub <- subset(gex_t, cluster == clu)
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
  
  saveRDS(gex_sub, file=paste0(dir,"/output/DEGsEnrichmentPathways_T_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  
}


gsea_df = data.frame()
for(clu in unique(gex_t@meta.data$cluster)){
  message(clu)
  gex_sub = readRDS(paste0(dir,"/output/DEGsEnrichmentPathways_T_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  gsea_df = rbind(gex_sub@tools$GSEA_Disease_wilcox$enrichment %>%
                    dplyr::mutate(cluster = clu),
                  gsea_df)
}


pdf(paste0(dir,"/output/figure/gsea_selectedPathways_byDisease_T.pdf"), width = 10, height = 3.8)
plot(
  ggplot(gsea_df %>%
           dplyr::filter(p.adjust<0.05) %>%
           dplyr::filter(grepl("T-",cluster)) %>%
           dplyr::filter(NES>0) %>%
           dplyr::filter(Description %in% c(
             unique(grep("inner mitochondrial membrane protein complex", gsea_df$Description,value = TRUE)),
             unique(grep("lysosome", gsea_df$Description,value = TRUE)),
             
             unique(grep("JAK-STAT", gsea_df$Description,value = TRUE)),
             unique(grep("NF-kappa B signaling", gsea_df$Description,value = TRUE)),
             unique(grep("Extrafollicular and follicular B cell activation", gsea_df$Description,value = TRUE)),
             unique(grep("NF−kappa B", gsea_df$Description,value = TRUE)),
             unique(grep("MAPK signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("TNF signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Interferon", gsea_df$Description,value = TRUE)),
             unique(grep("interferon", gsea_df$Description,value = TRUE)),
             unique(grep("^cytokine-mediated signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Epstein-Barr virus infection", gsea_df$Description,value = TRUE))
           )) %>%
           dplyr::mutate(Description = factor(Description, levels = c(
             unique(grep("inner mitochondrial membrane protein complex", gsea_df$Description,value = TRUE)),
             unique(grep("lysosome", gsea_df$Description,value = TRUE)),
             
             unique(grep("JAK-STAT", gsea_df$Description,value = TRUE)),
             unique(grep("NF-kappa B signaling", gsea_df$Description,value = TRUE)),
             unique(grep("Extrafollicular and follicular B cell activation", gsea_df$Description,value = TRUE)),
             unique(grep("NF−kappa B", gsea_df$Description,value = TRUE)),
             unique(grep("MAPK signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("TNF signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Interferon", gsea_df$Description,value = TRUE)),
             unique(grep("interferon", gsea_df$Description,value = TRUE)),
             unique(grep("^cytokine-mediated signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Epstein-Barr virus infection", gsea_df$Description,value = TRUE))
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

