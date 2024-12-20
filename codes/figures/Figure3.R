source("utils.R")
dir.create(paste0(dir,"/output/figure"), showWarnings = FALSE)


# Figure 3

## B/plasma
gex_b <- suppressMessages(readRDS(file = paste0(dir,"/output/SeuratObj_B_merged_GEX.rds")))
gex_b <- subset(gex_b, cells = doubletcells, invert = TRUE)
Idents(gex_b) = dplyr::left_join(data.frame(cluster = gex_b@meta.data$seurat_clusters),
                                 cluster_df[cluster_df$cell == "B/plasma",] %>%
                                   dplyr::mutate(cluster = factor(cluster)),
                                 by = "cluster") %>%
  dplyr::mutate(cluster = paste0("B/plasma-",cluster,":",clu_name),
                cluster = factor(cluster, levels = paste0("B/plasma-",cluster_df[cluster_df$cell=="B/plasma",]$cluster,":",cluster_df[cluster_df$cell=="B/plasma",]$clu_name))) %>%
  .$cluster
gex_b@meta.data$cluster = Idents(gex_b)

png(paste0(dir,"/output/figure/UMAP_B.png"), width = 6, height = 6, units = "in", res = 300)
CellDimPlot(
  srt = gex_b, group.by = c("seurat_clusters"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 0.5, pt.alpha = 0.8, 
) +
  scale_color_manual(labels = paste0(cluster_df[cluster_df$cell=="B/plasma",]$clu_name),
                     values = manual_colors,
                     name = "") + 
  guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                               title = "",
                               ncol = 1)
  )
dev.off()

markers = c(
  "CD19", "MS4A1", 
  "CD27","IGHD", 
  "IGHG1","IGHA1","SDC1", 
  "RPN2","XBP1","PRDX4", "MZB1", "CD38", 
  "ITGAX", 
  "TBX21", 
  "NR4A1", 
  "NFKBIA",
  "SLAMF7", 
  "MX1",
  "AICDA", 
  "CD86",  
  "CD5", 
  "BCL2", 
  "CXCR5", 
  "FAS", 
  "TLR9",
  "FCRL4", 
  "FCRL5", 
  "CD83", 
  "CD79B", 
  "CD40", 
  "MME", 
  "PAX5", 
  "ISG15",  
  "MKI67", 
  "PRDM1", 
  "HLA-DRB1",
  "ITGAM", 
  "MT-CO1"
)


pdf(paste0(dir,"/output/figure/AveExpression_B.pdf"), width = 20, height = 5)
plot_aveExp(gex_b, 
            markers,
            fontsize = 18,
            fontsize_row = 18,
            cellwidth = 18, cellheight = 18)
dev.off()

for(gene in markers){
  p = Plot_Density_Custom(gex_b, pt.size = 1, 
                          features = gene) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    guides(colour = "none")
  png(paste0(dir,"/output/figure/Expression_",gene,"_B.png"), width = 2.3, height = 2.5, units = "in", res = 300)
  plot(p)
  dev.off() 
}


data_type = "BCR"
contig_list <- list()
for (sample_id in sample_ids) {
  tryCatch({
    S <- read.csv(paste0(dir,"/data/",sample_id,"_",data_type,"/outs/filtered_contig_annotations.csv"))
    contig_list[[sample_id]] <- S
  }, error = function(e) {
    print(paste("Error reading CSV for sample_id =", sample_id))
  })
}
combined <- combineBCR(contig_list, 
                       samples = sample_ids, 
                       ID = sample_ids)
for(feature in colnames(meta)){
  combined <- addVariable(combined, 
                          name = feature, 
                          variables = meta[,feature])
}

seurat_b <- combineExpression(combined, 
                              gex_b, 
                              cloneCall = "strict", 
                              chain = "both",
                              group.by = "sample", 
                              proportion = FALSE, 
                              cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
pdf(paste0(dir,"/output/figure/ClonalSpaceHomeostasis_UMAP_B.pdf"), width = 4, height = 2)
plot(ggplot() +
       geom_point(
         data = cbind(seurat_b[["umap"]]@cell.embeddings,
                      seurat_b@meta.data) %>% 
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

b_clone = read.table(paste0(dir,"/output/BCR_contig_Frequency_Sequence_merged.txt"),sep="\t",header=TRUE)

top_clone_n = 20

for (disease_group in c("nonSS","pSS(SSA)","centSS","pSS(SSA/cent)")){
  disease_label = gsub(" ","_",gsub("\\/","_",gsub("\\)","_",gsub("\\(","_",disease_group))))
  top_clone = b_clone %>%
    dplyr::filter(Disease == disease_group) %>%
    dplyr::arrange(-Frequency.of.clonotypes) %>%
    dplyr::distinct(combination.of.the.nucleotide.and.gene.sequence.CTstrict., .keep_all=TRUE) %>%
    head(top_clone_n) %>%
    .$combination.of.the.nucleotide.and.gene.sequence.CTstrict.
  pdf(paste0(dir,"/output/figure/ExpandedClones_B_",disease_label,".pdf"), width = 3, height = 2)
  plot(
    b_clone %>%  
      dplyr::distinct(barcode, .keep_all = TRUE) %>%
      dplyr::filter(combination.of.the.nucleotide.and.gene.sequence.CTstrict. %in% top_clone) %>%
      dplyr::group_by(combination.of.the.nucleotide.and.gene.sequence.CTstrict.) %>%
      dplyr::mutate(Frequency = dplyr::n()) %>%
      ungroup() %>%
      merge(.,cluster_df[cluster_df$cell=="B/plasma",],by.x="clusters_B",by.y="cluster") %>%
      dplyr::arrange(-Frequency) %>%
      dplyr::mutate(value = 1) %>%
      ggplot(., aes(x = reorder(combination.of.the.nucleotide.and.gene.sequence.CTstrict.,-Frequency), 
                    y = value, 
                    fill = factor(clusters_B))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name ="cluster", 
                        labels = paste0("B-",cluster_df[cluster_df$cell=="B/plasma",]$cluster,":",cluster_df[cluster_df$cell=="B/plasma",]$clu_name),
                        values=manual_colors) +
      ggtitle(paste0("Cell clusters in expanded clones [", disease_group,"]")) +
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

seurat_b@meta.data$cluster = dplyr::left_join(seurat_b@meta.data,
                                              cluster_df[cluster_df$cell=="B/plasma",] %>%
                                                dplyr::mutate(seurat_clusters = cluster) %>%
                                                dplyr::select(-cluster),
                                              by="seurat_clusters") %>%
  dplyr::mutate(cluster = paste0(cell,"-",cluster,":",clu_name)) %>%
  .$cluster
cluster_colors = manual_colors[1:length(unique(seurat_b@meta.data$cluster))]
names(cluster_colors) = paste0("B/plasma-",names(cluster_colors),":",cluster_df[cluster_df$cell=="B/plasma",]$clu_name)

Idents(seurat_b) = dplyr::left_join(data.frame(cluster = seurat_b@meta.data$seurat_clusters),
                                    cluster_df[cluster_df$cell == "B/plasma",] %>%
                                      dplyr::mutate(cluster = factor(cluster)),
                                    by = "cluster") %>%
  dplyr::mutate(cluster = paste0("B/plasma-",cluster,":",clu_name)) %>%
  .$cluster
Idents(seurat_b) = factor(Idents(seurat_b),
                          levels = cluster_df[cluster_df$cell == "B/plasma",] %>%
                            dplyr::mutate(cluster = factor(cluster)) %>%
                            dplyr::mutate(cluster = paste0("B/plasma-",cluster,":",clu_name)) %>%
                            .$cluster
)
seurat_b@meta.data$cluster = Idents(seurat_b)
clu_COLORS <- manual_colors[1:length(levels(Idents(seurat_b)))]
names(clu_COLORS) = levels(Idents(seurat_b))
circles <- getCirclize(seurat_b,
                       group.by = "cluster")
#Graphing the chord diagram
pdf(paste0(dir,"/output/figure/B_Circlize.pdf"), width = 4, height = 4)
chordDiagram(circles, self.link = 1, grid.col = clu_COLORS)
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
db_airr = dplyr::left_join(db_airr,
                           cellranger_airr[,c("sequence",setdiff(colnames(cellranger_airr),colnames(db_airr)))],
                           by="sequence")
meta_cell = cbind(
  gex_b@meta.data[,c("seurat_clusters",colnames(gex_b@meta.data)[grepl("score",colnames(gex_b@meta.data))])],
  as.data.frame(t(gex_b@assays$RNA@scale.data[all_markers,]))) %>%
  dplyr::mutate(cell_id = stringr::str_split(rownames(.),pattern="_",simplify=TRUE)[,2])
db_airr = dplyr::left_join(db_airr,meta_cell,by = "cell_id")
# remove cells with multiple heavy chain
multi_heavy <- table(filter(db_airr, locus == "IGH")$cell_id)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]
db_airr <- filter(db_airr, !cell_id %in% multi_heavy_cells)
cat(paste("There are", nrow(db_airr), "rows in the data after filtering out cells with multiple heavy chains.\n"))


pdf(paste0(dir,"/output/figure/Diversity_B.pdf"), width = 2.5, height = 1.7)
alphaDiversity_custom("Disease", "B/plasma")
dev.off()

g = alphaDiversity_custom_fixed_list("cluster", "B/plasma", remove_clusters = NA)
pdf(paste0(dir,"/output/figure/Diversity_fixed_B_cluster.pdf"), width = 10, height = 1.7)
plot(g[[1]] + g[[3]] + g[[5]])
dev.off()

g = alphaDiversity_custom_fixed_list("Disease", "B/plasma", remove_clusters = NA)
pdf(paste0(dir,"/output/figure/Diversity_fixed_B_Disease.pdf"), width = 8, height = 1.7)
plot(g[[1]] + g[[3]] + g[[5]])
dev.off()

sample_curve <- alphaDiversity(db_airr[db_airr$Disease %in% c("nonSS","pSS(SSA)","centSS","pSS(SSA/cent)"),] %>%
                                 dplyr::left_join(.,
                                                  cluster_df[cluster_df$cell=="B/plasma",] %>%
                                                    dplyr::mutate(seurat_clusters = cluster),
                                                  by="seurat_clusters") %>%
                                 dplyr::mutate(cluster = paste0(clu_name)) %>%
                                 dplyr::mutate(Disease_cluster = paste0(Disease,"_",cluster)), 
                               group=c("Disease_cluster"), 
                               clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)

for(Q in c(0,1,2,3,4)){
  pdf(paste0(dir,"/output/figure/Diversity_fixed_DiseaseCluster_B_",Q,".pdf"), width = 5, height = 4)
  plotDiversityTest_nested(sample_curve, as.integer(Q), 
                           legend_title="Disease",
                           main_title = paste0("Diversity [q=",Q,")"),
                           colors = disease_colors)
  dev.off()
}




dbs = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome")

for(clu in unique(gex_b@meta.data$cluster)){
  gex_sub <- subset(gex_b, cluster == clu)
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
  
  saveRDS(gex_sub, file=paste0(dir,"/output/DEGsEnrichmentPathways_B_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  
}


gsea_df = data.frame()
for(clu in unique(gex_b@meta.data$cluster)){
  message(clu)
  gex_sub = readRDS(paste0(dir,"/output/DEGsEnrichmentPathways_B_",gsub(":","_",gsub("-","_",gsub("/","",clu))),".rds"))
  gsea_df = rbind(gex_sub@tools$GSEA_Disease_wilcox$enrichment %>%
                    dplyr::mutate(cluster = clu),
                  gsea_df)
}

cat(paste0("Total differentially up-regulated pathways: ", length(unique(gsea_df[gsea_df$p.adjust<0.05 & gsea_df$NES>0,]$Description))))
g = ggvenn::ggvenn(list(`nonSS up-pathways` = gsea_df[gsea_df$Groups=="nonSS" & gsea_df$p.adjust<0.05 & gsea_df$NES>0,]$Description,
                        `CENT+ up-pathways` = gsea_df[gsea_df$Groups=="centSS" & gsea_df$p.adjust<0.05 & gsea_df$NES>0,]$Description,
                        `SSA+ up-pathways` = gsea_df[gsea_df$Groups=="pSS(SSA)" & gsea_df$p.adjust<0.05 & gsea_df$NES>0,]$Description,
                        `SSA+CENT+ up-pathways` = gsea_df[gsea_df$Groups=="pSS(SSA/cent)" & gsea_df$p.adjust<0.05 & gsea_df$NES>0,]$Description), 
                   show_elements = FALSE, show_percentage = TRUE, 
                   digits = 1, fill_color = c("#E41A1C", "#FF7F00", "#984EA3", "#FFFF33"), 
                   fill_alpha = 0.5, stroke_color = "white", stroke_alpha = 0.5, 
                   stroke_size = 0, stroke_linetype = "solid", set_name_color = "black", 
                   set_name_size = 4, text_color = "black", text_size = 3.5, label_sep = ",") 
g
pdf(paste0(dir,"/output/figure/GEX_Venn_uppathways_B.pdf"), width = 3.5, height = 2.5)
plot(g)
dev.off()

pdf(paste0(dir,"/output/figure/gsea_selectedPathways_byDisease_B.pdf"), width = 10.3, height = 3.1)
plot(
  ggplot(gsea_df %>%
           dplyr::filter(p.adjust<0.05) %>%
           dplyr::filter(grepl("B/plasma",cluster)) %>%
           dplyr::filter(NES>0) %>%
           dplyr::filter(Description %in% c(
             unique(grep("Interferon alpha", gsea_df$Description,value = TRUE)),
             unique(grep("Interferon gamma", gsea_df$Description,value = TRUE)),
             unique(grep("Toll Like Receptor 7", gsea_df$Description,value = TRUE)),
             unique(grep("Toll Like Receptor 9", gsea_df$Description,value = TRUE)),
             unique(grep("RIG−I−like", gsea_df$Description,value = TRUE)),
             unique(grep("NOD−like", gsea_df$Description,value = TRUE)),
             unique(grep("NF−kappa B", gsea_df$Description,value = TRUE)),
             unique(grep("TNF signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Antigen processing and presentation", gsea_df$Description,value = TRUE)),
             unique(grep("MAPK signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Ras signaling", gsea_df$Description,value = TRUE)),
             unique(grep("helicase activity", gsea_df$Description,value = TRUE))
           )) %>%
           dplyr::mutate(Description = factor(Description, levels = c(
             unique(grep("Interferon alpha", gsea_df$Description,value = TRUE)),
             unique(grep("Interferon gamma", gsea_df$Description,value = TRUE)),
             unique(grep("Toll Like Receptor 7", gsea_df$Description,value = TRUE)),
             unique(grep("Toll Like Receptor 9", gsea_df$Description,value = TRUE)),
             unique(grep("RIG−I−like", gsea_df$Description,value = TRUE)),
             unique(grep("NOD−like", gsea_df$Description,value = TRUE)),
             unique(grep("NF−kappa B", gsea_df$Description,value = TRUE)),
             unique(grep("TNF signaling pathway 8", gsea_df$Description,value = TRUE)),
             unique(grep("Antigen processing and presentation", gsea_df$Description,value = TRUE)),
             unique(grep("MAPK signaling pathway", gsea_df$Description,value = TRUE)),
             unique(grep("Ras signaling", gsea_df$Description,value = TRUE)),
             unique(grep("helicase activity", gsea_df$Description,value = TRUE))
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



ab_df = read.csv(paste0(dir,"/R_codes/20240803_antibody_response.csv"))

antigen_colors = c("grey90","#984EA3","violet","#FF7F00")
names(antigen_colors) = c("unknown", "Ro60","Ro52","CENT")

pdf(paste0(dir,"/output/figure/antibody_response_B.pdf"), width = 10, height = 3)
ggplot(ab_df, aes(x = sample, y = proportion_cells_in_total, color = cells_in_node, fill = antigen)) +
  geom_bar(stat = "identity", position = "stack", size=0.05) +
  scale_fill_manual(values = antigen_colors) +
  scale_color_manual(values = c("top1" = "black",
                                "top2" = "black",
                                "top3" = "black",
                                "top4" = "black",
                                "top5" = "black"
                                )) +
  facet_wrap(~Disease, scales = "free", nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Antibody response to antigens", 
       y = "Proportion of cells in trees among total B/plasma cells (%)", x = "Sample", fill = "Antigen", color = "Top nodes")
dev.off()
