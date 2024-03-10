library(Seurat)
library(dplyr)

# Import counts matrix and set up seurat object
countsN1 <- Read10X_h5("../data/x1_negative_filtered_feature_bc_matrix.h5")
countsP1 <- Read10X_h5("../data/x1_positive_filtered_feature_bc_matrix.h5")
countsN2 <- Read10X_h5("../data/x2_negative_filtered_feature_bc_matrix.h5")
countsP2 <- Read10X_h5("../data/x2_positive_filtered_feature_bc_matrix.h5")

so <- CreateSeuratObject(cbind(countsN1, countsP1, countsN2, countsP2))
so$channel <- c(rep("negative", dim(countsN1)[2]), 
                rep("positive", dim(countsP1)[2]),
                rep("negative", dim(countsN2)[2]),
                rep("positive", dim(countsP2)[2])
)
so$donor <- c(rep("d1", dim(countsN1)[2]), 
                rep("d1", dim(countsP1)[2]),
                rep("d2", dim(countsN2)[2]),
                rep("d2", dim(countsP2)[2])
)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nCount_RNA >= 1000 & nFeature_RNA >= 500  & percent.mt < 5)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters(resolution = 0.4)


#### 
# make basic plots
####

mk_plot <- function(gene){
  pu <- FeaturePlot(so, features = c( gene),  
                    pt.size = 0.1, max.cutoff = "q90") + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot("DCN"),
    mk_plot("VWF"),
    mk_plot("FN1"),
    ncol = 2, scale = 1
  ), file = "../output/ffpe_markers_umap_main.png",width = 1.5*4, height = 1.5*4, dpi = 600)


pumapb <- DimPlot(so, label = FALSE, group.by = c( "channel"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  scale_color_manual(values = jdb_palette("corona")[c(4,5)])+
  theme_void() + ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(pumapb, file = "../output/umap_base.png", width = 3, height = 3, dpi = 400)

cor(t(data.matrix(so@assays$RNA@scale.data[c("VWF", "DCN", "FN1","COL1A1"),])), method = "pearson") %>%
  pheatmap::pheatmap()


####
# find top genes
####

fm <- FindMarkers(so, ident.1 = "positive", ident.2 = "negative", group.by = "channel", 
                  only.pos = FALSE, min.pct = 0.02, logfc.threshold = 0.05)
fm$gene <- rownames(fm)
pFC <- fm %>% arrange(desc(avg_log2FC)) %>% 
  mutate(rank = 1:n()) %>%
  mutate( cg = case_when(
    gene %in% c("DCN", "VWF", "FN1") ~ "panel", 
    grepl("^COL", gene) ~ "COL", 
    TRUE ~ "other")) %>%
  ggplot(aes(x = rank, y = avg_log2FC, color = cg)) +
  scale_x_log10() + 
  geom_point(size = 1) +
  scale_color_manual(values = jdb_palette("corona")[c(1,8,5)]) + 
  coord_cartesian(ylim = c(-4, 6)) + 
  pretty_plot(fontsize = 7) + 
  theme(legend.position = "none") + L_border() + 
  labs(x = "Genes ranked by logFC", y = "logFC Panel+ / Panel -")
pFC
cowplot::ggsave2(pFC, file = "../output/rank_FC_ffpe.pdf", width = 1.4, height = 1.4)


fm[c("DCN", "VWF", "FN1"),]

###
# CDF
###

so$on_target <- so@assays$RNA@counts["FN1", ] + so@assays$RNA@counts["DCN", ] + so@assays$RNA@counts["VWF", ]
p1 <- so@meta.data %>%
  ggplot(aes(x = on_target , color = channel)) + 
  stat_ecdf() + coord_cartesian(xlim = c(0, 30)) +
  scale_color_manual(values = jdb_palette("corona")[c(4,5)])+
  pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none") +
  labs(x = "Panel UMI count")
p1
cowplot::ggsave2(p1, file = "../output/ecdf_ffpe_count.pdf", width = 1.4, height = 1.4)


###########################
# Subset to the one cluster
#########################
so_ss <- subset(so, subset = seurat_clusters %in% c(5) & channel == "positive")
so_ss <- NormalizeData(so_ss) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30) 
so_ss <- so_ss %>% FindClusters(resolution = 0.6)

DimPlot(so_ss, group.by = "seurat_clusters", label = TRUE)

cowplot::ggsave2(
  DimPlot(so_ss, group.by = "seurat_clusters", pt.size = 0.1) + theme_void() + FontSize(main = 0.0001) + 
    theme_void() +
    theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  , file = "../output/subset_channel_umap_ffpe.png",width = 1.5*1, height = 1.5*1, dpi = 600)


fam <- FindAllMarkers(so_ss, only.pos = TRUE)
fam %>% arrange(desc(avg_log2FC)) %>% filter(cluster == 0) %>% head(20)

FeaturePlot(so_ss, c("FN1", "DCN", "VWF", "FCGBP", "APOD", "OGN", "CD74", "C1QC", "MKI67", "VEGFA",
                     "PLVAP","DLL4", "PDGFRB"))

mk_plot_ss <- function(gene){
  pu <- FeaturePlot(so_ss, features = c( gene),  
                    pt.size = 0.1, max.cutoff = "q90", sort = TRUE) + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}
cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot_ss("PDGFRB"),
    mk_plot_ss("PLVAP"),
    mk_plot_ss("MKI67"),
    mk_plot_ss("OGN"),
    ncol = 4, scale = 1
  ), file = "../output/subset_ffpe_markers_umap_main.png",width = 1.5*4, height = 1.5*1, dpi = 600)

