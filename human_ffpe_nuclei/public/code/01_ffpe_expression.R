library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(viridis)

# Import raw counts from 10x website
raw_ffpeA <- Read10X_h5("../../../../ps-large-data-files/other_data/4plex_human_lungcancer_glioblastoma_scFFPE_multiplex_Multiplex_count_raw_feature_bc_matrix.h5")

# BCs 3 and 4 are GBM
bc001 <- "ACTTTAGG"
bc002 <- "AACGGGAA"
bc003 <- "AGTAGGCT" 
bc004 <- "ATGTTGAC"

raw_ffpe_BC4 <- raw_ffpeA[,substr(colnames(raw_ffpeA), 17, 24) %in% c(bc004, bc003)]
boo <- colSums(raw_ffpe_BC4) >= 1500 & colSums(raw_ffpe_BC4 > 0) >= 750
raw_ffpe_cells <- raw_ffpe_BC4[,boo]
raw_ffpe_notcells <- raw_ffpe_BC4[,!boo]

# Do seurat things
so <- CreateSeuratObject(raw_ffpe_cells) %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30) 
so <- so %>% FindClusters()

DimPlot(so, label = TRUE) 
pumapb <- DimPlot(so, label = TRUE, group.by = c( "seurat_clusters"), 
                  shuffle = TRUE, seed = 3, pt.size = 0.01) +
  theme_void() + ggtitle("") + theme(legend.position = "right")
cowplot::ggsave2(pumapb, file = "../plots/pub_ffpe_umap_base.png", 
                 width = 4, height = 4.2, dpi = 600)

FeaturePlot(so, c("DCN","VWF","FN1"))
FeaturePlot(so, c("DCN","VWF","FN1"))

table(so$seurat_clusters)/length(so$seurat_clusters)


rowSums(raw_ffpe_cells)[c("DCN","FN1","VWF","CLDN5", "CD34")]

fm <- FindMarkers(so, 10, only.pos = TRUE)
fm %>% 
  mutate(stat = -log10(p_val_adj + 1E-300)*avg_log2FC) %>% arrange(desc(stat)) %>% head(20)

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
    mk_plot("COL4A2"),
    ncol = 2, scale = 1
  ), file = "../plots/selected_markers_umap.png", width = 3, height = 3, dpi = 600)



