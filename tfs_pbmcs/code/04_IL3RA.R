library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(viridis)
library(rjson)

# Import data
IL3RA_mat <- Read10X_h5("../data/counts/230601_IL3RAplus_sample_filtered_feature_bc_matrix.h5")

# Build seurat object
mat <- cbind(IL3RA_mat)
mdf <- data.frame(
  row.names = colnames(mat),
  sampleID = c(rep("Il3RA", dim(IL3RA_mat)[2]))
)

so <- CreateSeuratObject(counts = mat,
                         meta.data = mdf)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )


so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters()
DimPlot(so, label = TRUE, group.by = c("seurat_clusters", "sampleID"), shuffle = TRUE)

FeaturePlot(so, c("BCL11A", "SPI1", "AXL", "SIGLEC6", "CLEC4C","FCGR3A","CSF3R"))
FindMarkers(so,"3", "0")

# ASDC analysis
so_asdc <- subset(so, seurat_clusters == 11 )
so_asdc <- NormalizeData(so_asdc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_asdc <- so_asdc %>% FindClusters(resolution = 0.5)
dim(so_asdc)
DimPlot(so_asdc, group.by = c("seurat_clusters"), label = TRUE, shuffle = TRUE)

genes_plot <- c("SPI1", "BCL11A", "IL3RA", "MZB1", "IFI30", "ITGAX")
# Using log1p counts instead of scale data since not everything is a variable feature
cor(t(data.matrix(log1p(so_asdc@assays$RNA@counts[genes_plot,]))), method = "spearman") %>%
  pheatmap::pheatmap()

so$is_asdc <- so$seurat_clusters == 11
set.seed(2020)
pumapb <- DimPlot(so, label = FALSE, group.by = c( "is_asdc"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  scale_color_manual(values = c( "lightgrey","forestgreen")) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")

cowplot::ggsave2(pumapb, file = "../output/umap_base_il3ra_subset.png", width = 4, height = 4, dpi = 400)
