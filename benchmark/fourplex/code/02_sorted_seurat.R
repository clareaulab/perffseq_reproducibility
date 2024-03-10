library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(viridis)
library(rjson)

# Import data
sorted_mat <- Read10X_h5("../data/cd3bench_cr_CD3_positive_filtered_feature_bc_matrix.h5")
dim(sorted_mat)

# Build seurat object
so <- CreateSeuratObject(counts = sorted_mat)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters(resolution = 0.4)
DimPlot(so, label = TRUE, group.by = c("seurat_clusters"), shuffle = TRUE)

mean(so@meta.data$seurat_clusters %in% c(8, 9))
FeaturePlot(so, c("CD3E", "CD3D", "CD4", "CD8A", "MS4A1", "CEBPE", "NCAM1"))

so$CD3Dc <- so@assays$RNA@counts["CD3D", ]
so$CD3Ec <- so@assays$RNA@counts["CD3E", ]

so@meta.data %>%
  ggplot(aes(x = CD3Dc+CD3Ec)) + 
  stat_ecdf() + coord_cartesian(xlim = c(0, 30)) +
  scale_color_manual(values = c("firebrick", "grey","dodgerblue3"))+
  pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../output/ecdf_BCL11A_count.pdf", width = 1.4, height = 1.4)

