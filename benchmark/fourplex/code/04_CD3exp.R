library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(viridis)
library(rjson)

# Import data
flex <- Read10X_h5("../data/cd3bench_cr_Unsorted_filtered_feature_bc_matrix.h5")
nsnp <- Read10X_h5("../data/cd3bench_cr_No_sort_no_probe_filtered_feature_bc_matrix.h5")
nsyp <- Read10X_h5("../data/cd3bench_cr_No_sort_yes_probe_filtered_feature_bc_matrix.h5")
cd3m <- Read10X_h5("../data/cd3bench_cr_CD3_positive_filtered_feature_bc_matrix.h5")

mat <- cbind(flex, nsyp, nsnp, cd3m)
mdf <- data.frame(
  row.names = colnames(mat),
  sampleID = c(rep("flex", dim(flex)[2]),
               rep("nsnp", dim(nsnp)[2]),
               rep("nsyp", dim(nsyp)[2]),
               rep("CD3m", dim(cd3m)[2]))
)

so <- CreateSeuratObject(counts = mat,
                         meta.data = mdf)


so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )

so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters(resolution = 0.4)
DimPlot(so, label = TRUE, group.by = c("seurat_clusters", "sampleID"), shuffle = TRUE)

mean(so@meta.data$seurat_clusters %in% c(8, 9))
FeaturePlot(so, c("CD3E", "CD3D", "CD4", "CD8A", "MS4A1", "CD14", "NCAM1"))

so$CD3Dc <- log1p(so@assays$RNA@counts["CD3D", ])
so$CD3Ec <- log1p(so@assays$RNA@counts["CD3E", ])

so_t <- subset(so, seurat_clusters %in%c("0", "1", "2", "4", "5"))
library(ggpubr)
ggplot(so_t@meta.data, aes(x = CD3Dc, y = CD3Ec)) +
  geom_point() + facet_wrap(~sampleID) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson") + pretty_plot()


so@meta.data %>%
  ggplot(aes(x = CD3Dc+CD3Ec)) + 
  stat_ecdf() + coord_cartesian(xlim = c(0, 30)) +
  scale_color_manual(values = c("firebrick", "grey","dodgerblue3"))+
  pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../output/ecdf_BCL11A_count.pdf", width = 1.4, height = 1.4)

