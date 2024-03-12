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
  sampleID = c(rep("1flex", dim(flex)[2]),
               rep("2nsnp", dim(nsnp)[2]),
               rep("3nsyp", dim(nsyp)[2]),
               rep("4CD3m", dim(cd3m)[2]))
)

so <- CreateSeuratObject(counts = mat,
                         meta.data = mdf)


so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )

so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters(resolution = 0.4)

DimPlot(so, label = TRUE, group.by = c("sampleID", "seurat_clusters"), shuffle = TRUE) 

pu_split <- DimPlot(so, label = FALSE, group.by = c("sampleID"), shuffle = TRUE) +
  facet_wrap(~sampleID, nrow = 1) +
  scale_color_manual(values = jdb_palette("corona")[c(5,2,3,1)] ) +
  theme_void() + ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(pu_split, file = "../plots/4plex_umap.png", height= 2*2, width = 7.2*2)

pu_1 <- DimPlot(so, label = FALSE, group.by = c("sampleID"), shuffle = TRUE, pt.size = 0.1) +
  scale_color_manual(values = jdb_palette("corona")[c(5,2,3,1)] ) +
  theme_void() + ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(pu_1, file = "../plots/single_umap.png", height= 3*4, width = 3*4, dpi = 400)


so$CD3Dc <- (so@assays$RNA@counts["CD3D", ])
so$CD3Ec <- (so@assays$RNA@counts["CD3E", ])

so_t <- subset(so, seurat_clusters %in%c("0", "1", "2", "4", "5"))
library(ggpubr)
pcore <- ggplot(so_t@meta.data, aes(x = log1p(CD3Dc), y = log1p(CD3Ec))) +
  geom_point() + facet_wrap(~sampleID, nrow = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  pretty_plot(fontsize = 7)
cowplot::ggsave2(pcore, file = "../plots/cd3_core.pdf", height= 1.8, width = 7.2)

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
    mk_plot("CD3D"),
    mk_plot("CD8A"),
    mk_plot("CD4"),
    mk_plot("NKG7"),
    mk_plot("MS4A1"), 
    mk_plot("CEBPB"), 
    mk_plot("CD14"), 
    mk_plot("CLEC4C"), 
    ncol = 4, scale = 1
  ), file = "../plots/4plex_markers_umap_supp.png",width = 3*4, height = 1.5*4, dpi = 600)

