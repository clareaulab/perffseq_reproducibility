library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(viridis)

cd20_mat <- Read10X_h5("../data/CD20_filtered_feature_bc_matrix.h5")
cd3_mat <- Read10X_h5("../data/CD3Pos_CD4Neg_filtered_feature_bc_matrix.h5")
dp_mat <- Read10X_h5("../data/CD3Pos_CD4Pos_filtered_feature_bc_matrix.h5")
dn_mat <- Read10X_h5("../data/Double_Negative_filtered_feature_bc_matrix.h5")
mat <- cbind(cd20_mat, cd3_mat, dp_mat,dn_mat)
mdf <- data.frame(
  row.names = colnames(mat),
  sampleID = c(rep("CD20p", dim(cd20_mat)[2]),
               rep("CD3pCD4n", dim(cd3_mat)[2]),
               rep("CD3pCD4p", dim(dp_mat)[2]),
               rep("Negative", dim(dn_mat)[2]))
)

so <- CreateSeuratObject(counts = mat,
                         meta.data = mdf)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )

so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters(resolution = 1)
DimPlot(so, label = TRUE, group.by = c("seurat_clusters", "sampleID"), shuffle = TRUE)
FeaturePlot(so, "nCount_RNA")
so$CD4scale <- so@assays$RNA@scale.data["CD4", ]
so$CD8scale <- so@assays$RNA@scale.data["CD8A", ]
so$CD20scale <- so@assays$RNA@scale.data["MS4A1", ]


so$CD4x <- so@assays$RNA@counts["CD4", ]
so$CD8x <- so@assays$RNA@counts["CD8A", ]
so$CD3x <- so@assays$RNA@counts["CD3E", ] +so@assays$RNA@counts["CD3D", ] 
so$CD20x <- so@assays$RNA@counts["MS4A1", ]


cluster_prop_df <- so@meta.data %>% group_by(seurat_clusters) %>%
  summarize(
            cd8t = mean(sampleID == "CD3pCD4n"), 
            cd4t = mean(sampleID == "CD3pCD4p"), 
            bcell = mean(sampleID == "CD20p"), 
            neg = mean(sampleID == "Negative")
  )
cluster_prop_df$max_vec <- c("CD8T", "CD4T", "bcell", "neg")[max.col(cluster_prop_df[,c(2:5)])]
so$anno <- cluster_prop_df$max_vec[as.numeric(as.character(so$seurat_clusters))+1]

# VCompute counts
so@meta.data %>% group_by(sampleID, anno) %>%
  summarize(count = n()) %>%
  ungroup() %>% group_by(sampleID) %>%
  mutate(prop = count /sum(count)) %>% arrange(desc(prop)) %>% head(5) # for the other script

set.seed(2020)
pumapb <- DimPlot(so, label = FALSE, group.by = c( "sampleID"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  scale_color_manual(values = c("orange2", "firebrick","purple2", "lightgrey")) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(pumapb, file = "../output/umap_base.pdf", width = 3, height = 3.6)
cowplot::ggsave2(pumapb, file = "../output/umap_base.png", width = 3, height = 3.6, dpi = 400)


FeaturePlot(so, features = c("NCAM1", "CD3E", "GNLY"))



fm <- FindMarkers(so, group.by =  "sampleID", ident.1 =  "CD3pCD4p", ident.2 = "CD3pCD4n", 
                  logfc.threshold = 0.1, min.pct = 0.05)

fm %>% arrange(avg_log2FC)
p1 <- fm %>%
  mutate(p_val_adj= (p_val_adj + 1e-305) ) %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
  geom_point() +
  coord_cartesian(ylim = c(0, 350), xlim = c(-4, 2.2)) + 
  pretty_plot(fontsize = 7) + L_border() + 
  labs(x = "log2FC ()", y = "-log10padj")
#cowplot::ggsave2(p1, file = "../output/volcano_cd4_cd8.pdf", width = 2, height = 1.5)


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
    mk_plot("CD8A"),
    mk_plot("CD4"),
    mk_plot("CD3E"),
    mk_plot("CD3D"),
    mk_plot("NKG7"),
    mk_plot("MS4A1"), 
    mk_plot("CEBPB"), 
    ncol = 4, scale = 1
  ), file = "../output/markers_umap.png",width = 3*4, height = 1.5*4, dpi = 600)

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot("GNLY"),
    mk_plot("CD19"),
    mk_plot("FCGR3A"),
    mk_plot("CD14"),
    mk_plot("CD34"),
    mk_plot("PPBP"), 
    ncol = 3, scale = 1
  ), file = "../output/markers_umapSUP.png",width = 2.25*4, height = 1.5*4, dpi = 600)


FeaturePlot(so, features = c("CFP", "TNFAIP2", "PTPRE", "CLEC7A", "CEBPB"))

so_CD4 <- subset(so, sampleID == "CD3pCD4p")
so_CD4 <- NormalizeData(so_CD4) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_CD4 <- so_CD4 %>% FindClusters(resolution = 0.5)
DimPlot(so_CD4, label = TRUE)
FeaturePlot(so_CD4, c("CD4", "LEF1","FOXP3",  "IFIT1",  "GZMH"))

bonus <- c("CD8A", "CD4", "LEF1", "MKI67")

mk_plot_cd4 <- function(gene){
  pu <- FeaturePlot(so_CD4, features = c( gene),  
                    pt.size = 0.1, max.cutoff = "q90") + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}
cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot_cd4("CD4"),
    mk_plot_cd4("FOXP3"),
    mk_plot_cd4("IFIT1"),
    mk_plot_cd4("GZMH"),
    ncol = 2, scale = 1
  ), file = "../output/markers_umap_CD4.png",width = 1.6*4, height = 1.5*4, dpi = 600)

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot_cd4("LEF1"),
    mk_plot_cd4("TNF"),
    mk_plot_cd4("IFIT2"),
    mk_plot_cd4("NKG7"),
    mk_plot_cd4("MKI67"),
    mk_plot_cd4("CD8A"),
    ncol = 3, scale = 1
  ), file = "../output/markers_umap_CD4SUP.png",width = 2.25*4, height = 1.5*4, dpi = 600)


