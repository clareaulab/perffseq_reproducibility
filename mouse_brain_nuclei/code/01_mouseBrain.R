library(Seurat)
library(dplyr)
library(BuenColors)

countsP <- Read10X_h5("../data/240107_Mobpplus_filtered_feature_bc_matrix.h5")
countsN <- Read10X_h5("../data/240107_Mobpminus_filtered_feature_bc_matrix.h5")

# Initialize Seurat object
so <- CreateSeuratObject(cbind(countsP, countsN))
so$channel <- c(rep("positive", dim(countsP)[2]), rep("negative", dim(countsN)[2]))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
summary(so$percent.mt)
summary(so$nCount_RNA)
summary(so$nFeature_RNA)

# subset and normalize and cluster
so <- subset(so, subset = nCount_RNA >= 1000 & nFeature_RNA >= 500  & percent.mt < 5)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters(resolution = 0.2)
dim(so)

FindMarkers(so, ident.1 = "7", ident.2 = "3", group.by = "seurat_clusters", logfc.threshold = 0.5, only.pos = TRUE) %>% 
  head(20)

# Find DEGs
fm <- FindMarkers(so, ident.1 = "positive", ident.2 = "negative", group.by = "channel")
fm %>% arrange(desc(avg_log2FC)) %>% tail(20)
fm["Mobp",]


mk_plot <- function(gene){
  pu <- FeaturePlot(so, features = c( gene),  
                    pt.size = 0.1, max.cutoff = "q90") + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}

cowplot::ggsave2(
  mk_plot("Mobp")
  , file = "../output/mobp_umap.png",width = 1.5*4, height = 1.5*4, dpi = 600)

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot("Gabra6"),
    mk_plot("Rbfox3"),
    mk_plot("Gad1"),
    mk_plot("Itih3"),
    ncol = 2, scale = 1
  ), file = "../output/maincluster_supp_marker_umap.png", width = 1.5*4*2, height = 1.5*4*2, dpi = 600)


cowplot::ggsave2(
  DimPlot(so, group.by = "channel", pt.size = 0.1) + theme_void() + FontSize(main = 0.0001) + 
    scale_color_manual(values = c("dodgerblue3","firebrick"))+
    theme_void() +
    theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  , file = "../output/channel_umap.png",width = 1.5*4, height = 1.5*4, dpi = 600)


so$Mobp_count <- so@assays$RNA@counts["Mobp", ]

so@meta.data %>% group_by(channel) %>%
  summarize(mean(Mobp_count > 0))

p1 <- so@meta.data %>%
  ggplot(aes(x = Mobp_count , color = channel)) + 
  stat_ecdf() + coord_cartesian(xlim = c(0, 30)) +
  scale_color_manual(values = c("dodgerblue3","firebrick"))+
  pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none") +
  labs(x = "Mobp UMI count")
cowplot::ggsave2(p1, file = "../output/ecdf_mobp_count.pdf", width = 1.4, height = 1.4)


##################################################

FindMarkers(so,2, 0, logfc.threshold = 2)
FeaturePlot(so, c("Klk6", "S100b", "Ptgds"))

FeaturePlot(so, c("Mobp", "Adcy1", "Cbln3", "Cbln1", "Grm4", "Scn1b"))
FeaturePlot(so, c("Il33", "Il1rl1"))

FeaturePlot(so, "nCount_RNA" , split.by = "channel")

FindMarkers(so,11, 0, subset.ident = , gr)

#############################
# Now subset to Mobp+ cells
DimPlot(so, label = TRUE)
so_ss <- subset(so, subset = seurat_clusters %in% c(0,2) & channel == "positive")
so_ss <- NormalizeData(so_ss) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30) 
so_ss <- so_ss %>% FindClusters(resolution = 0.2)

cowplot::ggsave2(
  DimPlot(so_ss, group.by = "seurat_clusters", pt.size = 0.1) + theme_void() + FontSize(main = 0.0001) + 
    theme_void() +
    theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  , file = "../output/subset_channel_umap.png",width = 1.5*4, height = 1.5*4, dpi = 600)



FeaturePlot(so_ss, features = c("nCount_RNA", "nFeature_RNA"))

if(FALSE){
  fam <- FindAllMarkers(so_ss, only.pos = TRUE)
  fam %>% group_by(cluster) %>% top_n(50, wt = avg_log2FC) %>% data.frame() %>%
    write.table("../output/oligodendrocyte_marker_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

mk_plot_ss <- function(gene){
  pu <- FeaturePlot(so_ss, features = c( gene),  
                    pt.size = 0.1, max.cutoff = "q90", sort.cell = TRUE) + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot_ss("Il33"),
    mk_plot_ss("Klk6"),
    mk_plot_ss("Agt"),
    mk_plot_ss("Sncb"),
    ncol = 2, scale = 1),
  file = "../output/subcluster_main_marker_umap.png",width = 1.5*4, height = 1.5*4, dpi = 600)

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot_ss("Atp1a3"),
    mk_plot_ss("Atp1b1"),
    mk_plot_ss("Atp1a2"),
    mk_plot_ss("Atp1b2"),
    mk_plot_ss("Mobp"),
    mk_plot_ss("Olig1"),
    mk_plot_ss("Olig2"),
    ncol = 4, scale = 1
  ), file = "../output/subcluster_supp_marker_umap.png", width = 1.5*4*2, height = 1.5*4, dpi = 600)

