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
               rep("3nsyp", dim(nsyp)[2]),
               rep("2nsnp", dim(nsnp)[2]),
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

so_t@meta.data %>% 
  mutate(ratio = log2((CD3Ec+0.1)/(CD3Dc+0.1))) %>%
  ggplot(aes(x = sampleID, y = ratio)) + 
  geom_violin()
 

if(FALSE){
  library(ggpubr)
  pcore <- ggplot(so_t@meta.data, aes(x = log1p(CD3Dc), y = log1p(CD3Ec))) +
    geom_point() + facet_wrap(~sampleID, nrow = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    pretty_plot(fontsize = 7)
  cowplot::ggsave2(pcore, file = "../plots/cd3_core.pdf", height= 1.8, width = 7.2)
}

cpm <- function(vec){
  vec/sum(vec)*1000000
}

# Compute CPM
rs1 <- rowSums(so_t@assays$RNA@counts[,so_t@meta.data$sampleID=="1flex"]) %>% cpm
rs2 <- rowSums(so_t@assays$RNA@counts[,so_t@meta.data$sampleID=="2nsnp"]) %>% cpm
rs3 <- rowSums(so_t@assays$RNA@counts[,so_t@meta.data$sampleID=="3nsyp"]) %>% cpm
rs4 <- rowSums(so_t@assays$RNA@counts[,so_t@meta.data$sampleID=="4CD3m"]) %>% cpm
g2 <- c("CD3D", "CD3E")
data.frame(
  what = rep(c("1flex", "2nsnp", "3nsyp", "4CD3m"), each = 2),
  exp = round(c(rs1[g2], rs2[g2], rs3[g2], rs4[g2]),0),
  gene = rep(g2, 4)
) %>% mutate(l2e = log2(exp)) %>% ggplot(aes(x = what, y = log2(exp), fill = gene)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.8) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12)) + scale_fill_manual(values = c("lightblue", "grey"))+
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") -> pX
  
cowplot::ggsave2(pX, file = "../plots/bars_cd3_log2.pdf", width = 2, height = 2)

FindMarkers(so_t, group.by = "sampleID", ident.1 = "4CD3m", ident.2 = "1flex", logfc.threshold = 0.5)
FindMarkers(so_t, group.by = "sampleID", ident.1 = "3nsyp", ident.2 = "1flex", logfc.threshold = 0.5)
FindMarkers(so_t, group.by = "sampleID", ident.1 = "2nsnp", ident.2 = "1flex", logfc.threshold = 0.5)

so_t@meta.data %>% group_by(sampleID) %>%
  summarize(CD3D = mean(CD3Dc > 0), CD3E = mean(CD3Ec > 0)) %>%
  reshape2::melt(id.vars = "sampleID") %>%
  ggplot(aes(x = sampleID, y = value*100, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.8) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() + scale_fill_manual(values = c("lightblue", "grey"))+
  theme(legend.position = "none") -> pY
cowplot::ggsave2(pY, file = "../plots/bars_cd3_pctPos.pdf", width = 2, height = 2)


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

