library(data.table)
library(Seurat)
library(dplyr)
library(BuenColors)
library(viridis)
library(rjson)

# Import data
bcl_mat <- Read10X_h5("../data/counts/BCL11A_filtered_feature_bc_matrix.h5")
spi_mat <- Read10X_h5("../data/counts/SPI1_filtered_feature_bc_matrix.h5")
neg_mat <- Read10X_h5("../data/counts/Negative_filtered_feature_bc_matrix.h5")

# Import TF target data
json_data <- fromJSON(paste0(readLines('../data/targets/BCL11A-targets.json')))
bcl11a_targets <- sapply(json_data$associations, function(x){x$gene$symbol})
json_data <- fromJSON(paste0(readLines('../data/targets/SPI1-targets.json')))
spi_targets <- sapply(json_data$associations, function(x){x$gene$symbol})

# Build seurat object
mat <- cbind(bcl_mat,spi_mat, neg_mat)
mdf <- data.frame(
  row.names = colnames(mat),
  sampleID = c(rep("BCL11A", dim(bcl_mat)[2]),
               rep("SPI1", dim(spi_mat)[2]), 
               rep("Neg", dim(neg_mat)[2]))
)

so <- CreateSeuratObject(counts = mat,
                         meta.data = mdf)

so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )

table(so$sampleID)

so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so <- so %>% FindClusters()
DimPlot(so, label = TRUE, group.by = c("seurat_clusters", "sampleID"), shuffle = TRUE)


set.seed(2020)
pumapb <- DimPlot(so, label = FALSE, group.by = c( "sampleID"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  scale_color_manual(values = c("firebrick", "grey","dodgerblue3")) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(pumapb, file = "../output/umap_base.png", width = 3, height = 3.6, dpi = 400)


#FeaturePlot(so, c("nFeature_RNA"))
FeaturePlot(so, c("BCL11A", "FCGR3A", "CLEC4C", "MS4A1", "CD14", "CD3E", "CD34", "CD4", "CD8A", "NCAM1", "PPBP", "LEF1"))


doublets <- c(6, 13, 15, 20, 22)
bcells <- c(11)
cd34 <- c(21)
platlet <- c(18)
asdcs <- c(19)
pdcs <- c(0)
nCD4T <- c(2)
mCD4T <- 9
nCD8T <- c(5)
mCD8T <- c(7,14)
NKcell <- c(12,16)
cd16mono <- c(8,17)
cd14mono <- c(1,4,10)
cdc <- 3

so$CL_anno <- case_when(
  so$seurat_clusters %in% doublets ~ "doublets",
  so$seurat_clusters %in% bcells ~ "Bcell",
  so$seurat_clusters %in% cd34 ~ "HSPC",
  so$seurat_clusters %in% platlet ~ "Platelet",
  so$seurat_clusters %in% asdcs ~ "ASDCs",
  so$seurat_clusters %in% pdcs ~ "pDCs",
  so$seurat_clusters %in% nCD4T ~ "NaiveCD4T",
  so$seurat_clusters %in% nCD8T ~ "NaiveCD8T",
  so$seurat_clusters %in% mCD4T ~ "MemCD4T",
  so$seurat_clusters %in% mCD8T ~ "MemCD8T",
  so$seurat_clusters %in% NKcell ~ "NKcell",
  so$seurat_clusters %in% cd16mono ~ "CD16Mono",
  so$seurat_clusters %in% cd14mono ~ "CD14Mono",
  so$seurat_clusters %in% cdc ~ "cDC",
  TRUE ~ "help"
)

color_vec <- c(
  "purple3", "pink2",
  "dodgerblue2", "dodgerblue4", "navyblue","lightgrey",
  jdb_palette("corona")[c(9:14)], 
  "red2", "darkgrey")

names(color_vec) <- sort(unique(so$CL_anno))
DimPlot(so, group.by = "CL_anno", shuffle = TRUE) +
  scale_color_manual(values = color_vec)

so@meta.data %>% group_by(sampleID, CL_anno) %>% 
  summarize(count = n()) %>%
  group_by(sampleID) %>%
  mutate(prop = count / sum(count)) -> melt_prop_df 

melt_prop_df$sampleID_order <- factor(as.character(melt_prop_df$sampleID), c("Neg","BCL11A", "SPI1"))
pX <- melt_prop_df %>% 
  ggplot(aes(x = sampleID_order, y = prop*100, fill = CL_anno)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = color_vec) +
  pretty_plot(fontsize = 6) + 
  L_border() +theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "% of cells in library")
cowplot::ggsave2(pX, file = "../output/stacked_bars.pdf", width = 1.5, height = 1.1 )

set.seed(2020)
cumapb <- DimPlot(so, label = FALSE, group.by = c( "CL_anno"), shuffle = TRUE, seed = 3,pt.size = 0.01) +
  scale_color_manual(values = color_vec) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")
cowplot::ggsave2(cumapb, file = "../output/umap_base_clusters.png", width = 3, height = 3.6, dpi = 400)

if(FALSE){
  so$BCL11Acount <- so@assays$RNA@counts["BCL11A", ]
  so$SPI1count <- so@assays$RNA@counts["SPI1", ]
  
  so$BCL11Ascale <- so@assays$RNA@scale.data["BCL11A", ]
  so$SPI1scale <- so@assays$RNA@scale.data["SPI1", ]
  
  p1 <- so@meta.data %>%
    ggplot(aes(x = BCL11Acount , color = sampleID)) + 
    stat_ecdf() + coord_cartesian(xlim = c(0, 30)) +
    scale_color_manual(values = c("firebrick", "grey","dodgerblue3"))+
    pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none")
  cowplot::ggsave2(p1, file = "../output/ecdf_BCL11A_count.pdf", width = 1.4, height = 1.4)
  
  p1 <- so@meta.data %>%
    ggplot(aes(x = SPI1count , color = sampleID)) + 
    stat_ecdf() + coord_cartesian(xlim = c(0, 30)) +
    scale_color_manual(values = c("firebrick", "grey","dodgerblue3"))+
    pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none")
  cowplot::ggsave2(p1, file = "../output/ecdf_SPI1_count.pdf", width = 1.4, height = 1.4)
  
  so@meta.data %>% group_by(sampleID) %>%
    summarize(mean(BCL11Acount > 0)*100, mean(SPI1count > 0)*100)
  
  p2 <- so@meta.data %>%
    ggplot(aes(x = BCL11Ascale , color = sampleID)) + 
    stat_ecdf() + 
    scale_color_manual(values = c("firebrick", "grey","dodgerblue3"))+
    pretty_plot(fontsize = 6) + L_border() + theme(legend.position = "none")
  cowplot::ggsave2(p2, file = "../output/ecdf_BCL11A_scale.pdf", width = 1.6, height = 1.6)
  
  p2 <- so@meta.data %>%
    ggplot(aes(x = SPI1scale , color = sampleID)) + 
    stat_ecdf() + 
    scale_color_manual(values = c("firebrick", "grey","dodgerblue3"))+
    pretty_plot(fontsize = 6) + L_border() + theme(legend.position = "none")
  cowplot::ggsave2(p2, file = "../output/ecdf_SPI1_scale.pdf", width = 1.6, height = 1.6)
  
}

p1 <- so@meta.data %>% group_by(sampleID, CL_anno) %>% #seurat_clusters
  summarize(count = n()) %>%
  group_by(CL_anno) %>%
  mutate(prop = count / sum(count)) %>%
  reshape2::dcast(., CL_anno ~ sampleID, value.var = "prop", fill = 0) %>%
  ggplot(aes(x = log2((BCL11A+0.0001)/(Neg + 0.0001)), y = log2((SPI1+0.0001)/(Neg + 0.0001)), color = CL_anno, label = CL_anno)) +
  geom_point() + 
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = color_vec) + pretty_plot(fontsize = 6) +
  theme(legend.position = "none") + coord_cartesian(xlim = c(-6, 10), ylim = c(-6,10)) +
  L_border() + labs(x = "log2 BCL11A / Neg", y = "log2 SPI1 / Neg")
p1
cowplot::ggsave2(p1, file = "../output/population_freq.pdf", width = 1.1, height = 1.1)


ps <- so@meta.data %>% group_by(sampleID, seurat_clusters) %>% #seurat_clusters
  summarize(count = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(prop = count / sum(count)) %>%
  reshape2::dcast(., seurat_clusters ~ sampleID, value.var = "prop", fill = 0) %>%
  ggplot(aes(x = log2((BCL11A+0.0001)/(Neg + 0.0001)), y = log2((SPI1 +0.0001)/(Neg + 0.0001)), label = seurat_clusters)) +
  geom_point() + pretty_plot(fontsize = 6) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  L_border() + labs(x = "log2 BCL11A / Neg", y = "log2 SPI1 / Neg")
cowplot::ggsave2(ps, file = "../output/supplement_cluster_freq.pdf", width = 1.5, height = 1.5)


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
    mk_plot("CD4"),
    mk_plot("NKG7"),
    mk_plot("MS4A1"), 
    mk_plot("CEBPB"), 
    mk_plot("IL3RA"), 
    mk_plot("CLEC4C"), 
    mk_plot("AXL"),
    ncol = 4, scale = 1
  ), file = "../output/markers_umap_supp.png",width = 3*4, height = 1.5*4, dpi = 600)

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot("BCL11A"),
    mk_plot("SPI1"),
    mk_plot("CD3E"),
    ncol = 2, scale = 1
  ), file = "../output/markers_umap_twoTFs.png",width = 3*2, height = 3*2, dpi = 600)

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


########
# B cell
########

so_bcell <- subset(so, seurat_clusters == 11 )
so_bcell <- NormalizeData(so_bcell) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_bcell <- so_bcell %>% FindClusters(resolution = 0.5)
DimPlot(so_bcell, group.by = c("seurat_clusters", "sampleID"), label = TRUE, shuffle = TRUE)

# Now compute target expression
so_bcell <- AddModuleScore(so_bcell, list(bcl11a_targets))
so_bcell$BCL11a_counts <- (so@assays$RNA@counts["BCL11A",]) + 0.5
so_bcell$borno <- ifelse(so$sampleID == "BCL11A", "BCL11A", "other")

pA <- ggplot(so_bcell@meta.data, aes(x = borno, y = BCL11a_counts, fill = borno)) + 
  geom_violin(color = "black") + scale_y_log10() + scale_fill_manual(values = c("firebrick", "darkgrey")) +
  pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none") +
  labs(x = "", y = "BCL11A UMI counts")

pB <- ggplot(so_bcell@meta.data, aes(x = borno, y = Cluster1, fill = borno)) + 
  geom_violin(color = "black")  + scale_fill_manual(values = c("firebrick", "darkgrey")) +
  pretty_plot(fontsize = 5) + L_border() + theme(legend.position = "none") +
  labs(x = "", y = "BCL11A target genes")

cowplot::ggsave2(pA, file = "../output/violin_BCL11A_counts.pdf", width = 1.3, height = 1.3)
cowplot::ggsave2(pB, file = "../output/violin_BCL11A_targets.pdf", width = 1.3, height = 1.3)

wilcox.test(Cluster1 ~ borno, data = so_bcell@meta.data) %>% str()
wilcox.test(BCL11a_counts ~ borno, data = so_bcell@meta.data) %>% str()

##########
# ASDC 
#########

# ASDC analysis
so_asdc <- subset(so, seurat_clusters == 19 & sampleID %in% c("BCL11A", "SPI1"))
so_asdc <- NormalizeData(so_asdc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30)
so_asdc <- so_asdc %>% FindClusters(resolution = 0.5)
DimPlot(so_asdc, group.by = c("seurat_clusters", "sampleID"), label = TRUE, shuffle = TRUE)

pV <- VlnPlot(so_asdc, group.by = "sampleID", 
        features = c( "SPI1", "BCL11A","ITGAX","IL3RA",
                      "GZMB", "TLR9", "IL1B", "CD5"), ncol = 4) &
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) &
  theme_void() & theme(legend.position = "none") &
  FontSize(main = 0.0001, x.text = 4, y.text = 4) & labs(x = "", y = "")
pV
cowplot::ggsave2(pV, file = "../output/violin.pdf", width = 3.1, height = 1.2)

pV <- VlnPlot(so_asdc, group.by = "sampleID", 
              features = c( "IFI30", "MZB1"), ncol = 1) &
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) &
  theme_void() & theme(legend.position = "none") &
  FontSize(main = 0.0001, x.text = 4, y.text = 4) & labs(x = "", y = "")
pV
cowplot::ggsave2(pV, file = "../output/violin_supp.pdf", width = 1.2, height = 1.2)


fm <- FindMarkers(so_asdc, group.by = "sampleID", "BCL11A", "SPI1",
                  logfc.threshold = 0.025, min.pct = 0.05) 

fm %>% arrange(p_val_adj) %>%
  write.table("../output/asdc_DEG.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

gg <- c( "SPI1", "BCL11A","ITGAX","IL3RA",
         "GZMB", "TLR9", "IL1B", "CD5", "IFI30", "MZB1")
fm$gene <- rownames(fm)

fm[gg,]

p1 <- fm %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = gene %in% gg)) + 
  geom_point(size = 0.5) +
  pretty_plot(fontsize = 7) + L_border() + 
  labs(x = "log2FC ()", y = "-log10padj") +
  scale_color_manual(values = c("black", "firebrick")) +
  coord_cartesian(xlim = c(-3, 4), ylim = c(0, 18)) +
  theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../output/volcano_BCL11A_SPI.pdf", width = 2.2, height = 1.2)


mk_plot_asdc <- function(gene){
  pu <- FeaturePlot(so_asdc, features = c( gene),  
                    pt.size = 1, max.cutoff = "q90") + FontSize(main = 0.0001) + 
    theme_void() + theme(legend.position = "none") + ggtitle("") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_gradientn(colors = c("lightgrey", jdb_palette("solar_rojos")[c(2:9)]))
  return(pu)
}

cowplot::ggsave2(
  cowplot::plot_grid(
    mk_plot_asdc("AXL"),
    mk_plot_asdc("SIGLEC6"),
    ncol = 2, scale = 1
  ), file = "../output/markers_umap_sup.png",width = 1.6*4, height = 0.8*4, dpi = 600)

pumapb <- DimPlot(so_asdc, label = FALSE, group.by = c( "sampleID"), shuffle = TRUE, seed = 3,pt.size = 1) +
  scale_color_manual(values = c("firebrick", "dodgerblue3")) + 
  theme_void() + ggtitle("") + theme(legend.position = "none")
pumapb
cowplot::ggsave2(pumapb, file = "../output/umap_asdc.png", width = 3.2, height = 2.2, dpi = 600)




