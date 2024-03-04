library(Seurat)
library(dplyr)

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

# Find DEGs
fm <- FindMarkers(so, ident.1 = "positive", ident.2 = "negative", group.by = "channel")
fm %>% arrange(desc(avg_log2FC)) %>% head(20)
fm["Mobp",]

library(ggplot2)
DimPlot(so, group.by = c("seurat_clusters", "channel"), label = TRUE) &
  geom_vline(xintercept = 0)
DimPlot(so, group.by = c("channel"))
DimPlot(so, group.by = c("seurat_clusters"), label = TRUE)

FeaturePlot(so, c("Mobp"), split.by = "channel") &
  theme_void()

FindMarkers(so,2, 0, logfc.threshold = 2)
FeaturePlot(so, c("Klk6", "S100b", "Ptgds"))

FeaturePlot(so, c("Mobp", "Hopx", "Il33", "Cbln1", "C1qa", "Ppfibp1"), split.by = "channel")
FeaturePlot(so, c("Il33", "Il1rl1"))

FeaturePlot(so, "nCount_RNA" , split.by = "channel")

FindMarkers(so,11, 0, subset.ident = , gr)

# Now subset to Mobp+ cells

so_ss <- subset(so, subset = seurat_clusters %in% c(0,2) & channel == "positive")
so_ss <- NormalizeData(so_ss) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  FindNeighbors() %>% RunUMAP(dims = 1:30) 
so_ss <- so_ss %>% FindClusters(resolution = 0.2)

DimPlot(so_ss, group.by = "seurat_clusters", label = TRUE)

if(FALSE){
  fam <- FindAllMarkers(so_ss, only.pos = TRUE)
  fam %>% group_by(cluster) %>% top_n(50, wt = avg_log2FC) %>% data.frame() %>%
    write.table("../output/oligodendrocyte_marker_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

FeaturePlot(so_ss, c("Mobp", "Il33", "Ptgds",   "Klk6", 
                     "Sncb","Fgfr3", "Mt2", "Agt",
                     "Atp1b1","Atp1a3", "Atp1a2", "Atp1b2"), sort.cell = FALSE)
