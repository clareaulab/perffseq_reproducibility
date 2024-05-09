library(Seurat)
library(dplyr)
library(BuenColors)
library(harmony)
samples <- gsub(".filteredcounts.h5", "", list.files("../data/counts/"))
lapply(samples, function(dir_base){
  
  mat <- Read10X_h5(paste0("../data/counts/", dir_base, ".filteredcounts.h5"))
  colnames(mat) <- paste0(dir_base, colnames(mat))
  mat
}) %>% do.call("cbind",.) -> big_mat

lapply(samples, function(dir_base){
  dt <- readRDS(paste0("../output/azimuth_pbmc_",dir_base,".rds"))
  dt$barcode <- paste0(dir_base, rownames(dt))
  dt
}) %>% rbindlist() %>% data.frame() -> mdf
rownames(mdf) <- mdf$barcode
common <- intersect(mdf$barcode, colnames(big_mat))

# Create seurat object and do seurat things
so <- CreateSeuratObject(
  counts = big_mat, meta.data = mdf
)


so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )
dim(so)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  RunHarmony(., c("donor")) %>%
  RunUMAP(dims = 1:30, reduction = "harmony")

so <- so %>% FindNeighbors(reduction = "harmony")
so <- so %>% FindClusters(resolution = 0.3)

DimPlot(so, group.by = c("seurat_clusters", "donor", "sort"), shuffle = TRUE, 
        split.by = "sort")
DimPlot(so, group.by = c("predicted.celltype.l2"), label = TRUE)
table(so$name)

fm <- FindMarkers(so, ident.1 = "neg_bottom", ident.2 = "pos", group.by = "sort",)

fm2 <- FindMarkers(so, ident.1 = "neg_bottom", ident.2 = "pos", group.by = "sort",subset.ident = "2")
fm0 <- FindMarkers(so, ident.1 = "neg_bottom", ident.2 = "pos", group.by = "sort",subset.ident = "0")
fm7 <- FindMarkers(so, ident.1 = "neg_bottom", ident.2 = "pos", group.by = "sort",subset.ident = "7")


FeaturePlot(so, features = c("FOXP3", "CTLA4", "TCL1A"), split.by = "sort", sort.cell = TRUE)
