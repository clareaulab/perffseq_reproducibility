library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)

# import reference
ref_path <- "/Users/lareauc/Dropbox/main_papers/pearson/pearson_large_data_files/input/pbmc/pbmc_multimodal.h5seurat"

options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat(ref_path)

# Import
import_project_scRNAseq <- function(dir_base){
  pheno <- strsplit(dir_base, "_", 2)[[1]][1]
  data.path <- paste0("../data/", dir_base, "_filtered_feature_bc_matrix.h5")
  raw <- Read10X_h5(filename =  data.path)
  colnames(raw) <- paste0(substr(colnames(raw), 1, 16), "-1")
  
  # Remove crazy high and low expressors
  n_feature_rna <- colSums(raw > 0)
  n_total_rna <- colSums(raw)
  pct_mito <- colSums(raw[grepl("^MT", rownames(raw)), ])/n_total_rna * 100
  qc_cells <- colnames(raw)[pct_mito < 10 & n_total_rna > 1000 & n_feature_rna > 500]
  length(qc_cells)
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), qc_cells]
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- SCTransform(raw)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = raw,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  projected <- MapQuery(
    anchorset = anchors,
    query = raw,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  df <- data.frame(
    projected@meta.data,
    projected@reductions$ref.umap@cell.embeddings)
  df$name <- dir_base
  df$pheno <- pheno
  saveRDS(df, file = paste0("../output//_projected_PBMC_", dir_base, ".rds"))
  dir_base
}

samples <- gsub("_filtered_feature_bc_matrix.h5", "", list.files("../data/"))
lapply(samples, import_project_scRNAseq)


