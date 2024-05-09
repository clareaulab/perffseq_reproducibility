library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
library(stringr)
library(rhdf5)

# Import probes with CD3 association
probes <- h5read(molecule_info_file, "probes")$probe
cd3_idx <- which(grepl("CD3E|CD3D", h5read(molecule_info_file, "probes")$probe))
probes[cd3_idx]

process_sample_CD3idx <- function(molecule_info_file, cell_barcodes){
  #Import molecular counts
  probe_idx_per_molecule <- h5read(molecule_info_file, "probe_idx") + 1
  barcode_idx_per_molecule <- h5read(molecule_info_file, "barcode_idx")
  feature_per_molecule <- h5read(molecule_info_file, "feature_idx")
  features <- h5read(molecule_info_file, "features")[["name"]]
  gex_barcode_per_molecule <- h5read(molecule_info_file, "barcodes")[barcode_idx_per_molecule + 1]
  
  # Filter for in cells
  in_cells <- gex_barcode_per_molecule %in% cell_barcodes
  probe_idx_per_molecule_cells <- probe_idx_per_molecule[in_cells]
  gex_barcode_per_molecule_cells <- gex_barcode_per_molecule[in_cells]
  num_cells <- length(cell_barcodes)
  data.frame(
    gex_barcode_per_molecule_cells,
    probe_idx_per_molecule_cells
  ) %>%
    filter(probe_idx_per_molecule_cells %in% cd3_idx) %>%
    group_by(probe_idx_per_molecule_cells, gex_barcode_per_molecule_cells) %>%
    summarize(count = n()) %>%
    ungroup() %>% 
    group_by(probe_idx_per_molecule_cells) %>% 
    summarize(sum(count)/num_cells)
}  

# Import counts
mat_1 <- Seurat::Read10X_h5("../data/cd3bench_cr_Unsorted_filtered_feature_bc_matrix.h5")
molecule_info_file1 <- "../../../../ps-large-data-files/molecule_info/perffseq/pbmc/Unsorted_sample_molecule_info.h5"
cell_barcodes1 <- stringr::str_remove(colnames(mat_1), "-1$")

mat_2 <- Seurat::Read10X_h5("../data/cd3bench_cr_No_sort_no_probe_filtered_feature_bc_matrix.h5")
molecule_info_file2 <- "../../../../ps-large-data-files/molecule_info/perffseq/pbmc/No_sort_no_probe_sample_molecule_info.h5"
cell_barcodes2 <- stringr::str_remove(colnames(mat_2), "-1$")

mat_3 <- Seurat::Read10X_h5("../data/cd3bench_cr_No_sort_yes_probe_filtered_feature_bc_matrix.h5")
molecule_info_file3 <- "../../../../ps-large-data-files/molecule_info/perffseq/pbmc/No_sort_yes_probe_sample_molecule_info.h5"
cell_barcodes3 <- stringr::str_remove(colnames(mat_3), "-1$")

mat_4 <- Seurat::Read10X_h5("../data/cd3bench_cr_CD3_positive_filtered_feature_bc_matrix.h5")
molecule_info_file4 <- "../../../../ps-large-data-files/molecule_info/perffseq/pbmc/CD3_positive_sample_molecule_info.h5"
cell_barcodes4 <- stringr::str_remove(colnames(mat_4), "-1$")

df1 <- process_sample_CD3idx(molecule_info_file1, cell_barcodes1)
df4 <- process_sample_CD3idx(molecule_info_file4, cell_barcodes4)
# all seemingly impacted
