library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
library(stringr)
library(rhdf5)

# https://github.com/slowkow/saturation/blob/main/saturation.R
# This doesn't do any of the filtering
downsample_analysis <- function(pcr_dups_per_molecule_cells, feature_per_molecule_cells, gex_barcode_per_molecule_cells,
                       prob = 0.1, seed_value = 2020, what = "experiment") {
  set.seed(seed_value)
  my_read_counts <- rbinom(n = length(pcr_dups_per_molecule_cells), size = pcr_dups_per_molecule_cells, prob = prob)
  keep <- my_read_counts > 0
  sdt <- data.table(
    pcr_counts = my_read_counts[keep],
    barcodes = gex_barcode_per_molecule_cells[keep],
    gene = feature_per_molecule_cells[keep]
  )
  per_barcode_nUmis <- sdt[, .(N = sum(pcr_counts >= 1), readsPerCell = sum(pcr_counts)), by = list(barcodes)] 
  per_barcode_nGenes <- sdt[, .(N_both = sum(pcr_counts >= 1)), by = list(barcodes,gene)][, .(N = sum(N_both >= 1)), by = list(barcodes)] 
  
  data.frame(
    prob = prob,
    seed = seed_value,
    what = what,
    mean_reads_per_cell = mean(per_barcode_nUmis$readsPerCell),
    median_reads_per_cell = median(per_barcode_nUmis$readsPerCell),
    mean_umis_per_cell = mean(per_barcode_nUmis$N),
    median_umis_per_cell = median(per_barcode_nUmis$N),
    mean_genes_per_cell = mean(per_barcode_nGenes$N),
    median_genes_per_cell = median(per_barcode_nGenes$N)
  )
}

# Helper function to process stuff
process_sample <- function(molecule_info_file, cell_barcodes, what = "what", probs = c(0.01, 0.1, 0.25, 0.5, 1 , 3, 5, 7.5, 10)/10){
  #Import molecular counts
  pcr_dups_per_molecule <- h5read(molecule_info_file, "count")
  barcode_idx_per_molecule <- h5read(molecule_info_file, "barcode_idx")
  feature_per_molecule <- h5read(molecule_info_file, "feature_idx")
  features <- h5read(molecule_info_file, "features")[["name"]]
  gex_barcode_per_molecule <- h5read(molecule_info_file, "barcodes")[barcode_idx_per_molecule + 1]
  
  # Filter vectors per cell
  in_cells <- gex_barcode_per_molecule %in% cell_barcodes
  feature_per_molecule_cells <- feature_per_molecule[in_cells]
  gex_barcode_per_molecule_cells <- gex_barcode_per_molecule[in_cells]
  pcr_dups_per_molecule_cells <- pcr_dups_per_molecule[in_cells]

  # loop over probabilities
  ds_df <- lapply(probs, function(prob){
    print(prob)
    downsample_analysis(pcr_dups_per_molecule_cells, feature_per_molecule_cells, gex_barcode_per_molecule_cells,
                        prob = prob, seed_value = 2020, what = what)
  }) %>% rbindlist() %>% data.frame()
  ds_df
}

# Import data 
mat_1 <- Seurat::Read10X_h5("../data/240107_Mobpminus_filtered_feature_bc_matrix.h5")
molecule_info_file1 <- "../../../ps-large-data-files/molecule_info/perfseq/mousenuclei/240107_Mobpminus_cr_sample_molecule_info.h5"
cell_barcodes1 <- stringr::str_remove(colnames(mat_1), "-1$")

mat_2 <- Seurat::Read10X_h5("../data/240107_Mobpplus_filtered_feature_bc_matrix.h5")
molecule_info_file2 <- "../../../ps-large-data-files/molecule_info/perfseq/mousenuclei/240107_Mobpplus_cr_sample_molecule_info.h5"
cell_barcodes2 <- stringr::str_remove(colnames(mat_2), "-1$")

mat_3 <- Seurat::Read10X_h5("../../../ps-large-data-files/molecule_info/public_flex/mousenuclei/counts_matrix/Eye_Nuclei_BC1_sample_filtered_feature_bc_matrix.h5")
molecule_info_file3 <- "../../../ps-large-data-files/molecule_info/public_flex/mousenuclei/Eye_Nuclei_BC1_sample_molecule_info.h5"
cell_barcodes3 <- stringr::str_remove(colnames(mat_3), "-1$")


# import data
full_ds_df <- rbind(
  process_sample(molecule_info_file1, cell_barcodes1, what = "MobpM"),
  process_sample(molecule_info_file2, cell_barcodes2, what = "MobpP"),
  process_sample(molecule_info_file3, cell_barcodes3, what = "Public")
) 

pA <- ggplot(full_ds_df, aes(x = mean_reads_per_cell, y = median_umis_per_cell, color = what)) +
  geom_point(size = 0.3) + geom_line() +
  pretty_plot(fontsize = 5) + labs(x = 'Mean reads / cell', y = "Median umis per cell") + L_border() + 
  scale_color_manual(values = jdb_palette("corona")) + theme(legend.position = "none")
cowplot::ggsave2(pA, file = "../output/mouse_nuclei_benchmarking.pdf", width = 1.2 ,height = 0.9)

# compute downsampling probabilities to lowest comparator library
full_ds_df %>% group_by(what) %>% top_n(1, prob) %>%
  ungroup() %>%
  mutate(prop = min(mean_reads_per_cell)/mean_reads_per_cell) %>% data.frame()

compare_df <- rbind(
  process_sample(molecule_info_file1, cell_barcodes1, what = "MobpM", 0.852),
  process_sample(molecule_info_file2, cell_barcodes2, what = "MobpP", 0.56),
  process_sample(molecule_info_file3, cell_barcodes3, what = "Public", 1)
) %>% data.frame()
compare_df
