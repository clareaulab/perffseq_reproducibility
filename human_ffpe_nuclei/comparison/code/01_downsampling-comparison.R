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
mat_1 <- Seurat::Read10X_h5("../../perfseqanalyses/data/x1_positive_filtered_feature_bc_matrix.h5")
molecule_info_file1 <- "../../../../ps-large-data-files/molecule_info/perfseq/gbm/x1_positive_sample_molecule_info.h5"
cell_barcodes1 <- stringr::str_remove(colnames(mat_1), "-1$")

mat_2 <- Seurat::Read10X_h5("../../perfseqanalyses/data/x1_negative_filtered_feature_bc_matrix.h5")
molecule_info_file2 <- "../../../../ps-large-data-files/molecule_info/perfseq/gbm/x1_negative_sample_molecule_info.h5"
cell_barcodes2 <- stringr::str_remove(colnames(mat_2), "-1$")

# cell barcoding funky based on what was available on 10x website
# Import raw counts from 10x website
raw_ffpeA <- Read10X_h5("../../../../ps-large-data-files/other_data/4plex_human_lungcancer_glioblastoma_scFFPE_multiplex_Multiplex_count_raw_feature_bc_matrix.h5")

# BCs 3 and 4 are GBM
bc001 <- "ACTTTAGG"
bc002 <- "AACGGGAA"
bc003 <- "AGTAGGCT" 
bc004 <- "ATGTTGAC"

raw_ffpe_BC3 <- raw_ffpeA[,substr(colnames(raw_ffpeA), 17, 24) %in% c(bc003)]
raw_ffpe_BC4 <- raw_ffpeA[,substr(colnames(raw_ffpeA), 17, 24) %in% c(bc004)]

boo3 <- colSums(raw_ffpe_BC3) >= 1500 & colSums(raw_ffpe_BC3 > 0) >= 750
boo4 <- colSums(raw_ffpe_BC4) >= 1500 & colSums(raw_ffpe_BC4 > 0) >= 750

cell_barcodes3 <- colnames(raw_ffpe_BC3[,boo3]) %>% stringr::str_remove(., "-1$")
cell_barcodes4 <- colnames(raw_ffpe_BC4[,boo4]) %>% stringr::str_remove(., "-1$")

# path to molecular files
mat_3 <- Seurat::Read10X_h5("../../public/data/Glioblastoma_Manual_BC3_sample_filtered_feature_bc_matrix.h5")
molecule_info_file3 <- "../../../../ps-large-data-files/molecule_info/public_flex/gbm/Glioblastoma_Manual_BC3_sample_molecule_info.h5"
cell_barcodes3 <- stringr::str_remove(colnames(mat_3), "-1$")

mat_4 <- Seurat::Read10X_h5("../../public/data/Glioblastoma_Octo_BC4_sample_filtered_feature_bc_matrix.h5")
molecule_info_file4 <- "../../../../ps-large-data-files/molecule_info/public_flex/gbm/Glioblastoma_Octo_BC4_sample_molecule_info.h5"
cell_barcodes4 <- stringr::str_remove(colnames(mat_4), "-1$")


# import data
full_ds_df <- rbind(
  process_sample(molecule_info_file1, cell_barcodes1, what = "Positive"),
  process_sample(molecule_info_file2, cell_barcodes2, what = "Negative"),
  process_sample(molecule_info_file3, cell_barcodes3, what = "Public1"),
  process_sample(molecule_info_file4, cell_barcodes4, what = "Public2")
) 

pA <- ggplot(full_ds_df, aes(x = mean_reads_per_cell, y = median_umis_per_cell, color = what)) +
  geom_point(size = 0.3) + geom_line() +
  pretty_plot(fontsize = 5) + labs(x = 'Mean reads / cell', y = "Median umis per cell") + L_border() + 
  scale_color_manual(values = jdb_palette("corona")[c(4:7)]) + theme(legend.position = "none")
cowplot::ggsave2(pA, file = "../plots/gbm_ffpe_benchmarking.pdf", width = 1.2 ,height = 0.9)

ggplot(full_ds_df, aes(x = mean_reads_per_cell, y = median_umis_per_cell, color = what)) +
  geom_point(size = 0.3) + geom_line() +
  pretty_plot(fontsize = 5) + labs(x = 'Mean reads / cell', y = "Median umis per cell") + L_border() + 
  scale_color_manual(values = jdb_palette("corona")[c(4:7)])

# compute downsampling probabilities to lowest comparator library
full_ds_df %>% group_by(what) %>% top_n(1, prob) %>%
  ungroup() %>%
  mutate(prop = min(mean_reads_per_cell)/mean_reads_per_cell) %>% data.frame()

compare_df <- rbind(
  process_sample(molecule_info_file1, cell_barcodes1, what = "Positive", 0.5014276),
  process_sample(molecule_info_file2, cell_barcodes2, what = "Negative", 0.3296781),
  process_sample(molecule_info_file3, cell_barcodes3, what = "Public1", 1.0000000),
  process_sample(molecule_info_file4, cell_barcodes4, what = "Public2", 0.9349910)
) %>% data.frame()

compare_df


