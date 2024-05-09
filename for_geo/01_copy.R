library(data.table)
library(dplyr)

# Make a mapping from names in github to what will be sanitized on GEO
fig2 <- list.files("../benchmark/fourplex/data", full.names = TRUE)[-5]
fig3 <- list.files("../logic_gating/data", full.names = TRUE)
fig4 <- list.files("../tfs_pbmcs/data/counts", full.names = TRUE)
fig5m <- list.files("../mouse_brain_nuclei/data", full.names = TRUE)
fig5h <- list.files("../human_ffpe_nuclei/perffseqanalyses/data", full.names = TRUE)

data.frame(
  cp = "cp",
  old = c(fig2, fig3, fig4, fig5m, fig5h),
  geo = paste0(c("Benchmark_PERRF", "Benchmark_NPNS", "Benchmark_YPNS", "Benchmark_Flex",
                 "Multiplex_MS4A1p","Multiplex_CD3pCD4n","Multiplex_CD3pCD4p","Multiplex_DN",
                 "Validation_IL3RAn", "Validation_IL3RAp",
                 "TF_BCL11A","TF_DN","TF_SPI1",
                 "MobpN", "MobpP",
                 "PanelN_D1", "PanelP_D1", "PanelN_D2", "PanelP_D2"), "_filtered_feature_bc_matrix.h5")
) %>%
  write.table(sep = " ", row.names = FALSE, quote = FALSE)

# collect QC
library(Seurat)
allh5 <- list.files(".", pattern = "h5")

lapply(allh5, function(one){
  mat <- Read10X_h5(one)
  data.frame(one, celln = dim(mat)[2], genes = median(colSums(mat>0)), umis = median(colSums(mat)))
})%>% rbindlist() %>% data.frame()
