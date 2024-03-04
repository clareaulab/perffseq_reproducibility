library(dplyr)
library(BuenColors)
library(Seurat)

v0 <- sort(colSums(Read10X_h5("../seq_data_knee/20220930_CD20_sample_raw_probe_bc_matrix.h5")), decreasing = TRUE)
pbmc <- sort(colSums(Read10X_h5("../seq_data_knee/20220930_PBMC_sample_raw_probe_bc_matrix.h5")), decreasing = TRUE)

data.frame(
  rank = c(1:length(v0), 1:length(pbmc)),
  counts = c(v0, pbmc),
  what = c(rep("v0", length(v0)), rep("flex", length(pbmc)))
) %>% filter(counts > 0) -> knee_df

p1 <- ggplot(knee_df %>% filter(counts > 2), aes(x = rank, y = counts, color = what)) + 
  geom_point(size = 0.1) + scale_x_log10() + scale_y_log10() +
  pretty_plot(fontsize = 4) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick", "dodgerblue3"))
cowplot::ggsave2(p1, file = "../plots_for_paper/knee_plot_v0.pdf", width = 0.95, height = 0.8)
