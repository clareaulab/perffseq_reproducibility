library(dplyr)
library(BuenColors)
library(Seurat)

m1 <- Read10X_h5("../seq_data_knee/20220930_CD20_sample_raw_probe_bc_matrix.h5")
m2 <- Read10X_h5("../seq_data_knee/2023_03_27_CD3_cr_filtered_feature_bc_matrix.h5")

data.frame(
  umis = c(colSums(m1), colSums(m2)),
  genes = c(colSums(m1>0), colSums(m2>0)),
  what = c(rep("av0", dim(m1)[2]), rep("perfseq",dim(m2)[2]))
) %>% group_by(what) %>% arrange(desc(genes)) %>% mutate(rank = 1:n()) -> umi_df

p2 <- ggplot(umi_df %>% filter(rank < 3001), aes(x = what, y = genes, color = what)) + 
  geom_violin() + scale_y_log10() +
  pretty_plot(fontsize = 6) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("darkgrey", "darkgrey"))

p1 <- ggplot(umi_df %>% filter(rank < 3001), aes(x = what, y = umis, color = what)) + 
  geom_violin() + scale_y_log10() +
  pretty_plot(fontsize = 6) + L_border() + theme(legend.position = "none") +
  scale_color_manual(values = c("darkgrey", "darkgrey"))

cowplot::ggsave2(p2, file = "../plots_for_paper/gene_violin.pdf", width = 1.4, height = 1.0)
cowplot::ggsave2(p1, file = "../plots_for_paper/UMI_violin.pdf", width = 1.4, height = 1.0)

umi_df %>% filter(rank < 1001) %>% group_by(what) %>% summarize(median(umis), median(genes))
