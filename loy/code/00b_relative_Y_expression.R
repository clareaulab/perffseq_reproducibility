library(Seurat)
library(dplyr)
library(BuenColors)

# import and setup
txg <- Read10X_h5("../data/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex_count_raw_feature_bc_matrix.h5")[["Gene Expression"]]
txg <- txg[,colSums(txg>0) > 1000]
msy_genes <- c("ZFY", "TBL1Y", "USP9Y", "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "KDM5D", "EIF1AY") # no RPS4Y1 in flex
other_ff <- c("BCL11A", "CD3E", "MS4A1", "CD4", "SPI1", "CD8A")

data.frame(gene = rownames(txg), totalC = rowSums(txg), n = dim(txg)[2],
           n_pos = rowSums(txg > 0)) %>%
  mutate(pos_mean = (totalC)/(n_pos+0.001)) %>%
  filter(!grepl("DEPRECATED", gene)) %>%
  arrange(desc(pos_mean)) %>% mutate(rank = 1:n()) %>%
  mutate(color = case_when(
    gene %in% other_ff ~ "others", 
    gene %in% msy_genes ~ "msy", 
    TRUE ~ "blah"
  )) %>%
  arrange((color)) -> val_df

val_df %>% 
  filter(pos_mean > 0) %>%
  ggplot(aes(x = rank, y = pos_mean, color = color)) + 
  geom_point(size = 0.9) +
  scale_color_manual(values = c("lightgrey", "firebrick", "purple3"))+
  pretty_plot(fontsize = 8) + L_border() + scale_y_log10() +
  theme(legend.position = "none") + labs(x = "Rank ordered gene expression",
                                         y = "Mean UMIs / + cell") -> pD

cowplot::ggsave2(pD, file = "../plots/public_umi_range.pdf",
                 width = 1.4, height = 1.4)


val_df %>%
  filter(gene %in% c(msy_genes, other_ff)) %>% 
  group_by(color) %>% summarize(mean(pos_mean))


val_df %>%
  filter(gene %in% c(msy_genes)) %>% 
  summarize(tc = sum(totalC)/ max(n))
