library(data.table)
library(dplyr)
library(Seurat)
library(BuenColors)
library(pheatmap)
library(viridisLite)

d <- fread("../data/GSE94820_raw.expMatrix_deeper.characterization.set.submission.txt.gz")
data_mat <- data.matrix(d[,-1]); rownames(data_mat) <- d[["V1"]]

ss_df <- t(data_mat[c("SPI1", "BCL11A", "ITGAX","IFI30", "MZB1", "IL3RA"),]) %>% data.frame() 
cor_mat <- cor(log1p(ss_df))
pheatmap(cor_mat)
cell_subset <- ss_df[grepl("AXLSIGLEC6", colnames(data_mat)), ]
cor_mat <- cor(log1p(cell_subset), method = "spearman")
pheatmap(cor_mat)

cor(cell_subset)
ss_df %>% arrange(desc(SPI1))
ggplot(ss_df, aes(x = log1p(SPI1), y = log1p(BCL11A))) + 
  geom_point() 
