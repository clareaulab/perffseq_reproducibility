library(SingleR)
library(celldex)
library(Seurat)
library(dplyr)

# Import data
cd4_mat <- Read10X_h5("../data/CD3Pos_CD4Pos_filtered_feature_bc_matrix.h5"); cd4_mat <- log1p(cd4_mat[,colSums(cd4_mat) >= 1000])
cd8_mat <- Read10X_h5("../data/CD3Pos_CD4Neg_filtered_feature_bc_matrix.h5"); cd8_mat <- log1p(cd8_mat[,colSums(cd8_mat) >= 1000])
cd20_mat <- Read10X_h5("../data/CD20_filtered_feature_bc_matrix.h5"); cd20_mat <- log1p(cd20_mat[,colSums(cd20_mat) >= 1000])
dn_mat <- Read10X_h5("../data/Double_Negative_filtered_feature_bc_matrix.h5"); dn_mat <- log1p(dn_mat[,colSums(dn_mat) >= 1000])

# Import Monaco immune signature
MID <- MonacoImmuneData(ensembl=FALSE)
common <- intersect(rownames(cd20_mat), rownames(MID))

set.seed(2001)
trained <- trainSingleR(MID[common,], labels=MID$label.main, aggr.ref=TRUE)

# Annotate 
pred_cd20 <- classifySingleR(cd20_mat, trained, assay.type=1)
pred_cd4 <- classifySingleR(cd4_mat, trained, assay.type=1)
pred_cd8 <- classifySingleR(cd8_mat, trained, assay.type=1)
pred_dn <- classifySingleR(dn_mat, trained, assay.type=1)

cd_cd20 <- sum((table(pred_cd20$labels)/length(pred_cd20$labels))[c("B cells")])
cd_cd4 <- sum((table(pred_cd4$labels)/length(pred_cd4$labels))[c("CD4+ T cells")] )
cd_cd8 <- sum((table(pred_cd8$labels)/length(pred_cd8$labels))[c("CD8+ T cells", "T cells")] )
cd_dn <- sum((table(pred_dn$labels)/length(pred_dn$labels))[c("NK cells", "Progenitors", "Monocytes", "Dendritic cells")])


# Enrich based on azimuth
lapply(list.files("../output/", full.names = TRUE, pattern = "*projected*"), readRDS) %>%
  rbindlist() %>%
  group_by(predicted.celltype.l2, name) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(name) %>% mutate(prop = count/ sum(count)) -> prop_df

azi4 <- prop_df %>% filter(name == "CD3Pos_CD4Pos") %>%
  filter(grepl("^CD4",  predicted.celltype.l2)) %>% pull(prop) %>% sum()

aziB <- prop_df %>% filter(name == "CD20") %>%
  filter(grepl("^B",  predicted.celltype.l2)) %>% pull(prop) %>% sum()

azi8<- prop_df %>% filter(name == "CD3Pos_CD4Neg") %>%
  filter(grepl("^CD8",  predicted.celltype.l2)) %>% pull(prop) %>% sum() 

aziDN<- prop_df %>% filter(name == "Double_Negative") %>%
  filter(!grepl("^CD8",  predicted.celltype.l2) & !grepl("^B",  predicted.celltype.l2) & !grepl("^CD4",  predicted.celltype.l2)) %>% pull(prop) %>% sum() 


data.frame(
  method = c(rep("azimuth", 4), rep("celldex", 4), rep("manual", 4)),
  what = rep(c("CD20", "CD8T", "CD4T", "DN"), 3),
  perc = c(aziB, azi8, azi4, aziDN,
           cd_cd20, cd_cd8, cd_cd4,cd_dn,
           0.87, 0.878, 0.854, 0.74)*100
) -> plot_df
plot_df
plot_df %>%
  ggplot(aes(x = method, y = perc, fill = what)) + 
  geom_bar(stat = "identity", color = "black", position="dodge") +
  scale_y_continuous(expand = c(0,0), limits = c(0,105)) + pretty_plot(fontsize = 7) + L_border() +
  scale_fill_manual(values = c("orange2", "purple2","firebrick", "lightgrey")) +
  labs(x = "annotation method", y = "% cells matching identity") -> p0
cowplot::ggsave2(p0, file = "../output/stacked_bar.pdf", width = 2.4, height = 1.5)
