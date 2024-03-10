library(SingleR)
library(celldex)
library(Seurat)
library(dplyr)

# Import data
cd3_mat <- Read10X_h5("../data/cd3bench_cr_CD3_positive_filtered_feature_bc_matrix.h5"); cd3_mat <- log1p(cd3_mat[,colSums(cd3_mat) >= 1000])

# Import Monaco immune signature
MID <- MonacoImmuneData(ensembl=FALSE)
common <- intersect(rownames(cd3_mat), rownames(MID))

set.seed(2001)
trained <- trainSingleR(MID[common,], labels=MID$label.main, aggr.ref=TRUE)

# Annotate 
pred_cd3 <- classifySingleR(cd3_mat, trained, assay.type=1)
prop_vec <- (table(pred_cd3$labels)/length(pred_cd3$labels))
celldex_prop <- sum(prop_vec[grepl("T", names(prop_vec))])


# Enrich based on azimuth
lapply(list.files("../output/", full.names = TRUE, pattern = "*projected*"), readRDS) %>%
  rbindlist() %>%
  group_by(predicted.celltype.l1, name) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(name) %>% mutate(prop = count/ sum(count)) -> prop_df

azi <- prop_df %>% filter(name == "cd3bench_cr_CD3_positive") %>%
  filter(grepl("T",  predicted.celltype.l1)) %>% pull(prop) %>% sum()
azi


plot_df <- data.frame(
  method = c("azimuth", "celldex", "clustering"),
  perc = c(azi, celldex_prop,0.964)*100, 
  what = "cd3"
)
plot_df %>%
  ggplot(aes(x = method, y = perc, fill = what)) + 
  geom_bar(stat = "identity", color = "black", position="dodge", width = 0.7) +
  scale_y_continuous(expand = c(0,0), limits = c(0,105)) + pretty_plot(fontsize = 7) + L_border() +
  scale_fill_manual(values = c("firebrick")) + theme(legend.position = "none") + 
  labs(x = "annotation method", y = "% T cells") -> p0
cowplot::ggsave2(p0, file = "../plots/enrich_freq.pdf", width = 1.4, height = 1.3)
