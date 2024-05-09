library(Seurat)
library(dplyr)
library(BuenColors)
library(data.table)

samples <- gsub(".filteredcounts.h5", "", list.files("../data/counts/"))

samples <-c ("S2453_neg_bottom", "S2453_pos")

lapply(samples, function(dir_base){
  
  mat <- Read10X_h5(paste0("../data/counts/", dir_base, ".filteredcounts.h5"))
  colnames(mat) <- paste0(dir_base, colnames(mat))
  mat
}) %>% do.call("cbind",.) -> big_mat

lapply(samples, function(dir_base){
  dt <- readRDS(paste0("../output/azimuth_pbmc_",dir_base,".rds"))
  dt$barcode <- paste0(dir_base, rownames(dt))
  dt
}) %>% rbindlist() %>% data.frame() -> mdf
rownames(mdf) <- mdf$barcode
common <- intersect(mdf$barcode, colnames(big_mat))

if(FALSE){
  p1 <- data.frame(
    gene = rownames(big_mat),
    meanUMI = rowMeans(big_mat)
  ) %>% arrange(desc(meanUMI)) %>%
    mutate(msy = gene %in% msy_genes) %>%
    mutate(rank = 1:n()) %>%
    arrange(msy) %>%
    ggplot(aes(x = rank, y = meanUMI, color = msy)) +
    geom_point() + scale_color_manual(values = c("lightgrey", "firebrick")) +
    pretty_plot(fontsize = 8) + L_border() +
    theme(legend.position = "none") + labs(x = "Rank ordered gene expression", y = "Mean UMI expression (MSY+ & MSY-)")
  cowplot::ggsave2(p1, file = "../plots/rank_knee.pdf", width = 1.8, height = 1.8)
  
}

# Create seurat object and do seurat things
so <- CreateSeuratObject(
  counts = big_mat, meta.data = mdf
)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 &
               nCount_RNA > 1000 )
dim(so)
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% FindNeighbors()
so <- so %>% FindClusters(resolution = 0.3)

DimPlot(so, group.by = c("seurat_clusters", "sort2"), shuffle = TRUE, 
        split.by = "sort2")
DimPlot(so, group.by = c("seurat_clusters"), label = TRUE, split.by = "sort2")
table(so$name)

fm0 <- FindMarkers(so, ident.1 = "neg", ident.2 = "pos", group.by = "sort2", subset.ident = "0", logfc.threshold = 0.1)
fm0

pm <- presto::wilcoxauc(so, group_by = "sort2")
pm %>% filter( group == "neg") %>%
  mutate(color = feature %in% msy_genes) %>%
  arrange((color)) %>%
  ggplot(aes(y = logFC, x = pct_out, color = color)) +
  geom_point() + scale_color_manual(values = c("lightgrey", "firebrick")) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") + labs(x = "% expressed in MSY+", y = "log(MSY-/MSY+)") -> p2
cowplot::ggsave2(p2, file = "../plots/MA.pdf", width = 1.8, height = 1.8)


so@meta.data %>% group_by(predicted.celltype.l2, sort2) %>%
  summarize(count = n()) %>% ungroup() %>% group_by(sort2) %>%
  mutate(prop = count / sum(count)*100) %>%
  ggplot(aes(x = sort2, y = prop, fill = predicted.celltype.l2)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") + 
  scale_fill_manual(values = (c(jdb_palette("corona"), "salmon"))) +
  pretty_plot(fontsize = 6) + L_border() + scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") -> pX
cowplot::ggsave2(pX, file = "../plots/stacked_bar.pdf", width = 0.7, height = 1.3)

px <- DimPlot(so, group.by = c("predicted.celltype.l2"), shuffle = TRUE, 
        split.by = "sort2", label = FALSE) +
  scale_color_manual(values = (c(jdb_palette("corona"), "salmon"))) &
  theme_void() & theme(legend.position = "none")
cowplot::ggsave2(px, file = "../plots/umap_split.png", width = 4*1, height = 2*1, dpi = 400)

# make plots of cell type proportions
ct <- sort(unique(so$predicted.celltype.l2))
cols_vec <- jdb_palette("corona")[1:26]
names(cols_vec) <- ct

# Reimport QC for
stat_df %>%
  filter(donor == "y51") %>%
  group_by(predicted.celltype.l2, sort2) %>%
  summarize(count = n()) %>%
  ungroup() %>% group_by(sort2) %>%
  ungroup() %>% reshape2::dcast(predicted.celltype.l2 ~ sort2, value.var = "count", fill = (0)) %>%
  mutate(neg_prop = neg/sum(neg), pos_prop = pos/sum(pos)) -> counts_df 

counts_df %>%  
  filter(neg >= 5) %>% mutate(rat = neg_prop / pos_prop) %>% arrange(desc(rat)) %>%
  ggplot(aes(x = neg_prop*100, y = log2(rat), label = predicted.celltype.l2, color = predicted.celltype.l2)) + 
  geom_point(size = 0.5) + scale_x_log10() +  
  pretty_plot(fontsize = 8) + L_border() +
  geom_hline(yintercept = 0, linetype = 2) + scale_color_manual(values = cols_vec) +
  geom_text_repel(size = 2) + labs(x = "% MSY- cells", y = "log2 MSY-/MSY+ proportions") +
  theme(legend.position = "none") -> p1a

counts_df %>% 
  filter(neg < 1) %>%
  arrange(desc(pos_prop)) %>%
  mutate(rank = 1:n())  %>%
  ggplot(aes(x = rank, y = pos_prop*100, color = predicted.celltype.l2, label = predicted.celltype.l2)) +
  geom_point(size = 0.5) +  geom_text_repel(size = 2) +
  labs(x = "Rank ordered cell types unobserved in MSY-", y = "% cell type abundance in MSY+") +
  theme(legend.position = "none") + scale_color_manual(values = cols_vec) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") -> p2a

cowplot::ggsave2(cowplot::plot_grid(p1a, p2a, nrow = 1),
                 file = '../plots/dots.pdf', width = 2.8, height = 1.4)

library(msigdbr)
library(fgsea)
m_df<- msigdbr(species = "Homo sapiens", category = "H")
unique(m_df$gs_name)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fm0$feature <- rownames(fm0)
cluster0.genes<- fm0 %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(feature, avg_log2FC)

vec <- cluster0.genes$avg_log2FC; names(vec) <- cluster0.genes$feature
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = vec)

p1 <- fgseaRes %>% arrange((padj)) %>%
  ggplot(aes(x = ES, y = -log10(padj))) + 
  geom_point() + labs(x = "Enrichment score", y = "-log10 padj") +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p1, file = "../plots/volcano_gsea.pdf", width = 2, height = 1.4)


