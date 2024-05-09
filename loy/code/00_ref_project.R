library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)
library(ggrepel)
msy_genes <- c("ZFY", "TBL1Y", "USP9Y", "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "KDM5D", "EIF1AY") # no RPS4Y1 in flex

# import reference
ref_path <- "../../../ps-large-data-files/other_data/pbmc_multimodal.h5seurat"

options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat(ref_path)

# Import
import_project_scRNAseq <- function(dir_base){
  donor <- strsplit(dir_base, "_", 2)[[1]][1]
  data.path <- paste0("../data/counts/", dir_base, ".filteredcounts.h5")
  raw_counts <- Read10X_h5(filename =  data.path)
  
  raw <- CreateSeuratObject(counts = raw_counts, project = "RNA")
  raw <- SCTransform(raw)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = raw,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  projected <- MapQuery(
    anchorset = anchors,
    query = raw,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  df <- data.frame(
    projected@meta.data,
    projected@reductions$ref.umap@cell.embeddings)
  df$name <- dir_base
  df$donorID <- donor
  
  df$donor <- case_when(
    donor == "S2453" ~ "y51",
    donor == "S5499" ~ "y44",
    donor == "S7361" ~ "y50"
  )
  
  df$sort3 <- case_when(
    grepl("neg_bottom", dir_base) ~ "neg_bottom", 
    grepl("pos", dir_base) ~ "pos", 
    grepl("neg_top", dir_base) ~ "neg_top"
  )
  
  df$sort2 <- case_when(
    grepl("neg_bottom", dir_base) ~ "neg", 
    grepl("pos", dir_base) ~ "pos", 
    grepl("neg_top", dir_base) ~ "neg"
  )
  
  # enumerate totals
  df$total_y <- colSums(raw_counts[msy_genes,])
  df <- cbind(df, t(data.matrix(raw_counts[msy_genes,])))
  saveRDS(df, file = paste0("../output/azimuth_pbmc_", dir_base, ".rds"))
  dir_base
}

samples <- gsub(".filteredcounts.h5", "", list.files("../data/counts/"))
lapply(samples, import_project_scRNAseq)

# Make a QC table
lapply(samples, function(s1){
  readRDS(paste0("../output/azimuth_pbmc_", s1,".rds"))
}) %>% rbindlist() %>%
  group_by(name) %>%
  summarize(count = n(), median(nCount_RNA), median(nFeature_RNA))


# plot summarizing the monocyte percentage to justify not looking at neg top anymore
pd <- lapply(samples, function(s1){
  readRDS(paste0("../output/azimuth_pbmc_", s1,".rds"))
}) %>% rbindlist() %>%
  group_by(donor, sort3, predicted.celltype.l1) %>%
  summarize(count = n()) %>% ungroup() %>%
  group_by(donor, sort3) %>%
  summarize(propMono = sum((predicted.celltype.l1 == "Mono")*count)/sum(count)*100) %>%
  ggplot(aes(x = donor, y = propMono, fill = sort3)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = jdb_palette("corona")) +
  scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 8)+ L_border() +
  labs(fill = "population", x = "Donor", y = "% monocytes from azimuth L1")
cowplot::ggsave2(pd, file = "../plots/supp_bargraphs.pdf", width = 2.8, height = 1.9)


lapply(samples, function(s1){
  readRDS(paste0("../output/azimuth_pbmc_", s1,".rds"))
}) %>% rbindlist() %>%
  filter(sort3 != "neg_top") %>%
  mutate(pctY = total_y/nCount_RNA*100) -> stat_df

# count 0 umi cells
stat_df %>%
  group_by(donor, sort2) %>%
  summarize(mean(total_y <=0)*100)

stat_df %>%
  ggplot(aes(color = sort2, x = pctY)) +
  stat_ecdf() + 
  pretty_plot(fontsize = 7) +
  facet_wrap(~donor) +
  coord_cartesian(xlim = c(0, 0.1)) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1)*10) +
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "% UMIs mapping to Y chromosome", y = "Cumulative fraction") +
  theme(legend.position = "none") -> p1
cowplot::ggsave2(p1, file = "../plots/ecdf_MSY.pdf", width = 3.5, height = 1.3)

lapply(unique(stat_df$donor), function(od){
  xd <- stat_df %>% filter(donor == od) 
  kso <- ks.test(xd$pctY[xd$sort2 == "pos"],
                 xd$pctY[xd$sort2 == "neg"])
  data.frame(p = kso$p.value, donor = od)
}) %>% rbindlist() %>% data.frame()


