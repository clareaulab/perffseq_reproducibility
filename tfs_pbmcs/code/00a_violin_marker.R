library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
library(SeuratDisk)

# import reference
ref_path <- "../../../ps-large-data-files/other_data/pbmc_multimodal.h5seurat"

options(future.globals.maxSize = 4000 * 1024^2)
reference <- LoadH5Seurat(ref_path)
#reference$BCL11Ascale <- reference@assays$SCT$data["BCL11A",]
#reference$SPI1scale <- reference@assays$SCT$data["SPI1",]

VLN <- VlnPlot(reference, group.by = "celltype.l2", features = c("BCL11A", "SPI1"),
        pt.size=0) & pretty_plot(fontsize = 7) & L_border() & theme(legend.position = "none") &
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
cowplot::ggsave2(VLN, file = "../output/violin_azimuth_2.pdf", width = 7.5, height = 2)
