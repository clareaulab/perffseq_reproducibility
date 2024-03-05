library(SingleR)
library(celldex)
library(Seurat)
library(dplyr)

# Enrich based on azimuth
lapply(list.files("../output/", full.names = TRUE, pattern = "*projected*"), readRDS) %>%
  rbindlist() %>%
  group_by(predicted.celltype.l2, name) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(name) %>% mutate(prop = count/ sum(count)) -> prop_df

ggplot(prop_df, aes(x = name, y = prop, fill = predicted.celltype.l2)) + 
  geom_bar(stat = "identity") 

lapply(list.files("../output/", full.names = TRUE, pattern = "*projected*"), readRDS) %>%
  rbindlist() %>%
  ggplot(aes(x = refUMAP_1, y = refUMAP_2, color = name)) +
  geom_point()
