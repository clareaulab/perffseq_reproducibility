library(bioanalyzeR)
library(dplyr)
library(BuenColors)

make_one <- function(sample){
  rb <- read.bioanalyzer(sample)
  rbd <- rb$data %>%
    mutate(what = paste0(rb$data$sample.index , " - ", rb$samples$sample.name[rb$data$sample.index]))
  
  p1 <- ggplot(rbd, aes(x = length, y = fluorescence)) +
    geom_line(color = "firebrick") + facet_wrap(~what, scales = "free_y") +
    scale_x_log10() + ggtitle(sample) + pretty_plot() 
  
  cowplot::ggsave2(p1, file = paste0("output/", sample, ".png"), width = 8, height = 6)
}
lapply(list.files(".", pattern = "xml$"), make_one)
