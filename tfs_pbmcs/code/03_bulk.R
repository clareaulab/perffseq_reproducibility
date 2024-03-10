library(BuenColors)
library(Seurat)
library(dplyr)
library(data.table)
library(forcats)
library(ggbeeswarm)

hdf <- fread("../data/haemopedia-genes.txt")
order_p2 <- hdf %>% arrange(Pop3) %>% pull(Pop2) %>% unique()

p1 <- hdf %>%
  mutate(Pop2 = factor(as.character(Pop2), levels = order_p2)) %>%
  ggplot(aes(x = Pop2, y = BCL11A)) +
  geom_bar(data = hdf %>% group_by(Pop3, Pop2) %>% summarize(BCL11A = mean(BCL11A)) %>% 
             mutate(Pop2 = factor(as.character(Pop2), levels = order_p2)), stat = "identity", aes(fill = Pop3), color = "black", width = 0.7) + 
  pretty_plot(fontsize = 5) + theme(legend.position = "none") + 
  geom_quasirandom(size = 0.5) + coord_flip() + labs(x = "", y = "log2 BCL11A TPM") + L_border() + scale_y_continuous(expand = c(0,0))

p2 <- hdf %>%
  mutate(Pop2 = factor(as.character(Pop2), levels = order_p2)) %>%
  ggplot(aes(x = Pop2, y = IL3RA)) +
  geom_bar(data = hdf %>% group_by(Pop3, Pop2) %>% summarize(IL3RA = mean(IL3RA)) %>% 
             mutate(Pop2 = factor(as.character(Pop2), levels = order_p2)), stat = "identity", aes(fill = Pop3), color = "black", width = 0.7) + 
  pretty_plot(fontsize = 5) + theme(legend.position = "none") + 
  geom_quasirandom(size = 0.5) + coord_flip() + labs(x = "", y = "log2 IL3RA TPM") + L_border() + scale_y_continuous(expand = c(0,0))

cowplot::ggsave2(p1, file = "../output/haemopedia_BCL11A.pdf", width = 1.5, height = 1.25)
cowplot::ggsave2(p2, file = "../output/haemopedia_IL3RA.pdf", width = 1.5, height = 1.25)

