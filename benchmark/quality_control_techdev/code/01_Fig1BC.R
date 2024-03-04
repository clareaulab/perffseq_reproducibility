library(bioanalyzeR)
library(dplyr)
library(BuenColors)

rb <- read.bioanalyzer("../xml_data/210915TA_High Sensitivity DNA Assay_DE54107210_2022-09-15_10-16-59.xml")
rbd <- rb$data[,c("length", "fluorescence")] %>%
  mutate(what = paste0(rb$data$sample.index , " - ", rb$samples$sample.name[rb$data$sample.index]))
u_s <- unique(rbd$what)
rbd_0 <- data.frame(
  length = rep(c(20, 34, 10380, 15000), length(u_s)),
  fluorescence = rep(0, 4*length(u_s)),
  what = rep(u_s, each = 4)
)

p1 <- ggplot(rbind(rbd_0, rbd) %>%
               filter(what %in% c("3 - CD20", "4 - PBMC")), aes(x = length, y = fluorescence)) +
  geom_line(color = "firebrick") + facet_wrap(~what, scales = "free_y", nrow = 1) +
  scale_x_log10(breaks = c(35, 75, 150, 300, 1000, 3000, 10380)) +  pretty_plot(fontsize = 6)+
  scale_y_continuous(expand = c(0,0)) + L_border()

cowplot::ggsave2(p1, file = "../plots_for_paper/initial_experiment.pdf", width = 3.5, height = 1.1)


### panel c/d
data_df0 <- data.frame(
  Cond = c("A", "B"),
  internal = c("20220930_PBMC_cr_20220930_PBMC",  "20220930_CD20_cr_20220930_CD20"),
  conf = c(84, 6.3),
  half = c(8.9, 55.5)
)

pB0 <- ggplot(data_df0 %>% reshape2::melt(id.vars = c("Cond", "internal")),
              aes(x = Cond, y = value, fill = variable)) +
  geom_bar(stat = "identity", position=position_dodge(), color = "black", width = 0.5) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "", y = "% of molecules") +
  pretty_plot(fontsize = 6) + L_border() + scale_fill_manual(values = c("blue", "lightgrey"))

cowplot::ggsave2(pB0, file = "../plots_for_paper/initial_experiment_qc.pdf", width = 2, height = 0.9)

