library(bioanalyzeR)
library(dplyr)
library(BuenColors)

rb <- read.bioanalyzer("../xml_data/240303-meril-nostrip.xml")
rbd <- rb$data[,c("length", "fluorescence")] %>%
  mutate(what = paste0(rb$data$sample.index , " - ", rb$samples$sample.name[rb$data$sample.index]))
u_s <- unique(rbd$what)
rbd_0 <- data.frame(
  length = rep(c(20, 34, 10380, 15000), length(u_s)),
  fluorescence = rep(0, 4*length(u_s)),
  what = rep(u_s, each = 4)
)

p1 <- ggplot(rbind(rbd_0, rbd) %>%
               filter(what %in% c("1 - sample 1")), aes(x = length, y = fluorescence)) +
  geom_line(color = "firebrick") + 
  scale_x_log10(breaks = c(35, 75, 150, 300, 1000, 3000, 10380)) +  pretty_plot(fontsize = 6)+
  scale_y_continuous(expand = c(0,0)) + L_border()

cowplot::ggsave2(p1, file = "../plots_for_paper/meril_final_corrob.pdf", width = 2, height = 1.1)



##################################

rb1 <- read.bioanalyzer("../xml_data/230301TA_High Sensitivity DNA Assay_DE54107210_2023-03-01_16-46-26.xml")
rbd1 <- rb1$data[,c("length", "fluorescence")] %>%
  mutate(what = paste0(rb1$data$sample.index , " - ", rb1$samples$sample.name[rb1$data$sample.index]))

rb2 <- read.bioanalyzer("../xml_data/230314TA_High Sensitivity DNA Assay_DE54107210_2023-03-13_14-34-59.xml")
rbd2 <- rb2$data[,c("length", "fluorescence")] %>%
  mutate(what = paste0(rb2$data$sample.index , " - ", rb2$samples$sample.name[rb2$data$sample.index]))

rb3 <- read.bioanalyzer("../xml_data/230405_High Sensitivity DNA Assay_DE54107210_2023-04-05_14-53-19.xml")
rbd3 <- rb3$data[,c("length", "fluorescence")] %>%
  mutate(what = paste0(rb3$data$sample.index , " - ", rb3$samples$sample.name[rb3$data$sample.index]))

rbd <- rbind(
  rbd1 %>% filter(what == "9 - CTD"),
  rbd2 %>% filter(what == "2 - sample 2"),
  rbd3 %>% filter(what == "2 - CD3")
)

u_s <- unique(rbd$what)
rbd_0 <- data.frame(
  length = rep(c(20, 34, 10382, 15000), length(u_s)),
  fluorescence = rep(0, 4*length(u_s)),
  what = rep(u_s, each = 4)
)

p1 <- ggplot(rbind(rbd_0, rbd) %>%
               filter(what == "9 - CTD"), aes(x = length, y = fluorescence)) +
  geom_line(color = "firebrick") + 
  scale_x_log10(breaks = c(35, 75, 150, 300, 1000, 3000, 10380)) +  pretty_plot(fontsize = 6)+
  scale_y_continuous(expand = c(0,0)) + L_border()


p2 <- ggplot(rbind(rbd_0, rbd) %>%
               filter(what == "2 - sample 2"), aes(x = length, y = fluorescence)) +
  geom_line(color = "firebrick") + 
  scale_x_log10(breaks = c(35, 75, 150, 300, 1000, 3000, 10380)) +  pretty_plot(fontsize = 6)+
  scale_y_continuous(expand = c(0,0)) + L_border()


p3 <- ggplot(rbind(rbd_0, rbd) %>%
               filter(what == "2 - CD3"), aes(x = length, y = fluorescence)) +
  geom_line(color = "firebrick") + 
  scale_x_log10(breaks = c(35, 75, 150, 300, 1000, 3000, 10380)) +  pretty_plot(fontsize = 6)+
  scale_y_continuous(expand = c(0,0)) + L_border()

al <- cowplot::align_plots(p1, p2, p3, align = "v")

cowplot::ggsave2(cowplot::plot_grid(al[[1]], al[[2]], al[[3]], nrow = 3, ncol = 1), 
                 file = "../plots_for_paper/stack3.pdf", width = 1.8, height = 2.2)
 