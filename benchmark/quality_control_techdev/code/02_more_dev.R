library(dplyr)
library(BuenColors)

#### Panel e
data_df <- data.frame(
  Cond = c("A", "B", "C", "D"),
  internal = c("unstained_221118",  "HCR_221104", "HYPR_221104","unstained_hairpin_221118"),
  conf = c(70.3, 1.6, 1.2, 76.2),
  half = c(13.9, 51.1, 51.1, 11.9)
)

pB1 <- ggplot(data_df %>% reshape2::melt(id.vars = c("Cond", "internal")),
       aes(x = Cond, y = value, fill = variable)) +
  geom_bar(stat = "identity", position=position_dodge(), color = "black", width = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 90)) + labs(x = "", y = "% of molecules") +
  pretty_plot(fontsize = 7) + L_border() + scale_fill_manual(values = c("blue", "lightgrey")) +
  theme(legend.position = "none")

cowplot::ggsave2(pB1, file = "../plots_for_paper/decompress_qc.pdf", width = 1.8, height = 1.3)

######3
data_dfG <- data.frame(
  Cond = c("E", "F", "G", "H", "I"),
  internal = c("20220930_CD20_cr_20220930_CD20",
               "20221129_EtOH_CD3_cr_2022_11_29_EtOH_CD3", "20221129_EtOH_unstained_cr_2022_11_29_EtOH_unstained",
               "20230307_None-Roche_cr_20230307_None-Roche", "20230307_FISH-Thermo_cr_20230307_FISH-Thermo" ),
  conf = c(6.3, 0.25, 0.46, 73.13, 66.34),
  half = c(55.5, 52.61, 50.7, 11.68,14.35)
)

pBS <- ggplot(data_dfG %>% reshape2::melt(id.vars = c("Cond", "internal")),
              aes(x = Cond, y = value, fill = variable)) +
  geom_bar(stat = "identity", position=position_dodge(), color = "black", width = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 90)) + labs(x = "", y = "% of molecules") +
  pretty_plot(fontsize = 7) + L_border() + scale_fill_manual(values = c("blue", "lightgrey")) +
  theme(legend.position = "none")

cowplot::ggsave2(pBS, file = "../plots_for_paper/stripping_qc.pdf", width = 2.2, height = 1.3)


##### h
data_dfSORT <- data.frame(
  Cond = c("J", "K", "L"),
  internal = c("230320_CD3_sort_cr_230320_CD3_sort",  "230320_None_cr_230320_None", "230320_CD3_sort_cr_230320_CD3_sort"),
  conf = c(33.07, 18.6, 89.9),
  half = c(25.4, 30.1, 4.2)
)

pB1 <- ggplot(data_dfSORT %>% reshape2::melt(id.vars = c("Cond", "internal")),
              aes(x = Cond, y = value, fill = variable)) +
  geom_bar(stat = "identity", position=position_dodge(), color = "black", width = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 90)) + labs(x = "", y = "% of molecules") +
  pretty_plot(fontsize = 7) + L_border() + scale_fill_manual(values = c("blue", "lightgrey")) +
  theme(legend.position = "none")

cowplot::ggsave2(pB1, file = "../plots_for_paper/sortbuffer_qc.pdf", width = 1.4, height = 1.3)
