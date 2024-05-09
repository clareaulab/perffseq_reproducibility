library(dplyr)
library(BuenColors)

data.frame(
  age = c(20, 44, 50, 51),
  loy_pct = c(0.6, 0.91, 0.88, 1.91)
) %>% 
  ggplot(aes(x = as.character(age), y = loy_pct)) + 
  geom_bar(stat = "identity", fill = "lightgrey", color = "black") +
  scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Donor age", y = "% MSY- cells") -> p1

cowplot::ggsave2(p1, file = "../plots/msyNeg_pct.pdf", width = 1.3, height = 1.3)