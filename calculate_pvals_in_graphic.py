library(ggpubr)
library(rstatix)

stat_pvalue <- df %>%
  rstatix::t_test(x ~ y) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position() %>%
  mutate(y.position = seq(min(y.position), max(y.position),length.out = n())) 

p <- p + ggpubr::stat_pvalue_manual(stat_pvalue, label = 'p.signif', inherit.aes = F)
