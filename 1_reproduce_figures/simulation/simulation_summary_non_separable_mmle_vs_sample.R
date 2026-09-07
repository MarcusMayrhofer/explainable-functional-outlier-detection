# =============================================================================
# simulation_summary_non_separable_mmle_vs_sample.R
#
# One supplement figure: the MMLE against the sample covariance on clean
# non-separable data.
#
# IN : ../stored_results/simulation_non_separable_p10.RData
# OUT: plots/non_separable/sample_vs_mmle.pdf
# =============================================================================

library(tidyverse)
library(ggh4x)
library(stringr)
library(knitr)
library(kableExtra)
library(dplyr)
library(latex2exp)
library(ggridges)
library(gridExtra)
require(this.path)
setwd(this.path::this.dir())

load("../stored_results/simulation_non_separable_p10.RData")
res_non_separable <- res %>% mutate_all(.funs = function(x) unname(unlist(x)))

levels(factor(res_non_separable %>% mutate(method = paste0(method,"_",n_basis)) %>% pull(method)))

res_non_separable_filtered <- res_non_separable %>% 
  filter(method %in% c("sample_clean", "mmle_clean"), eps == 0.05, eigen_shift == 15) %>% 
  distinct() %>% 
  mutate(method = ifelse(method == "mmle_clean", "MMLE", "Sample covariance"),
         n = ifelse(n == 1500, "n = 1500", "n = 5000"))


p1 <- ggplot(res_non_separable_filtered, aes(x = factor(kron_sum_components), y = md_cor_spearman, color = method)) + 
  geom_boxplot() + 
  facet_grid(~n) +
  theme_bw() +
  labs(x = "Sum of separable components", y = "Spearman correlation of MD", color = "Method")+
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),  
        panel.grid.minor.x = element_blank()  
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81"), guide = "none") +
  scale_y_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6))

p2 <- ggplot(res_non_separable_filtered, aes(x = factor(kron_sum_components), y = md_error, color = method)) + 
  geom_boxplot() + 
  facet_grid(~n) +
  theme_bw() +
  labs(x = "Sum of separable components", y = "Covariance estimation error", color = "Method")+
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),  
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_y_log10()

p3 <- grid.arrange(p1, p2, nrow = 2, heights = c(0.45,0.55))
ggsave("plots/non_separable/sample_vs_mmle.pdf", p3, height = 6, width = 9)
