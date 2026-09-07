# =============================================================================
# shapley_simulation_summary.R
#
# One main-paper figure: Shapley localisation of the outlying coordinates and
# time intervals.
#
# IN : ../stored_results/shapley_simulation_res_subset.RData
#      built by 2_run_simulations/build_shapley_subset.R
# OUT: plots/shapley/shapley_main.pdf
# =============================================================================

library(tidyverse)
library(ggh4x)
library(stringr)
library(knitr)
library(kableExtra)
library(dplyr)
library(latex2exp)
library(ggridges)
library(lvplot)
require(this.path)
setwd(this.path::this.dir())

ggsave <- function(...) ggplot2::ggsave(..., device = cairo_pdf)


load("../stored_results/shapley_simulation_res_subset.RData")

res_plot <- res_plot %>% mutate(outlier_shift = paste0("Outlier shift = ",outlier_shift), 
                                cov_function = paste0("Covariance function = ",cov_function))

ggplot(res_plot %>% 
         filter(type %in% c("coord", "time"), 
                outlier_shift %in% c("Outlier shift = 1"), eps == 0.1, n_interval %in% c(0,5,10,100), !(method %in% c("Shapley MMCD 10","Shapley MMCD 20"))) %>% 
         mutate(type = factor(ifelse(type == "coord", "Coordinate", "Time"))),
       aes(x = AUC, y = method, color = factor(n_interval))) +
  geom_violin(draw_quantiles = c(0.5), quantile.linewidth = 1, scale = "width", bw = 0.03) +
  facet_nested(cov_function ~ type, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "AUC", color = "Number of time intervals")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank()  
  ) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  scale_color_manual(values = c("black", "#0F65FF", "#FFA90F", "#00C2A0")) +
  coord_flip()
ggsave(filename = "plots/shapley/shapley_main.pdf", height = 6, width = 9)
