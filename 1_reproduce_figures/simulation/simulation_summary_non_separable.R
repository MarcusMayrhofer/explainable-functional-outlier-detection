# =============================================================================
# simulation_summary_non_separable.R
#
# Figures for the non-separable settings at p = 10: Figure 4 of the main paper
# and the four supplement figures.
#
# IN : ../stored_results/simulation_non_separable_p10.RData
#      ../stored_results/simulation_non_separable_p10_FIF.RData
#      (the MFIF file holds all three dictionaries: Brownian, Gaussian
#      wavelets, Self)
# OUT: plots/non_separable/
#
# MFIF gives no outlyingness cut-off, so it is dropped from the three
# supplement figures that report count-based scores.
# =============================================================================


library(tidyverse)
library(ggh4x)
library(stringr)
library(knitr)
library(kableExtra)
library(dplyr)
library(latex2exp)
library(ggridges)
require(this.path)
setwd(this.path::this.dir())
ggsave <- function(...) ggplot2::ggsave(..., device = cairo_pdf)
load("../stored_results/simulation_non_separable_p10.RData")
res_non_separable <- res %>% mutate_all(.funs = function(x) unname(unlist(x)))
load("../stored_results/simulation_non_separable_p10_FIF.RData")
res_non_separable_FIF <- res %>% mutate_all(.funs = function(x) unname(unlist(x)))
res_non_separable <- rbind(res_non_separable, res_non_separable_FIF)
levels(factor(res_non_separable %>% mutate(method = paste0(method,"_",n_basis)) %>% pull(method)))
method_labels <- c("actual raw", 
                   "fDO", #dprojdepth 
                   "MFIF brown",
                   "MFIF gauss",
                   "MFIF self",
                   "MMCD smooth 10", 
                   "MMCD smooth 20", 
                   "MMCD smooth 30", 
                   "MMCD raw", 
                   "MMLE clean raw",
                   "MMLE smooth 10", 
                   "MMLE smooth 20", 
                   "MMLE smooth 30", 
                   "MMLE raw", 
                   "MS", 
                   "fSDO", #projdepth
                   "sample clean raw",
                   "fAO")#sprojdepth
method_labels_reordered <- method_labels[c(9, 6, 7, 8, 14, 11, 12, 13, 16, 18, 2, 15, 17, 10, 1, 3, 4, 5)]
##############################################################################
# FIGURES FOR THE MAIN PAPER
##############################################################################
res_non_separable_filtered <- res_non_separable %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered)) %>%
  filter(!is.na(time),
         !method %in% c("sample clean raw", "MMLE clean raw", "actual raw"),
         n == 1500) %>%  
  mutate(kron_sum_components = paste0("Sum of ",kron_sum_components, " sep. processes"),
         eigen_shift = factor(eigen_shift, labels = c("Small shift", "Large shift")))
ggplot(res_non_separable_filtered %>% filter(n_basis != 20),
       aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eigen_shift ~ kron_sum_components, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  # keep y major grid lines
        panel.grid.minor.y = element_blank(), # remove y minor grid lines
        panel.grid.major.x = element_line(),  # keep x major grid lines
        panel.grid.minor.x = element_blank()  # optional: remove x minor) +
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/non_separable/AUC_non_separable_p10.pdf", height = 6, width = 9)
##############################################################################
# FIGURES FOR THE SUPPLEMENT
##############################################################################
res_non_separable_supplement <- res_non_separable %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered)) %>%
  filter(!is.na(time), 
         !method %in% c("sample clean raw", "MMLE clean raw", "actual raw")) %>%  
  mutate(kron_sum_components = paste0("Sum of ",kron_sum_components, " sep. processes"),
         eigen_shift = factor(eigen_shift, labels = c("Small shift", "Large shift")),
         n = ifelse(n == 1500, "n = 1500", "n = 5000"))
ggplot(res_non_separable_supplement %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = precision, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(n + eigen_shift ~ kron_sum_components, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "Precision", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  # keep y major grid lines
        panel.grid.minor.y = element_blank(), # remove y minor grid lines
        panel.grid.major.x = element_line(),  # keep x major grid lines
        panel.grid.minor.x = element_blank()  # optional: remove x minor) +
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/non_separable/precision_non_separable_p10.pdf", height = 11, width = 9)
ggplot(res_non_separable_supplement %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = recall, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(n + eigen_shift ~ kron_sum_components, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "Recall", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  # keep y major grid lines
        panel.grid.minor.y = element_blank(), # remove y minor grid lines
        panel.grid.major.x = element_line(),  # keep x major grid lines
        panel.grid.minor.x = element_blank()  # optional: remove x minor) +
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/non_separable/recall_non_separable_p10.pdf", height = 11, width = 9)
ggplot(res_non_separable_supplement %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(n + eigen_shift ~ kron_sum_components, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "F-score", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  # keep y major grid lines
        panel.grid.minor.y = element_blank(), # remove y minor grid lines
        panel.grid.major.x = element_line(),  # keep x major grid lines
        panel.grid.minor.x = element_blank()  # optional: remove x minor) +
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/non_separable/fscore_non_separable_p10.pdf", height = 11, width = 9)
ggplot(res_non_separable_supplement,
       aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(n + eigen_shift ~ kron_sum_components, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  # keep y major grid lines
        panel.grid.minor.y = element_blank(), # remove y minor grid lines
        panel.grid.major.x = element_line(),  # keep x major grid lines
        panel.grid.minor.x = element_blank()  # optional: remove x minor) +
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/non_separable/AUC_non_separable_p10_all.pdf", height = 11, width = 9)
