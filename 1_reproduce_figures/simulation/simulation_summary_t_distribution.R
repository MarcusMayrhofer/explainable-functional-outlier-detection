# =============================================================================
# simulation_summary_t_distribution.R
#
# Figures for the heavy-tailed (t-distribution) settings at p = 10.
#
# IN : ../stored_results/simulation_res_p10_t_distribution_shift.RData
#      (contamination magnitudes eigen_shift = 1 and 10)
#      ../stored_results/simulation_res_p10_t_distribution_shift_FIF.RData
# OUT: plots/t_distribution/
#
# The method labels cover 13 methods, i.e. all three MFIF dictionaries. The
# stored MFIF file holds Brownian only, so this script needs a re-run of
# 2_run_simulations/simulation_t_distribution_p10_shift_FIF.R first.
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
source("../../functions/nemenyi_function.R")


load("../stored_results/simulation_res_p10_t_distribution_shift.RData")
res_t_shift <- res %>% mutate_all(.funs = function(x) unname(unlist(x)))
load("../stored_results/simulation_res_p10_t_distribution_shift_FIF.RData")
res_t_FIF <- res %>% mutate_all(.funs = function(x) unname(unlist(x)))
res_t <- rbind(res_t_shift, res_t_FIF)

excluded_methods2 <- c("mmle_clean",
                       "actual",
                       "dprojdepth_coef", "dprojdepth_smooth",
                       "ms_plot_coef", "ms_plot_smooth",
                       "projdepth_coef", "projdepth_smooth",
                       "sprojdepth_coef", "sprojdepth_smooth")

levels(factor(res_t %>% filter(!(method %in% excluded_methods2)) %>% mutate(method = paste0(method,"_",n_basis)) %>% pull(method)))

method_labels <- c("fDO", #dprojdepth
                   "MFIF brown",
                   "MFIF gauss",
                   "MFIF self",
                   "MMCD 10",
                   "MMCD 30",
                   "MMCD raw",
                   "MMLE 10",
                   "MMLE 30",
                   "MMLE raw",
                   "MS",
                   "fSDO", #projdepth
                   "fAO")#sprojdepth

method_labels_reordered <- method_labels[c(7,5,6,10,8,9,12,13,1,11,2,3,4)]

res_shift_filtered <- res_t %>% filter(eigen_id == 1) %>% 
  filter(!(method %in% excluded_methods2),
         df %in% c(3,5,8,12,15)) %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered),
         eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                levels = paste0(ceiling(c(0.1,0.5,1)*10)," outlying coords.")),
         df = factor(paste0("df = ", df), levels = c("df = 3", "df = 5", "df = 8", "df = 12", "df = 15")),
         eigen_shift = factor(ifelse(eigen_shift == 1, "Small shift", "Large shift"), levels = c("Small shift", "Large shift"))) %>%
  drop_na()

##############################################################################
# PERFORMANCE MEASURES
##############################################################################

ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = precision, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec + eigen_shift ~ df, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "Precision", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(),  
        panel.grid.minor.x = element_blank()  
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/t_distribution/t_precision.pdf", height = 11, width = 9)

ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = recall, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec + eigen_shift ~ df, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "Recall", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(),  
        panel.grid.minor.x = element_blank()  
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/t_distribution/t_recall.pdf", height = 11, width = 9)

ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec + eigen_shift ~ df, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "F-Score", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(),  
        panel.grid.minor.x = element_blank()  
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/t_distribution/t_fscroe.pdf", height = 11, width = 9)

ggplot(res_shift_filtered,
       aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec + eigen_shift ~ df, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(),  
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(),  
        panel.grid.minor.x = element_blank()  
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
ggsave("plots/t_distribution/t_AUC.pdf", height = 11, width = 10)

ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
       aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.1, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec + eigen_shift ~ df, scales = "free") +
  theme_bw() +
  labs(y = element_blank(), x = "Covariance error", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_fill_manual(values = c("#0F65FF", "#FFA90F", "#00C2A0", "#FF3D81")) +
  scale_x_log10() +
  coord_flip()
ggsave("plots/t_distribution/t_covariance.pdf", height = 11, width = 9)


#############################################################################################################


##############################################################################
# NEMENYI TEST
##############################################################################
cols = c("precision", "recall", "F.score", #"AUC",
         "score_mu", 
         "score_cov_function_multi", "score_cov_function_multi_k5", 
         "score_cov_function_multi_k10", "score_cov_function_multi_k20")
signs <- c(-1,-1,-1,1,1,1,1,1)


test_filtered <- nemenyi_function(res = res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")),
                                  cols = cols,
                                  vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                  signs = signs)
test_filtered_AUC <- nemenyi_function(res = res_shift_filtered,
                                  cols = "AUC",
                                  vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                  signs = c(-1))

test_filtered_df3 <- nemenyi_function(res = res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")) %>% filter(df == "df = 3"),
                                      cols = cols,
                                      vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                      signs = signs)
test_filtered_df3_AUC <- nemenyi_function(res = res_shift_filtered %>% filter(df == "df = 3"),
                                      cols = "AUC",
                                      vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                      signs = c(-1))

test_filtered_df5 <- nemenyi_function(res = res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")) %>% filter(df == "df = 5"),
                                      cols = cols,
                                      vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                      signs = signs)
test_filtered_df5_AUC <- nemenyi_function(res = res_shift_filtered %>% filter(df == "df = 5"),
                                      cols = "AUC",
                                      vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                      signs = c(-1))

test_filtered_df8 <- nemenyi_function(res = res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")) %>% filter(df == "df = 8"),
                                      cols = cols,
                                      vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                      signs = signs)
test_filtered_df8_AUC <- nemenyi_function(res = res_shift_filtered %>% filter(df == "df = 8"),
                                      cols = "AUC",
                                      vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                      signs = c(-1))

test_filtered_df12 <- nemenyi_function(res = res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")) %>% filter(df == "df = 12"),
                                       cols = cols,
                                       vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                       signs = signs)
test_filtered_df12_AUC <- nemenyi_function(res = res_shift_filtered %>% filter(df == "df = 12"),
                                       cols = "AUC",
                                       vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                       signs = c(-1))

test_filtered_df15 <- nemenyi_function(res = res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF gauss", "MFIF self")) %>% filter(df == "df = 15"),
                                       cols = cols,
                                       vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                       signs = signs)
test_filtered_df15_AUC <- nemenyi_function(res = res_shift_filtered %>% filter(df == "df = 15"),
                                       cols = "AUC",
                                       vars = c("eps", "run", "cov_function", "eps_coord_vec", "df", "eigen_shift"),
                                       signs = c(-1))


test_filtered_grouped <- rbind(test_filtered_df3 %>% mutate(df = "df = 3"), 
                               test_filtered_df5 %>% mutate(df = "df = 5"), 
                               test_filtered_df8 %>% mutate(df = "df = 8"), 
                               test_filtered_df12 %>% mutate(df = "df = 12"), 
                               test_filtered_df15 %>% mutate(df = "df = 15"),
                               test_filtered_df3_AUC %>% mutate(df = "df = 3"), 
                               test_filtered_df5_AUC %>% mutate(df = "df = 5"), 
                               test_filtered_df8_AUC %>% mutate(df = "df = 8"), 
                               test_filtered_df12_AUC %>% mutate(df = "df = 12"), 
                               test_filtered_df15_AUC %>% mutate(df = "df = 15")) %>%
  mutate(df = factor(df, levels = c("df = 3", "df = 5", "df = 8", "df = 12", "df = 15")))

highlighted_test_grouped_summary <- test_filtered_grouped  %>%
  group_by(score, df) %>%
  mutate(is_best = ifelse(lower <= min(upper, na.rm = TRUE), TRUE, FALSE)) %>%
  ungroup()

test_summary_highlighted <- left_join(test_filtered_grouped, highlighted_test_grouped_summary) %>%
  filter(!score %in% c("score_mu", "score_cov_function_multi_k5", "score_cov_function_multi_k10", "score_cov_function_multi_k20")) %>%
  mutate(score = factor(score, labels = c("AUC", "F-Score", "Precision", "Recall", "Covariance")),
         method = factor(method, levels = method_labels_reordered))

p_test <- ggplot(test_summary_highlighted, aes(x = method, y = avg_rank, color = is_best)) +
  geom_point(position = position_dodge(0.5), size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.3, position = position_dodge(0.5)) + 
  facet_grid(df ~ score, scales = "free") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05)) + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "darkgray")) + 
  guides(color = "none") + 
  labs(x = "Method", y = "Mean ranks")
p_test
ggsave("plots/t_distribution/t_test.pdf", height = 7, width = 10)

