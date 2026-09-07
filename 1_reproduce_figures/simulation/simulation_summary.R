# =============================================================================
# simulation_summary.R
#
# Figures for the separable Gaussian settings (p = 3 and p = 50), in three
# sections: main paper, computation time, supplement.
#
# IN : ../stored_results/simulation_res_p3.RData
#      ../stored_results/simulation_res_p50.RData
#      ../stored_results/simulation_res_p50_FIF.RData
#      all three built by 2_run_simulations/merge_results.R
# OUT: plots/gaussian/
#
# MFIF is excluded from the rank-based comparison, which needs balanced
# replications: at p = 50 it was run 10 times per configuration, not 100.
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
ggsave <- function(...) ggplot2::ggsave(..., device = cairo_pdf)
source("../../functions/nemenyi_function.R")
load("../stored_results/simulation_res_p3.RData")
load("../stored_results/simulation_res_p50.RData")
load("../stored_results/simulation_res_p50_FIF.RData")
res_p50 <- rbind(res_p50, res_p50_FIF)
res <- rbind(res_p3, res_p50) %>%
    mutate(kappa = factor(kappa, labels = c("low", "med", "high")),
           lambda = factor(lambda, labels = c("low", "med", "high")),
           nu = factor(nu, labels = c("low", "med", "high")))
unique_values <- res %>% dplyr::select(method, n, eps, eps_coord_vec, run, p, q, cov_function, eigen_id, eigen_shift, nu, kappa, lambda) %>% lapply(unique)
unique_values <- lapply(unique_values, function(x) as.character(x[!is.na(x)]))
summary_table <- data.frame(
  "config" = sapply(unique_values, function(x) paste(x, collapse = ", "))
)
summary_table
####################################################################################################################################
####################################################################################################################################
# new code after revision  
####################################################################################################################################
####################################################################################################################################
excluded_methods2 <- c("mmle_clean",
                       "actual",
                       "dprojdepth_coef", "dprojdepth_smooth",
                       "ms_plot_coef", "ms_plot_smooth",
                       "projdepth_coef", "projdepth_smooth",
                       "sprojdepth_coef", "sprojdepth_smooth")
levels(factor(res %>% filter(!(method %in% excluded_methods2)) %>% mutate(method = paste0(method,"_",n_basis)) %>% pull(method)))
method_labels <- c("fDO", #dprojdepth 
                   "MFIF brown",
                   "MFIF gauss",
                   "MFIF self",
                   "MMCD 10", 
                   "MMCD 20", 
                   "MMCD 30", 
                   "MMCD raw", 
                   "MMLE 10", 
                   "MMLE 20", 
                   "MMLE 30", 
                   "MMLE raw", 
                   "MS", 
                   "fSDO", #projdepth
                   "fAO")#sprojdepth
method_labels_reordered <- method_labels[c(8,5,6,7,12,9,10,11,14,15,1,13,2,4,3)]
res_shift_filtered <- res %>% filter(eigen_id == 1) %>% 
  filter(!(method %in% excluded_methods2)) %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered)) %>%
  drop_na()
subs <- rowSums(is.na(res_shift_filtered)) > 0
sum(subs)
res_shape_filtered <- res %>% filter(eigen_id == 10) %>% 
  filter(!(method %in% excluded_methods2)) %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered)) %>%
  drop_na()
subs <- rowSums(is.na(res_shift_filtered)) > 0
sum(subs)
res_cov_filtered <- res %>% filter(!is.na(kappa)) %>% 
  filter(!(method %in% excluded_methods2)) %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("cov_function", "eigen_id", "eigen_shift", "lambda")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered)) %>%
  drop_na()
subs <- rowSums(is.na(res_cov_filtered)) > 0
sum(subs)
res_iso_filtered <- res %>% filter(!is.na(lambda))  %>% 
  filter(!(method %in% excluded_methods2)) %>% 
  mutate(method = paste0(method,"_", n_basis)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("eigen_id", "eigen_shift", "nu", "kappa")) %>% 
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered)) %>%
  drop_na()
subs <- rowSums(is.na(res_iso_filtered)) > 0
sum(subs)
##############################################################################
# FIGURES FOR THE MAIN PAPER
##############################################################################
res_gaussian_main <- rbind(
  res_shift_filtered %>% filter(eigen_shift == "med", cov_function == "M") %>%
    select(-c(eigen_shift, eigen_id, cov_function)) %>% mutate(type = "Shift"),
  res_shape_filtered %>% filter(eigen_shift == "med", cov_function == "M") %>%
    select(-c(eigen_shift, eigen_id, cov_function)) %>% mutate(type = "Shape"),
  res_iso_filtered %>% filter(lambda == "med", cov_function == "M") %>%
    select(-c(lambda, cov_function)) %>% mutate(type = "Isolated"),
  res_cov_filtered %>% filter(kappa == "med", nu == "med") %>%
    select(-c(kappa, nu)) %>% mutate(type = "Covariance-induced")
) %>% mutate(type = factor(type, levels = c("Shift", "Shape", "Isolated", "Covariance-induced"))) %>%
  filter(n == 1000, p == 50, n_basis != 20, eps_coord_vec %in% c(0.5, 1)) %>%
  select(type, method, n, eps, eps_coord_vec, p, q, n_basis, cov_matrix, n_flagged, F.score, AUC, score_cov_function_multi)
plt_main_fscore <- ggplot(res_gaussian_main %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>%
         mutate(eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                       levels = paste0(ceiling(c(0.1,0.5,1)*50)," outlying coords."))),
       aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec ~ type, scales = "free") +
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
plt_main_AUC <- ggplot(res_gaussian_main %>%
                            mutate(eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                                          levels = paste0(ceiling(c(0.1,0.5,1)*50)," outlying coords."))),
                          aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(eps_coord_vec ~ type, scales = "free") +
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
plt_main_cov <- ggplot(res_gaussian_main %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>%
                            mutate(eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                                          levels = paste0(ceiling(c(0.1,0.5,1)*50)," outlying coords."))),
                          aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.05, rel_min_height = 0.001, panel_scaling = FALSE)+
  facet_nested(eps_coord_vec ~ type, scales = "free") +
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
plt_main <- gridExtra::grid.arrange(plt_main_fscore + theme(legend.position = "none"), 
                                    plt_main_AUC + theme(legend.position = "none"), 
                                    plt_main_cov, heights = c(0.9,0.9,1))
ggsave("plots/gaussian/gaussian_main.pdf", plt_main, height = 11, width = 9)

cols = c("precision", "recall", "F.score", "AUC",
         "score_mu", 
         "score_cov_function_multi", "score_cov_function_multi_k5", 
         "score_cov_function_multi_k10", "score_cov_function_multi_k20")
signs <- c(-1,-1,-1,-1,1,1,1,1,1)
test_shift_filtered_p50 <- nemenyi_function(res = res_shift_filtered %>% filter(p == 50, n_basis != 20, !method %in% c("MFIF brown", "MFIF self", "MFIF gauss")),
                                            cols = cols,
                                            vars = c("eps", "n", "run", "p", "eigen_shift", "cov_function", "eps_coord_vec"),
                                            signs = signs)
test_shape_filtered_p50 <- nemenyi_function(res = res_shape_filtered %>% filter(p == 50, n_basis != 20, !method %in% c("MFIF brown", "MFIF self", "MFIF gauss")),
                                            cols = cols,
                                            vars = c("eps", "n", "run", "p", "eigen_shift", "cov_function", "eps_coord_vec"),
                                            signs = signs)
test_cov_filtered_p50 <- nemenyi_function(res = res_cov_filtered %>% filter(p == 50, n_basis != 20, !method %in% c("MFIF brown", "MFIF self", "MFIF gauss")),
                                          cols = cols,
                                          vars = c("eps", "n", "run", "p", "kappa", "nu", "eps_coord_vec"),
                                          signs = signs)
test_iso_filtered_p50 <- nemenyi_function(res = res_iso_filtered %>% filter(p == 50, n_basis != 20,!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")),
                                          cols = cols,
                                          vars = c("eps", "n", "run", "p", "lambda", "cov_function", "eps_coord_vec"),
                                          signs = signs)
test_summary_p50 <- rbind(cbind(type = "Shift", test_shift_filtered_p50),
                          cbind(type = "Shape", test_shape_filtered_p50),
                          cbind(type = "Isolated", test_iso_filtered_p50),
                          cbind(type = "Cov-induced", test_cov_filtered_p50)) %>% 
  mutate(method = factor(method, levels = method_labels_reordered))
# Step 2: Identify the method with the highest median for each facet
highlighted_test_summary_p50 <- test_summary_p50 %>%
  group_by(type, score) %>%
  mutate(is_best = ifelse(lower <= min(upper, na.rm = TRUE), TRUE, FALSE)) %>%
  ungroup()
# Step 3: Join this information back to the original data
test_summary_highlighted_p50 <- left_join(test_summary_p50, highlighted_test_summary_p50) %>%
  filter(!score %in% c("score_mu", "score_cov_function_multi_k5", "score_cov_function_multi_k10", "score_cov_function_multi_k20"))%>%
  mutate(score = factor(score, labels = c("AUC", "F-Score", "Precision", "Recall", "Covariance")),
         method = factor(method, levels = method_labels_reordered),
         type = factor(type, levels = c("Shift", "Shape", "Isolated", "Cov-induced")))
p_test_p50 <- ggplot(test_summary_highlighted_p50, aes(x = method, y = avg_rank, color = is_best)) +
  geom_point(position = position_dodge(0.5), size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.3, position = position_dodge(0.5)) + 
  facet_grid(type ~ score) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05)) + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "darkgray")) + 
  guides(color = "none") + 
  labs(x = "Method", y = "Mean ranks")
p_test_p50
ggsave("plots/gaussian/test_p50.pdf", height = 8, width = 9)

##############################################################################
# FIGURES FOR THE SUPPLEMENT
##############################################################################

## COMPUTATION TIME
res_time <- rbind(
  res_shift_filtered %>% transmute(method, n, p, time),
  res_shape_filtered %>% transmute(method, n, p, time),
  res_iso_filtered %>% transmute(method, n, p, time),
  res_cov_filtered %>% transmute(method, n, p, time)
)

# Custom labeling function
custom_labels <- function(x) {
  ifelse(x < 1, scales::comma_format()(x), as.character(x))
}
p_time <- ggplot(res_time %>% 
                   mutate(p = paste0("p = ",p),
                          n = paste0("n = ",n)), 
                 aes(x = as.numeric(time), y = method)) + 
  geom_boxplot() +
  # geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.1, rel_min_height = 0.001)+
  facet_grid(n ~p, scales = "free") + 
  labs(x = "Computation time in seconds", y = "Method") +
  theme_bw() + 
  scale_x_log10(labels = custom_labels)
p_time
ggsave("plots/gaussian/res_time.pdf", p_time, height = 5, width = 9)
res_time %>% 
  filter(method %in% c("MMCD raw", "MMCD 10", "MMCD 20", "MMCD 30")) %>%
  group_by(p,n,method) %>% 
  summarize(avg = median(time)) %>% 
  ungroup() %>%
  group_by(p, n) %>% 
  mutate(rel_speed = avg/min(avg))

## Shift
plt_precision_shift <- ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
         filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
         mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                cov_function = paste0("covariance = ", cov_function),
                eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                       levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                p = paste0("p = ", p)),
       aes(x = precision, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shift outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Precision", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_precision_shift
ggsave("plots/gaussian/shift/precision_shift.pdf", plt_precision_shift, height = 11, width = 9)
plt_recall_shift <- ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
         filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
         mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                cov_function = paste0("covariance = ", cov_function),
                eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                       levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                p = paste0("p = ", p)),
       aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shift outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Recall", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_recall_shift
ggsave("plots/gaussian/shift/recall_shift.pdf", plt_recall_shift, height = 11, width = 9)
plt_fscore_shift <- ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
         filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
         mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                cov_function = paste0("covariance = ", cov_function),
                eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                       levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                p = paste0("p = ", p)),
       aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shift outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "F-Score", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_fscore_shift
ggsave("plots/gaussian/shift/fscore_shift.pdf", plt_fscore_shift, height = 11, width = 9)
plt_AUC_shift <- ggplot(res_shift_filtered %>% 
                             filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shift outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_AUC_shift
ggsave("plots/gaussian/shift/AUC_shift.pdf", plt_AUC_shift, height = 11, width = 10)
plt_cov_shift <- ggplot(res_shift_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.08, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shift outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Covariance error", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_log10() +
  coord_flip()
plt_cov_shift
ggsave("plots/gaussian/shift/cov_shift.pdf", plt_cov_shift, height = 11, width = 9)
## Shape
#############
plt_precision_shape <- ggplot(res_shape_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shape outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Precision", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_precision_shape
ggsave("plots/gaussian/shape/precision_shape.pdf", plt_precision_shape, height = 11, width = 9)
plt_recall_shape <- ggplot(res_shape_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shape outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Recall", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_recall_shape
ggsave("plots/gaussian/shape/recall_shape.pdf", plt_recall_shape, height = 11, width = 9)
plt_fscore_shape <- ggplot(res_shape_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shape outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "F-Score", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_fscore_shape
ggsave("plots/gaussian/shape/fscore_shape.pdf", plt_fscore_shape, height = 11, width = 9)
plt_AUC_shape <- ggplot(res_shape_filtered %>% 
                          filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                          mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                 cov_function = paste0("covariance = ", cov_function),
                                 eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                        levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                 p = paste0("p = ", p)),
                        aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shape outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_AUC_shape
ggsave("plots/gaussian/shape/AUC_shape.pdf", plt_AUC_shape, height = 11, width = 10)
plt_cov_shape <- ggplot(res_shape_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                          filter(eigen_shift %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                          mutate(eigen_shift = factor(ifelse(eigen_shift == "low", "Small shift", "Large shift"), levels = c("Small shift", "Large shift")),
                                 cov_function = paste0("covariance = ", cov_function),
                                 eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                        levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                 p = paste0("p = ", p)),
                        aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.08, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Shape outliers" + cov_function + eigen_shift, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Covariance error", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_log10() +
  coord_flip()
plt_cov_shape
ggsave("plots/gaussian/shape/cov_shape.pdf", plt_cov_shape, height = 11, width = 9)
## Isolated
#############
plt_precision_iso <- ggplot(res_iso_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(lambda %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(lambda = factor(ifelse(lambda == "low", "Small isolated peaks", "Large isolated peaks"), levels = c("Small isolated peaks", "Large isolated peaks")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Isolated outliers" + cov_function + lambda, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Precision", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_precision_iso
ggsave("plots/gaussian/iso/precision_iso.pdf", plt_precision_iso, height = 11, width = 9)
plt_recall_iso <- ggplot(res_iso_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(lambda %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(lambda = factor(ifelse(lambda == "low", "Small isolated peaks", "Large isolated peaks"), levels = c("Small isolated peaks", "Large isolated peaks")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Isolated outliers" + cov_function + lambda, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Recall", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_recall_iso
ggsave("plots/gaussian/iso/recall_iso.pdf", plt_recall_iso, height = 11, width = 9)
plt_fscore_iso <- ggplot(res_iso_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                             filter(lambda %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                             mutate(lambda = factor(ifelse(lambda == "low", "Small isolated peaks", "Large isolated peaks"), levels = c("Small isolated peaks", "Large isolated peaks")),
                                    cov_function = paste0("covariance = ", cov_function),
                                    eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                           levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                    p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Isolated outliers" + cov_function + lambda, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "F-Score", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_fscore_iso
ggsave("plots/gaussian/iso/fscore_iso.pdf", plt_fscore_iso, height = 11, width = 9)
plt_AUC_iso <- ggplot(res_iso_filtered %>% 
                        filter(lambda %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                        mutate(lambda = factor(ifelse(lambda == "low", "Small isolated peaks", "Large isolated peaks"), levels = c("Small isolated peaks", "Large isolated peaks")),
                                 cov_function = paste0("covariance = ", cov_function),
                                 eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                        levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                 p = paste0("p = ", p)),
                        aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Isolated outliers" + cov_function + lambda, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_AUC_iso
ggsave("plots/gaussian/iso/AUC_iso.pdf", plt_AUC_iso, height = 11, width = 10)
plt_cov_iso <- ggplot(res_iso_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                          filter(lambda %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                          mutate(lambda = factor(ifelse(lambda == "low", "Small isolated peaks", "Large isolated peaks"), levels = c("Small isolated peaks", "Large isolated peaks")),
                                 cov_function = paste0("covariance = ", cov_function),
                                 eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                        levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                 p = paste0("p = ", p)),
                        aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.08, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Isolated outliers" + cov_function + lambda, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Covariance error", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_log10() +
  coord_flip()
plt_cov_iso
ggsave("plots/gaussian/iso/cov_iso.pdf", plt_cov_iso, height = 11, width = 9)
## Covariance
#############
plt_precision_cov <- ggplot(res_cov_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                           filter(kappa %in% c("low", "high"), nu %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                           mutate(kappa = factor(ifelse(kappa == "low", "Small tau", "Large tau"), levels = c("Small tau", "Large tau")),
                                  nu = factor(ifelse(nu == "low", "Small nu", "Large nu"), levels = c("Small nu", "Large nu")),
                                  eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                         levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                  p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Covariance-induced outliers" + kappa + nu, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Precsion", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_precision_cov
ggsave("plots/gaussian/cov/precision_cov.pdf", plt_precision_cov, height = 11, width = 9)
plt_recall_cov <- ggplot(res_cov_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                           filter(kappa %in% c("low", "high"), nu %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                           mutate(kappa = factor(ifelse(kappa == "low", "Small tau", "Large tau"), levels = c("Small tau", "Large tau")),
                                  nu = factor(ifelse(nu == "low", "Small nu", "Large nu"), levels = c("Small nu", "Large nu")),
                                  eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                         levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                  p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Covariance-induced outliers" + kappa + nu, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Recall", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_recall_cov
ggsave("plots/gaussian/cov/recall_cov.pdf", plt_recall_cov, height = 11, width = 9)
plt_fscore_cov <- ggplot(res_cov_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                           filter(kappa %in% c("low", "high"), nu %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                           mutate(kappa = factor(ifelse(kappa == "low", "Small tau", "Large tau"), levels = c("Small tau", "Large tau")),
                                  nu = factor(ifelse(nu == "low", "Small nu", "Large nu"), levels = c("Small nu", "Large nu")),
                                  eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                         levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                                  p = paste0("p = ", p)),
                           aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Covariance-induced outliers" + kappa + nu, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "F-Score", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_fscore_cov
ggsave("plots/gaussian/cov/fscore_cov.pdf", plt_fscore_cov, height = 11, width = 9)
plt_AUC_cov <- ggplot(res_cov_filtered %>% 
                        filter(kappa %in% c("low", "high"), nu %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                        mutate(kappa = factor(ifelse(kappa == "low", "Small tau", "Large tau"), levels = c("Small tau", "Large tau")),
                               nu = factor(ifelse(nu == "low", "Small nu", "Large nu"), levels = c("Small nu", "Large nu")),
                               eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                      levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                               p = paste0("p = ", p)),
                        aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Covariance-induced outliers" + kappa + nu, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "AUC", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_continuous(limits = c(-0.05,1.05), breaks = seq(from = 0, to = 1, length.out = 6)) +
  coord_flip()
plt_AUC_cov
ggsave("plots/gaussian/cov/AUC_cov.pdf", plt_AUC_cov, height = 11, width = 10)
plt_cov_cov <- ggplot(res_cov_filtered %>% filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                        filter(kappa %in% c("low", "high"), nu %in% c("low", "high"), !(n_basis == 20), n == 1000) %>% 
                        mutate(kappa = factor(ifelse(kappa == "low", "Small tau", "Large tau"), levels = c("Small tau", "Large tau")),
                               nu = factor(ifelse(nu == "low", "Small nu", "Large nu"), levels = c("Small nu", "Large nu")),
                               eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."), 
                                                      levels = paste0(ceiling(c(c(0.1,0.5,1)*3, c(0.1,0.5,1)*50))," outlying coords.")),
                               p = paste0("p = ", p)),
                        aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) + 
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.08, rel_min_height = 0.001)+
  facet_nested(p + eps_coord_vec ~ "Covariance-induced outliers" + kappa + nu, scales = "free") + 
  theme_bw() +
  labs(y = element_blank(), x = "Covariance error", fill = "Fraction of contaminated samples", color = "Fraction of contaminated samples")+ 
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05),
        panel.grid.major.y = element_line(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(), 
        panel.grid.minor.x = element_blank() 
  ) +
  scale_color_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_fill_manual(values = c("#00C2A0", "#0F65FF", "#FF3D81", "#FFA90F")) + 
  scale_x_log10() +
  coord_flip()
plt_cov_cov
ggsave("plots/gaussian/cov/cov_cov.pdf", plt_cov_cov, height = 11, width = 9)

## Nemenyi test, p = 3
cols = c("precision", "recall", "F.score",# "AUC",
         "score_mu", 
         "score_cov_function_multi", "score_cov_function_multi_k5", 
         "score_cov_function_multi_k10", "score_cov_function_multi_k20")
signs <- c(-1,-1,-1,1,1,1,1,1)
test_shift_filtered <- nemenyi_function(res = res_shift_filtered  %>% 
                                          filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss"))%>% 
                                          filter(p == 3, n_basis != 20),
                                        cols = cols,
                                        vars = c("eps", "n", "run", "p", "eigen_shift", "cov_function", "eps_coord_vec"),
                                        signs = signs)
test_shift_filtered_AUC <- nemenyi_function(res = res_shift_filtered %>% 
                                          filter(p == 3, n_basis != 20),
                                        cols = "AUC",
                                        vars = c("eps", "n", "run", "p", "eigen_shift", "cov_function", "eps_coord_vec"),
                                        signs = c(-1))
test_shape_filtered <- nemenyi_function(res = res_shape_filtered %>%
                                          filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                                          filter(p == 3, n_basis != 20),
                                        cols = cols,
                                        vars = c("eps", "n", "run", "p", "eigen_shift", "cov_function", "eps_coord_vec"),
                                        signs = signs)
test_shape_filtered_AUC <- nemenyi_function(res = res_shape_filtered %>% 
                                          filter(p == 3, n_basis != 20),
                                        cols = "AUC",
                                        vars = c("eps", "n", "run", "p", "eigen_shift", "cov_function", "eps_coord_vec"),
                                        signs = c(-1))
test_cov_filtered <- nemenyi_function(res = res_cov_filtered %>% 
                                        filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                                        filter(p == 3, n_basis != 20),
                                      cols = cols,
                                      vars = c("eps", "n", "run", "p", "kappa", "nu", "eps_coord_vec"),
                                      signs = signs)
test_cov_filtered_AUC <- nemenyi_function(res = res_cov_filtered %>% 
                                        filter(p == 3, n_basis != 20),
                                      cols = "AUC",
                                      vars = c("eps", "n", "run", "p", "kappa", "nu", "eps_coord_vec"),
                                      signs = c(-1))
test_iso_filtered <- nemenyi_function(res = res_iso_filtered %>% 
                                        filter(!method %in% c("MFIF brown", "MFIF self", "MFIF gauss")) %>% 
                                        filter(p == 3, n_basis != 20),
                                      cols = cols,
                                      vars = c("eps", "n", "run", "p", "lambda", "cov_function", "eps_coord_vec"),
                                      signs = signs)
test_iso_filtered_AUC <- nemenyi_function(res = res_iso_filtered%>% 
                                        filter(p == 3, n_basis != 20),
                                      cols = "AUC",
                                      vars = c("eps", "n", "run", "p", "lambda", "cov_function", "eps_coord_vec"),
                                      signs = c(-1))
test_summary <- rbind(cbind(type = "Shift", test_shift_filtered),
                      cbind(type = "Shape", test_shape_filtered),
                      cbind(type = "Isolated", test_iso_filtered),
                      cbind(type = "Cov-induced", test_cov_filtered),
                      cbind(type = "Shift", test_shift_filtered_AUC),
                      cbind(type = "Shape", test_shape_filtered_AUC),
                      cbind(type = "Isolated", test_iso_filtered_AUC),
                      cbind(type = "Cov-induced", test_cov_filtered_AUC)) %>% 
  mutate(method = factor(method, levels = method_labels_reordered))
# Step 2: Identify the method with the highest median for each facet
highlighted_test_summary <- test_summary %>%
  group_by(type, score) %>%
  mutate(is_best = ifelse(lower <= min(upper, na.rm = TRUE), TRUE, FALSE)) %>%
  ungroup()
# Step 3: Join this information back to the original data
test_summary_highlighted <- left_join(test_summary, highlighted_test_summary) %>%
  filter(!score %in% c("score_mu", "score_cov_function_multi_k5", "score_cov_function_multi_k10", "score_cov_function_multi_k20"))%>%
  mutate(score = factor(score, labels = c("AUC", "F-Score", "Precision", "Recall", "Covariance")),
         method = factor(method, levels = method_labels_reordered),
         type = factor(type, levels = c("Shift", "Shape", "Isolated", "Cov-induced")))
p_test <- ggplot(test_summary_highlighted, aes(x = method, y = avg_rank, color = is_best)) +
  geom_point(position = position_dodge(0.5), size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.3, position = position_dodge(0.5)) + 
  facet_grid(type ~ score, scales = "free") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1.1, vjust = 1.05)) + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "darkgray")) + 
  guides(color = "none") + 
  labs(x = "Method", y = "Mean ranks")
p_test
ggsave("plots/gaussian/test_p3.pdf", height = 9, width = 10.5)

## Depth-based competitors
# MFIF is excluded here: at p = 50 it was run with 10 replications rather than
# 100, so it is not evaluated on the same datasets as the depth-based methods.
res_depth_based <- res %>% filter(!(method %in% c("mmle", "mmcd", "mmle_fda", "mmcd_fda", "actual", "mmle_clean",
                                                 "MFIF_Brownian", "MFIF_gaussian_wavelets", "MFIF_Self"))) %>%
  mutate(method = recode(method, 
                         "dprojdepth" = "fDO", 
                         "dprojdepth_coef" = "fDO coef", 
                         "dprojdepth_smooth" = "fDO smooth", 
                         "ms_plot" = "MS", 
                         "ms_plot_coef" = "MS coef", 
                         "ms_plot_smooth" = "MS smooth", 
                         "projdepth" = "fSDO", 
                         "projdepth_coef" = "fSDO coef", 
                         "projdepth_smooth" = "fSDO smooth", 
                         "sprojdepth" = "fAO", 
                         "sprojdepth_coef" = "fAO coef", 
                         "sprojdepth_smooth" = "fAO smooth"),
         method2 = sub(" .*", "", method),
         method = ifelse(n_basis == "raw", method, paste(method, n_basis)))
res_depth_based_shift_filtered <- res_depth_based %>% filter(eigen_id == 1) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  drop_na()
subs <- rowSums(is.na(res_depth_based_shift_filtered)) > 0
sum(subs)
res_depth_based_shape_filtered <- res_depth_based %>% filter(eigen_id == 10) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("kappa", "nu", "lambda")) %>% 
  drop_na()
subs <- rowSums(is.na(res_depth_based_shift_filtered)) > 0
sum(subs)
res_depth_based_cov_filtered <- res_depth_based %>% filter(!is.na(kappa)) %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("cov_function", "eigen_id", "eigen_shift", "lambda")) %>% 
  drop_na()
subs <- rowSums(is.na(res_depth_based_cov_filtered)) > 0
sum(subs)
res_depth_based_iso_filtered <- res_depth_based %>% filter(!is.na(lambda))  %>% 
  mutate(precision = if_else(TP == 0 & FP == 0, 0, precision)) %>% 
  dplyr::select(-c("eigen_id", "eigen_shift", "nu", "kappa")) %>% 
  drop_na()
subs <- rowSums(is.na(res_depth_based_iso_filtered)) > 0
sum(subs)
res_depth_based_gaussian <- rbind(
  res_depth_based_shift_filtered %>% filter(eigen_shift == "med", cov_function == "M") %>%
    select(-c(eigen_shift, eigen_id, cov_function)) %>% mutate(type = "Shift"),
  res_depth_based_shape_filtered %>% filter(eigen_shift == "med", cov_function == "M") %>%
    select(-c(eigen_shift, eigen_id, cov_function)) %>% mutate(type = "Shape"),
  res_depth_based_iso_filtered %>% filter(lambda == "med", cov_function == "M") %>%
    select(-c(lambda, cov_function)) %>% mutate(type = "Isolated"),
  res_depth_based_cov_filtered %>% filter(kappa == "med", nu == "med") %>%
    select(-c(kappa, nu)) %>% mutate(type = "Covariance-induced")
) %>% mutate(type = factor(type, levels = c("Shift", "Shape", "Isolated", "Covariance-induced"))) %>%
  filter(n == 1000, p == 50, n_basis != 20, eps_coord_vec %in% c(1)) %>%
  select(type, method, method2, n, eps, eps_coord_vec, p, q, n_basis, cov_matrix, n_flagged, F.score, AUC, score_cov_function_multi)
plt_depth_based_fscore <- ggplot(res_depth_based_gaussian %>%
                            mutate(eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                                          levels = paste0(ceiling(c(0.1,0.5,1)*50)," outlying coords."))),
                          aes(x = F.score, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(type ~ method2, scales = "free") +
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
plt_depth_based_fscore
ggsave("plots/gaussian/depth_based/gaussian_depth_based_fscore.pdf", plt_depth_based_fscore, height = 11, width = 9)
plt_depth_based_AUC <- ggplot(res_depth_based_gaussian %>%
                                   mutate(eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                                                 levels = paste0(ceiling(c(0.1,0.5,1)*50)," outlying coords."))),
                                 aes(x = AUC, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.02, rel_min_height = 0.001)+
  facet_nested(type ~ method2, scales = "free") +
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
plt_depth_based_AUC
ggsave("plots/gaussian/depth_based/gaussian_depth_based_AUC.pdf", plt_depth_based_AUC, height = 11, width = 9)
plt_depth_based_cov <- ggplot(res_depth_based_gaussian %>%
                                mutate(eps_coord_vec = factor(paste0(ceiling(eps_coord_vec*p)," outlying coords."),
                                                              levels = paste0(ceiling(c(0.1,0.5,1)*50)," outlying coords."))),
                              aes(x = score_cov_function_multi, y = method, fill = factor(eps), color = factor(eps))) +
  geom_density_ridges(alpha = 0.5, scale = 1, bandwidth = 0.05, rel_min_height = 0.001)+
  facet_nested(type ~ method2, scales = "free") +
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
plt_depth_based_cov
ggsave("plots/gaussian/depth_based/gaussian_depth_based_cov.pdf", plt_depth_based_cov, height = 11, width = 9)

