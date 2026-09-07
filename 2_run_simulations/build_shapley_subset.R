# =============================================================================
# build_shapley_subset.R
#
# Builds shapley_simulation_res_subset.RData (object res_plot), which
# shapley_simulation_summary.R loads, from shapley_simulation_res.RData
# (object res, ~70.6M rows) produced by shapley_simulation.R.
#
# The recipe previously existed only as a commented block at the top of the
# summary script. Only the subset is needed to draw the figures; the full result
# is kept because it is the only way to rebuild it.
# =============================================================================

suppressMessages(library(tidyverse))
require(this.path)
setwd(this.path::this.dir())

load("raw_results/shapley_simulation_res.RData")

res_filterd <- res %>%
  filter(method == "mmcd_fda" | n_basis == 10) %>%
  filter(!absolute_score) %>%
  mutate(method = ifelse(method == "mmcd_fda", paste0(method, "_", n_basis), method))
rm(res); gc()

method_labels <- c("Abs. deviation",
                   "Shapley ",
                   "fDO",              # dprojdepth
                   "Shapley MMCD 10",
                   "Shapley MMCD 20",
                   "Shapley MMCD 30",
                   "MS",
                   "fSDO",             # projdepth
                   "fAO")              # sprojdepth
method_labels_reordered <- method_labels[c(1, 2, 4, 5, 6, 8, 9, 3, 7)]

res_plot <- res_filterd %>%
  mutate(method = factor(method, labels = method_labels),
         method = fct_relevel(method, method_labels_reordered),
         n_interval = as.numeric(ifelse(n_interval == "all", 100,
                              ifelse(n_interval == "none", 0, n_interval))))

save(res_plot, file = "../1_reproduce_figures/stored_results/shapley_simulation_res_subset.RData", compress = "xz")
cat(sprintf("written shapley_simulation_res_subset.RData: %d rows, %.1f MB\n",
            nrow(res_plot), file.size("../1_reproduce_figures/stored_results/shapley_simulation_res_subset.RData") / 1e6))
