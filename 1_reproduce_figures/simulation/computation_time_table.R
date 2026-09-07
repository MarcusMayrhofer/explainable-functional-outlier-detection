# Median computation time per run, for every method and every simulation
# setting. Complements plots/gaussian/res_time.pdf, which covers the Gaussian
# settings only. Prints a plain table and a LaTeX version.

library(dplyr)
library(tidyr)
setwd(this.path::this.dir())
SR <- "../stored_results"

get1 <- function(f) {
  e <- new.env(); load(file.path(SR, f), envir = e)
  x <- get(ls(e)[1], envir = e); rm(e); x
}

# paper labels for the raw-data methods; smoothed variants get " <m>" appended
raw_label <- c(mmcd = "MMCD raw", mmle = "MMLE raw", projdepth = "fSDO",
               sprojdepth = "fAO", dprojdepth = "fDO", ms_plot = "MS",
               MFIF_Brownian = "MFIF brown", MFIF_gaussian_wavelets = "MFIF gauss",
               MFIF_Self = "MFIF self")
smooth_label <- c(mmcd_fda = "MMCD", mmle_fda = "MMLE",
                  projdepth_smooth = "fSDO smooth", projdepth_coef = "fSDO coef",
                  sprojdepth_smooth = "fAO smooth", sprojdepth_coef = "fAO coef",
                  dprojdepth_smooth = "fDO smooth", dprojdepth_coef = "fDO coef",
                  ms_plot_smooth = "MS smooth", ms_plot_coef = "MS coef")

label_of <- function(method, n_basis) {
  m <- as.character(method); b <- as.character(n_basis)
  out <- ifelse(m %in% names(raw_label), raw_label[m],
         ifelse(m %in% names(smooth_label), paste(smooth_label[m], b), NA))
  out
}

collect <- function(files, setting) {
  d <- bind_rows(lapply(files, function(f) {
    x <- get1(f)
    tibble(method = as.character(x$method), n_basis = as.character(x$n_basis),
           n = as.numeric(x$n), time = x$time)
  }))
  d %>% filter(!is.na(time), time > 0,
               !method %in% c("actual", "mmle_clean", "sample_clean")) %>%
    mutate(label = label_of(method, n_basis), setting = setting) %>%
    filter(!is.na(label))
}

all <- bind_rows(
  collect(c("simulation_res_p3.RData"), "Gaussian p=3"),
  collect(c("simulation_res_p50.RData", "simulation_res_p50_FIF.RData"), "Gaussian p=50"),
  collect(c("simulation_non_separable_p10.RData",
            "simulation_non_separable_p10_FIF.RData"), "Non-separable p=10"),
  collect(c("simulation_res_p10_t_distribution_shift.RData",
            "simulation_res_p10_t_distribution_shift_FIF.RData"), "Heavy-tailed p=10"))

# Columns are single (setting, n) configurations. Pooling over n would be
# misleading, because n differs between settings: the non-separable runs use
# n = 1500, 5000 while the Gaussian ones use n = 300, 1000.
all <- all %>% mutate(col = paste0(setting, ", n=", n))
COLS <- c("Gaussian p=3, n=300", "Gaussian p=3, n=1000",
          "Gaussian p=50, n=300", "Gaussian p=50, n=1000",
          "Non-separable p=10, n=1500", "Non-separable p=10, n=5000",
          "Heavy-tailed p=10, n=1000")

tab <- all %>%
  group_by(col, label) %>%
  summarise(med = median(time), .groups = "drop") %>%
  pivot_wider(names_from = col, values_from = med) %>%
  select(label, any_of(COLS)) %>%
  arrange(desc(`Gaussian p=50, n=1000`))

cat("\n=== median seconds per run, per (setting, n) ===\n")
print(as.data.frame(tab %>% mutate(across(where(is.numeric), ~round(., 1)))),
      row.names = FALSE)

# the subset that appears in the manuscript discussion
key <- c("MFIF brown", "MFIF gauss", "MFIF self", "MMCD raw",
         "MMCD 10", "MMCD 20", "MMCD 30", "fDO", "fSDO", "fAO", "MS")
cat("\n=== LaTeX, methods discussed in the text ===\n")
sub <- tab %>% filter(label %in% key) %>%
  mutate(label = factor(label, levels = key)) %>% arrange(label)
fmt <- function(v) ifelse(is.na(v), "--", formatC(v, format = "f", digits = 1))
present <- intersect(COLS, names(sub))
setting_of <- sub("(.*), n=.*", "\\1", present)
n_of <- paste0("$n=", sub(".*, n=", "", present), "$")
grp <- rle(setting_of)                      # group the columns by setting

cat("\\begin{tabular}{l", strrep("r", length(present)), "}\n\\toprule\n", sep = "")
pos <- 1; h1 <- character(0); cm <- character(0)
for (j in seq_along(grp$lengths)) {
  k <- grp$lengths[j]
  nm <- gsub("p=", "$p=", grp$values[j]); nm <- paste0(nm, "$")
  h1 <- c(h1, if (k > 1) sprintf("\\multicolumn{%d}{c}{%s}", k, nm) else nm)
  cm <- c(cm, sprintf("\\cmidrule(lr){%d-%d}", pos + 1, pos + k))
  pos <- pos + k
}
cat(" & ", paste(h1, collapse = " & "), " \\\\\n", sep = "")
cat(paste(cm, collapse = ""), "\n", sep = "")
cat("Method & ", paste(n_of, collapse = " & "), " \\\\\n\\midrule\n", sep = "")
for (i in seq_len(nrow(sub)))
  cat(sprintf("%s & %s \\\\\n", sub$label[i],
              paste(vapply(present, function(cc) fmt(sub[[cc]][i]), ""),
                    collapse = " & ")))
cat("\\bottomrule\n\\end{tabular}\n")

cat("\n=== total CPU time per setting ===\n")
print(as.data.frame(all %>% group_by(setting) %>%
  summarise(cpu_years = round(sum(time) / 3600 / 24 / 365.25, 2), .groups = "drop")),
  row.names = FALSE)
cat(sprintf("\ntotal over all settings: %.2f CPU-years\n",
            sum(all$time) / 3600 / 24 / 365.25))
cat("(the Shapley simulation records no timings and is not included)\n")
