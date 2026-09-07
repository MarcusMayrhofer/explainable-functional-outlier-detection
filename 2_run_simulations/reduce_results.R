# =============================================================================
# reduce_results.R
#
# Shrinks the stored result files so they can be distributed, by
#   (1) dropping columns that no figure script uses, and
#   (2) rounding the remaining numeric columns to 6 significant digits.
#
# Rewrites the .RData files in place. Neither step changes any figure.
# =============================================================================

require(this.path)
setwd(this.path::this.dir())

# Columns that appear in no figure script: the covariance errors before the
# multivariate (Kronecker) version, and the errors relative to the clean-data
# benchmark. They remain reproducible from the simulation scripts.
DROP <- c(
  "FN",                                     # implied by TP and recall
  "score_cov_matrix",                       # row-covariance error, superseded by
  "score_cov_function",                     #   score_cov_function_multi
  "score_cov_function_k5", "score_cov_function_k10", "score_cov_function_k20",
  "score_mu_rel", "score_cov_matrix_rel",   # errors relative to the ML estimate
  "score_cov_function_rel", "score_cov_function_rel_k5",
  "score_cov_function_rel_k10", "score_cov_function_rel_k20",
  "score_score_cov_function_multi_rel", "score_score_cov_function_multi_rel_k5",
  "score_score_cov_function_multi_rel_k10", "score_score_cov_function_multi_rel_k20")

DIGITS <- 6
# design columns are left untouched, so a setting label can never be perturbed
DESIGN <- c("method", "n", "eps", "eps_coord_vec", "run", "p", "q", "n_basis",
            "cov_matrix", "cov_function", "eigen_id", "eigen_shift", "nu",
            "kappa", "lambda", "df", "kron_sum_components", "outlier_shift",
            "out_numb", "type", "absolute_score", "n_interval")

DIRS <- c("raw_results", "../1_reproduce_figures/stored_results")
files <- sort(unlist(lapply(DIRS, function(d)
  list.files(d, pattern = "[.]RData$", full.names = TRUE))))

for (f in files) {
  cat(f, "\n")
  e <- new.env(); load(f, envir = e); objs <- ls(e)
  for (nm in objs) {
    x <- get(nm, envir = e); if (!is.data.frame(x)) next
    x <- x[, setdiff(names(x), DROP), drop = FALSE]
    num <- setdiff(names(x)[vapply(x, is.double, logical(1))], DESIGN)
    for (k in num) x[[k]] <- signif(x[[k]], DIGITS)
    assign(nm, x, envir = e)
  }
  save(list = objs, envir = e, file = f, compress = "xz")
  rm(e); gc()
}
