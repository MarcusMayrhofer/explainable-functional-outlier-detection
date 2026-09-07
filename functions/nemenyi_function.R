# =============================================================================
# nemenyi_function.R
#
# Average ranks of the methods per performance measure, with the confidence
# intervals of the Friedman and post-hoc Nemenyi tests (99% level).
#
#   res   a results table with a `method` column
#   cols  the performance measures to rank
#   vars  the design columns identifying one replication
#   signs +1 where larger is better, -1 where smaller is better
#
# Used by simulation_summary.R and simulation_summary_t_distribution.R.
# =============================================================================

nemenyi_function <- function(res, cols, vars, signs){
  ranks <- data.frame()
  for(i in 1:length(cols)){
    tmp <- res %>%
      dplyr::select(method, vars, "score" = cols[i]) %>%
      mutate(score = score^signs[i]) %>%
      pivot_wider(names_from = method, values_from = score) %>%
      dplyr::select(-all_of(vars)) %>%
      tsutils::nemenyi(plottype = "none", conf.level = 0.99)
    tmp2 <- t(rbind("avg_rank" = tmp$means, "lower"= tmp$intervals[1,], "upper"= tmp$intervals[2,]))
    ranks <- rbind(ranks,data.frame("score" = cols[i], rownames_to_column(data.frame(tmp2), var = "method")))
  }
  ranks
}
