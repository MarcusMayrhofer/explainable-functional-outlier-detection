# =============================================================================
# merge_results.R
#
# Combines the per-setting simulation results into the two merged files that
# simulation_summary.R loads.
#
# IN  : raw_results/simulation_res_p50_{shift,shape,iso,cov}.RData      (+ _FIF)
#       raw_results/simulation_res_p3_{shift,shape,iso,cov}.RData       (+ _FIF)
#       raw_results/simulation_non_separable_p10.RData                  (+ _FIF)
#       raw_results/simulation_res_p10_t_distribution_shift.RData       (+ _FIF)
# OUT : stored_results/simulation_res_p50.RData      (object res_p50)
#       stored_results/simulation_res_p3.RData       (object res_p3)
#       stored_results/simulation_res_p50_FIF.RData  (object res_p50_FIF)
#       stored_results/simulation_non_separable_p10.RData               (+ _FIF)
#       stored_results/simulation_res_p10_t_distribution_shift.RData    (+ _FIF)
#
# The non-separable and heavy-tailed settings are each produced by a single run,
# so nothing is combined for them; they are flattened and carried over so that
# every simulation script writes to raw_results/ and this is the only script
# that writes to stored_results/.
#
# Run after the simulation_main_*.R scripts; reads and writes this folder.
# =============================================================================

suppressMessages(library(dplyr))
require(this.path)
setwd(this.path::this.dir())

# the shift/shape settings carry the eigenvalue shift as a labelled factor
lab3 <- function(x) mutate(x, eigen_shift = factor(eigen_shift,
                                                   labels = c("low", "med", "high")))

get1 <- function(f) {
  e <- new.env(); load(f, envir = e); x <- get(ls(e)[1], envir = e); rm(e)
  x %>% mutate_all(.funs = function(z) unname(unlist(z)))
}

# ---- p = 50 : shift, shape, iso, cov ---------------------------------------
res_p50 <- rbind(
  lab3(get1("raw_results/simulation_res_p50_shift.RData")),
  lab3(get1("raw_results/simulation_res_p50_shape.RData")),
  get1("raw_results/simulation_res_p50_iso.RData"),
  get1("raw_results/simulation_res_p50_cov.RData"))
save(res_p50, file = "../1_reproduce_figures/stored_results/simulation_res_p50.RData", compress = "xz")
rm(res_p50); gc()

# ---- p = 3 : the four MFIF files first, then the four distance/depth files --
res_p3_main <- rbind(
  lab3(get1("raw_results/simulation_res_p3_shape.RData")),
  lab3(get1("raw_results/simulation_res_p3_shift.RData")),
  get1("raw_results/simulation_res_p3_cov.RData"),
  get1("raw_results/simulation_res_p3_iso.RData"))
res_p3 <- rbind(
  lab3(get1("raw_results/simulation_res_p3_shift_FIF.RData")),
  lab3(get1("raw_results/simulation_res_p3_shape_FIF.RData")),
  get1("raw_results/simulation_res_p3_iso_FIF.RData"),
  get1("raw_results/simulation_res_p3_cov_FIF.RData"),
  res_p3_main)
rm(res_p3_main); gc()
save(res_p3, file = "../1_reproduce_figures/stored_results/simulation_res_p3.RData", compress = "xz")
rm(res_p3); gc()

# ---- p = 50, MFIF : merged, but kept out of simulation_res_p50.RData --------
# MFIF was run 10 times per configuration at p = 50 rather than 100, so it must
# stay separable from the other methods in the rank-based comparison.
res_p50_FIF <- rbind(
  lab3(get1("raw_results/simulation_res_p50_shift_FIF.RData")),
  lab3(get1("raw_results/simulation_res_p50_shape_FIF.RData")),
  get1("raw_results/simulation_res_p50_iso_FIF.RData"),
  get1("raw_results/simulation_res_p50_cov_FIF.RData"))
save(res_p50_FIF, file = "../1_reproduce_figures/stored_results/simulation_res_p50_FIF.RData", compress = "xz")
rm(res_p50_FIF); gc()

# ---- non-separable and heavy-tailed : one file per setting, nothing to merge -
# Each of these is produced by a single run, so they are only flattened and
# carried over, which keeps every simulation script writing to raw_results/ and
# this script the only place that writes to stored_results/.
for (f in c("simulation_non_separable_p10.RData",
            "simulation_non_separable_p10_FIF.RData",
            "simulation_res_p10_t_distribution_shift.RData",
            "simulation_res_p10_t_distribution_shift_FIF.RData")) {
  res <- get1(file.path("raw_results", f))
  save(res, file = file.path("../1_reproduce_figures/stored_results", f), compress = "xz")
  rm(res); gc()
}
