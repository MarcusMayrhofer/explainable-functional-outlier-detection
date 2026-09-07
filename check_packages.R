# =============================================================================
# check_packages.R
#
# Reports which of the packages used anywhere in this code are missing, and
# prints the install.packages() call for them. Installs nothing itself.
#
#   Rscript check_packages.R
# =============================================================================

pkgs <- c(
  # the method itself
  "robustmatrix", "covsep",
  # infrastructure
  "this.path", "reticulate",
  # data handling and figures
  "tidyverse", "dplyr", "ggplot2", "stringr", "lubridate", "scales",
  "ggh4x", "ggnewscale", "ggrepel", "ggridges", "geomtextpath", "gridExtra",
  "latex2exp", "lvplot", "Polychrome", "plot.matrix", "knitr", "kableExtra",
  # functional data and the competing detectors
  "fda", "fdaoutlier", "mrfDepth", "pROC", "robustbase",
  # simulation
  "abind", "cellWise", "MixMatrix", "rSPDE", "tsutils",
  "foreach", "doParallel", "doSNOW",
  # real-data examples
  "countrycode", "zoo",
  # the El Nino map: rnaturalearthdata is only a *suggested* dependency of
  # rnaturalearth but ne_countries() needs it, so it must be installed too
  "sf", "rnaturalearth", "rnaturalearthdata"
)

have <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)

cat("R", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat(sum(have), "of", length(pkgs), "packages available\n\n")

if (any(!have)) {
  cat("MISSING:\n")
  for (p in pkgs[!have]) cat("  ", p, "\n")
  cat("\ninstall.packages(c(",
      paste0('"', pkgs[!have], '"', collapse = ", "), "))\n\n", sep = "")
} else {
  cat("All R packages present.\n\n")
}

# MFIF runs in Python through reticulate and imports numpy. We only look for an
# interpreter on the PATH.
py <- Sys.which(c("python", "python3"))
py <- py[nzchar(py)]
if (length(py)) {
  cat("Python on PATH:", py[[1]], "\n")
  cat("  MFIF additionally needs numpy in that interpreter.\n")
} else {
  cat("No Python found on PATH.\n")
  cat("  Needed only for the MFIF comparisons (Python + numpy, via reticulate).\n")
}
