library(pROC)
library(robustbase)  
library(robustmatrix)
library(fda)
library(latex2exp)
library(tidyverse)
library(covsep)

setwd(this.path::this.dir())

source("../../functions/multivariate_functional_pca_function.R")
source("../../functions/multivariate_functional_data_functions.R")

# Setup of Multivariate Functional Isolation Forest (Python implementation accessed via reticulate)
##################################
require(reticulate)
mfif_path <- "../../MFIF.py"
np <- import("numpy")
gc_py <- import("gc")
source_python(mfif_path)
###################################

# Resistance spot welding data: column 1 holds the outlier label, the remaining columns
# are the 5 dynamic-resistance curves concatenated end-to-end (5 * 150 = 750 values per sample).
Xtot <- read.table("data/X_RSW.txt")
labels <- Xtot[, 1]
X_raw <- Xtot[, -1]


# Quick sanity check: all curves in black, labelled outliers in colour
matplot(t(X_raw), type = "l", col = "black", lty = 1)
matlines(t(X_raw[labels == 1,]), type = "l", col = 2:100, lty = 1)

n <- dim(X_raw)[1]
p <- dim(X_raw)[2]

# Reshape into a 3D array with dimensions (5 welding spots, 150 time points, n samples)
X <- aperm(array(t(X_raw), dim = c(p/5,5,n)),c(2,1,3))

dims <- dim(X)
df <- expand_grid(coordinate = 1:dims[1], time = 1:dims[2],sample = 1:dims[3]) %>%
  mutate(value = as.vector(as.matrix((X_raw))))  %>%
  mutate(coordinate = paste(coordinate = "Welding spot", coordinate),
         # Observations 75 and 108 are known outliers and are highlighted in the plots.
         highlight = ifelse(sample %in% c(75,108), "yes", "no"))


# B-spline smoothing of each coordinate curve (n_basis basis functions, no roughness penalty).
# The resulting coefficient array has dimensions (n_basis, 5 coordinates, n samples).
n_basis <- 30
val_seq <- 1:dim(X)[2]
val_range <- range(val_seq)
basis <- create.bspline.basis(rangeval = val_range, nbasis = n_basis, norder = 4)
fdParobj = fdPar(fdobj=basis, Lfdobj=0, lambda=0.0)
X_smooth_basis <- smooth.basis(argvals = val_seq, y = aperm(X, perm = c(2,1,3)), fdParobj = fdParobj)
str(coef(X_smooth_basis))
coeffs <- coef(X_smooth_basis)
X_smoothed <- fd(coef = coef(X_smooth_basis), basisobj = basis)

# Matrix-variate location/scatter on the basis coefficients: classical MMLE and robust MMCD
par_fd_mmle <- mmle(coeffs)
set.seed(1)
t_start <- Sys.time()
par_fd_mmcd <- mmcd(coeffs, scale_consistency = "mmd_med", nsamp = 500)
t_mmcd_smooth <- Sys.time() - t_start

# Matrix Mahalanobis distances from each estimator (the MMCD distances are returned directly)
mmd_fd <- mmd(coeffs, par_fd_mmle$mu, par_fd_mmle$cov_row, par_fd_mmle$cov_col)
mmd_fd_rob <- par_fd_mmcd$md

# ROC curves of the distances against the known labels
roc_curve_mmle <- roc(
  labels,
  mmd_fd,
  plot = FALSE,
  print.auc = TRUE,
  main = "",
  quiet = TRUE
)
roc_curve_mmcd <- roc(
  labels,
  mmd_fd_rob,
  plot = FALSE,
  print.auc = TRUE,
  main = "",
  quiet = TRUE
)


##################################################################################################
# Test separability of the covariance (coordinate vs. time) on raw and smoothed evaluations
##################################################################################################

Z <- eval.fd(val_seq, fdobj = X_smoothed)
covsep::clt_test(aperm(X,c(3,1,2)), 1:5, 1:5)
covsep::clt_test(aperm(Z,c(3,1,2)), 1:5, 1:5)

#######################################
# QQ plot of the robust squared Mahalanobis distances against chi^2 with df = 5 * n_basis
#######################################

qq_all <- ggplot(data = NULL, aes(sample = par_fd_mmcd$md)) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") + 
  geom_qq(distribution = stats::qchisq, dparams = list("df" = 5 * n_basis)) +
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_classic()
qq_all
ggsave(filename = "plots/qq_plot.pdf", qq_all, height = 5, width = 5)


# Kolmogorov-Smirnov goodness-of-fit against chi^2 for the full sample and the reweighted h-subset
ks.test(par_fd_mmcd$md, "pchisq", df = n_basis*5)
ks.test(par_fd_mmcd$md[par_fd_mmcd$h_subset_reweighted], "pchisq", df = n_basis*5)

##################################################################################################
# Cell-wise Shapley decomposition of the MMCD distance (5 coordinates x 10 time intervals)
##################################################################################################

shv_cell <- shapley.mfd(fdobj = X_smooth_basis, par = par_fd_mmcd,
                           time_intervals = 10, sum_coordinates = FALSE, n_pc_kernel = NULL)$shv
dimnames(shv_cell)[[1]] <- 1:5
# Interval labels look like "1 to 16", "16 to 31", ...; keep the left endpoint so the
# time axis becomes numeric for downstream joining and plotting.
dimnames(shv_cell)[[2]] <- sub(" .*", "", dimnames(shv_cell)[[2]])
dimnames(shv_cell)[[3]] <- paste("Observation", 1:n)
shv_cell_long <- bind_rows(
  apply(shv_cell, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "var") %>%
      pivot_longer(-var, names_to = "time", values_to = "shv") %>%
      mutate(shv = sign(shv)*sqrt(abs(shv)),
             time = as.numeric(time))
    x_long
  }), .id = "sample")

# Evaluate the smoothed curves and the robust mean curve on the original time grid
eval_arg <- 1:150

mu_eval <- eval.fd(evalarg = eval_arg, fdobj = fd(coef = par_fd_mmcd$mu, basisobj = basis))
X_eval <- eval.fd(evalarg = eval_arg, fdobj = X_smoothed)

dimnames(X_eval)[[3]] <- dimnames(shv_cell)[[3]]

# Curves centred at the robust mean (used as input for the Shapley overlay)
X_eval_center <- array(dim = dim(X_eval), dimnames = dimnames(X_eval))
X_eval_center[] <- apply(X_eval, MARGIN = 3, function(x) x - mu_eval)

X_eval_long <-  bind_rows(
  apply(X_eval, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "time") %>%
      pivot_longer(-time, names_to = "var", values_to = "value") %>% 
      mutate(time = as.numeric(time))
    x_long
  }), .id = "sample")

# Join Shapley values onto the evaluated curves; carry each interval value forward in time.
X_shapley_and_eval <- left_join(X_eval_long, shv_cell_long) %>%
  group_by(sample, var) %>%
  fill(shv, .direction = "down") %>%
  ungroup()

# Normalise each observation's Shapley values to proportions and clip to [0, 0.5] for the colour scale.
X_shapley_and_eval_prop <- X_shapley_and_eval %>%
  group_by(sample) %>%
  mutate(shv_prop = shv/sum(shv)) %>%
  ungroup() %>%
  group_by(var, time) %>%
  mutate(value_mean = mean(value),
         value_center = value - value_mean) %>%
  ungroup() %>%
  mutate(shv_prop = ifelse(shv_prop < 0, 0, ifelse(shv_prop > 0.5, 0.5,shv_prop)),
         label = ifelse(time == max(time), var, NA))


# Curves coloured by Shapley contribution, with the two known outliers (#75 dotted, #108 dashed) overlaid.
plt_shapley_fd <- ggplot() +
  geom_line(data = X_shapley_and_eval_prop,
            mapping = aes(x = time, y = value, group = sample),
            color = "gray", alpha = 0.5) + 
  geom_line(data = X_shapley_and_eval_prop %>% filter(sample %in% c(paste0("Observation ", c(75,108)))),
            mapping = aes(x = time, y = value, group = sample, color = shv_prop),
            alpha = 0.5, size = 2) +
  geom_line(data = X_shapley_and_eval_prop %>% filter(sample %in% c(paste0("Observation ", c(75)))),
            mapping = aes(x = time, y = value, group = sample),
            color = "black", size = 0.5, linetype = "dotted") + 
  geom_line(data = X_shapley_and_eval_prop %>% filter(sample %in% c(paste0("Observation ", c(108)))),
            mapping = aes(x = time, y = value, group = sample),
            color = "black", size = 0.5, linetype = "dashed") +
  scale_color_gradient(low = "transparent", high = "red", na.value = "transparent", guide = "none") + 
  facet_wrap(~ var, nrow = 1) +
  theme_classic() +
  labs(
    x = "",
    y = latex2exp::TeX("m\\Omega")
  )
plt_shapley_fd
ggsave(filename = "plots/welding_curves_and_shapley.pdf", plt_shapley_fd, height = 3, width = 9)


############################################
############################################
############################################

##################################
# Run Multivariate Functional Isolation Forest with three different dictionaries
##################################

# Convert the data and time grid into numpy arrays expected by the Python implementation
X_py <- np$array(
  aperm(X,c(3,1,2)),
  dtype = "float64"
)

time_py <- np$array(
  1:150,
  dtype = "float64"
)

# Dictionary "Self": split nodes using samples from the data itself
MFIF_Self = MFIForest(X_py,
                      time=time_py,
                      ntrees = 100L,
                      D = 'Self',
                      Dsize = 1000L,
                      innerproduct='auto1',
                      alpha = 0.5)
scores_Self_py <- MFIF_Self$compute_paths()

# Dictionary "gaussian_wavelets"
MFIF_gaussian_wavelets = MFIForest(X_py,
                                   time=time_py,
                                   ntrees = 100L,
                                   D = 'gaussian_wavelets',
                                   Dsize = 1000L,
                                   innerproduct='auto1',
                                   alpha = 0.5)
scores_gaussian_wavelets_py <- MFIF_gaussian_wavelets$compute_paths()


# Dictionary "Brownian"
MFIF_Brownian = MFIForest(X_py,
                          time=time_py, 
                          ntrees = 100L, 
                          D = 'Brownian',
                          Dsize = 1000L, 
                          innerproduct='auto1', 
                          alpha = 0.5)
scores_Brownian_py <- MFIF_Brownian$compute_paths()

# ROC curves for each MFIF variant
roc_mfif_self <- roc(
  labels,
  scores_Self_py,
  plot = FALSE,
  print.auc = TRUE,
  main = "MFIF self",
  quiet = TRUE
) 

roc_mfif_gauss <- roc(
  labels,
  scores_gaussian_wavelets_py,
  plot = FALSE,
  print.auc = TRUE,
  main = "MFIF gauss",
  quiet = TRUE
) 

roc_mfif_brown <- roc(
  labels,
  scores_Brownian_py,
  plot = FALSE,
  print.auc = TRUE,
  main = "MFIF brown",
  quiet = TRUE
) 

############################################
# Competing functional outlier scores: fAO/fSDO/fDO from mrfDepth and directional outlyingness (MS) from fdaoutlier
############################################

f_out_fAO <- mrfDepth::fOutl(x = aperm(X,c(2,3,1)), type = "fAO", diagnostic = TRUE)
f_out_fSDO <- mrfDepth::fOutl(x = aperm(X,c(2,3,1)), type = "fSDO", diagnostic = TRUE)
f_out_fDO <- mrfDepth::fOutl(x = aperm(X,c(2,3,1)), type = "fDO", diagnostic = TRUE)
dout <- fdaoutlier::dir_out(aperm(X,c(3,2,1)))
roc_fAO <- roc(
  labels,
  f_out_fAO$fOutlyingnessX,
  plot = FALSE,
  print.auc = TRUE,
  main = "f_out_fAO",
  quiet = TRUE
)
roc_fSDO <- roc(
  labels,
  f_out_fSDO$fOutlyingnessX,
  plot = FALSE,
  print.auc = TRUE,
  main = "f_out_fSDO",
  quiet = TRUE
)
roc_fDO <- roc(
  labels,
  f_out_fDO$fOutlyingnessX,
  plot = FALSE,
  print.auc = TRUE,
  main = "f_out_fDO",
  quiet = TRUE
)
roc_dout <- roc(
  labels,
  dout$distance,
  plot = FALSE,
  print.auc = TRUE,
  main = "dout",
  quiet = TRUE
)


# Collect all ROC curves and their AUCs into a single tidy frame, then order the methods for the legend
roc_df <- rbind(
  data.frame(
    method = "MFIF self",
    specificity = rev(roc_mfif_self$specificities),
    sensitivity = rev(roc_mfif_self$sensitivities),
    auc = as.numeric(auc(roc_mfif_self))
  ),
  data.frame(
    method = "MFIF brown",
    specificity = rev(roc_mfif_brown$specificities),
    sensitivity = rev(roc_mfif_brown$sensitivities),
    auc = as.numeric(auc(roc_mfif_brown))
  ),
  data.frame(
    method = "MFIF gauss",
    specificity = rev(roc_mfif_gauss$specificities),
    sensitivity = rev(roc_mfif_gauss$sensitivities),
    auc = as.numeric(auc(roc_mfif_gauss))
  ),
  data.frame(
    method = "fAO",
    specificity = rev(roc_fAO$specificities),
    sensitivity = rev(roc_fAO$sensitivities),
    auc = as.numeric(auc(roc_fAO))
  ),
  data.frame(
    method = "fSDO",
    specificity = rev(roc_fSDO$specificities),
    sensitivity = rev(roc_fSDO$sensitivities),
    auc = as.numeric(auc(roc_fSDO))
  ),
  data.frame(
    method = "fDO",
    specificity = rev(roc_fDO$specificities),
    sensitivity = rev(roc_fDO$sensitivities),
    auc = as.numeric(auc(roc_fDO))
  ),
  data.frame(
    method = "MS",
    specificity = rev(roc_dout$specificities),
    sensitivity = rev(roc_dout$sensitivities),
    auc = as.numeric(auc(roc_dout))
  ),
  data.frame(
    method = "MMLE smooth",
    specificity = rev(roc_curve_mmle$specificities),
    sensitivity = rev(roc_curve_mmle$sensitivities),
    auc = as.numeric(auc(roc_curve_mmle))
  ),
  data.frame(
    method = "MMCD smooth",
    specificity = rev(roc_curve_mmcd$specificities),
    sensitivity = rev(roc_curve_mmcd$sensitivities),
    auc = as.numeric(auc(roc_curve_mmcd))
  )
) %>% mutate(method = factor(method, levels = c("MMCD smooth", "MMLE smooth", "fSDO", "fAO", "fDO", "MS", "MFIF brown", "MFIF gauss", "MFIF self")), 
             label = factor(paste0(method, "\n(AUC = ", round(auc,2),")")), 
             label = factor(label, levels = levels(label)[c(7,5,8,3,2,1,4,9,6)]))


# Combined ROC plot across all outlier-detection methods
plt_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = label)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_line(size = 1) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Method"
  ) +
  guides(color=guide_legend(
    keywidth=0.3,
    keyheight=0.35,
    default.unit="inch")
  ) +  
  scale_color_brewer(palette = "Set1")
plt_roc
ggsave(filename = "plots/welding_roc.pdf", plt_roc, height = 5, width = 7)

