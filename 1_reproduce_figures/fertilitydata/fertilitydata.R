library(fda)
library(tidyverse)
library(plot.matrix)
library(abind)
library(robustmatrix)
library(countrycode)
library(gridExtra)
library(Polychrome)
library(covsep)

setwd(this.path::this.dir())

source("../../functions/multivariate_functional_pca_function.R")
source("../../functions/multivariate_functional_data_functions.R")

# Age-specific fertility rates (ASFR) from the Human Fertility Database in long format.
fertilitydata <- read.csv("data/asfrRR.txt", sep = "", skip = 2)
str(fertilitydata)


# Dimension definitions: years (rows), ages 15-45 (columns), countries/regions (observations)
row_dim <- sort(unique(fertilitydata$Year))
col_dim <- 15:45
obs_dim <- unique(fertilitydata$Code)
p <- length(row_dim)
q <- length(col_dim)
n <- length(obs_dim)

# Reshape long data into a 3D array with dimensions (year, age, country)
fertilitydata_wide <- fertilitydata %>%
  filter(Age %in% col_dim) %>%
  pivot_wider(names_from = Year, values_from = ASFR) %>%
  dplyr::select(all_of(c("Code", "Age", row_dim))) %>%
  group_by(Code) %>%
  group_split()

fertilitydata_list <- lapply(fertilitydata_wide, function(x) as.matrix(x)[,-c(1,2)])

X <- aperm(abind(fertilitydata_list, along = 3), c(2,1,3))
X <- X[,,order(obs_dim)]
dimnames(X) <- list(row_dim, col_dim, sort(obs_dim))
class(X) <- "numeric"
str(X)

# Aggregate the years into twelve consecutive 5-year groups
year_agg <- list(70:74,75:79,80:84,85:89,90:94,95:99,100:104,105:109,110:114,115:119,120:124,125:129)
year_groups <- sapply(year_agg, function(sub) paste(range(row_dim[sub]), collapse = ":"))
# Average ASFR within each 5-year window and count how many missing values were aggregated.
X_grouped <- abind(apply(X,3,function(x) t(sapply(year_agg, function(sub) colMeans(x[sub,], na.rm = TRUE))), simplify = "array"),along = 3)
numb_miss_val <- apply(abind(apply(X,3,function(x) t(sapply(year_agg, function(sub) colSums(is.na(x[sub,])))), simplify = "array"),along = 3),3,max)
# Keep only countries with full coverage across all year groups
X_grouped <- X_grouped[,,names(numb_miss_val[numb_miss_val == 0])]

dimnames(X_grouped)[[1]] <- year_groups

X_good <- X_grouped
dimnames(X_good)
apply(X_good,3,function(x) sum(is.na(x)))
X_sub <- X_good[,,apply(X_good,3,function(x) sum(is.na(x))) == 0]
nsub <- dim(X_sub)[[3]]
nsub


# Long-format data frame used for the raw ASFR plot
X_sub_long <- bind_rows(
  apply(X_sub, 1, function(x){
    x_long <- data.frame(x) %>%
      rownames_to_column(var = "Age") %>%
      pivot_longer(-Age) %>%
      mutate(Age = as.numeric(Age), name = factor(name))
    x_long
  }), .id = "Year")

# One distinct colour per country/region
color_vec <- Polychrome::palette36.colors(n = nsub)
names(color_vec) <- NULL
ggplot(X_sub_long, aes(x = Age, y = value, group = name, color = name)) +
  geom_line() + 
  facet_wrap(~Year) + 
  scale_color_manual(values = color_vec)

############################################################################################
############################################################################################
############################################################################################
#### FDA
############################################################################################
############################################################################################
############################################################################################

# Unpenalised B-spline smoothing on the log scale (ASFR is strictly positive).
n_basis <- 10
basis <- create.bspline.basis(rangeval = range(as.numeric(dimnames(X_sub)[[2]])), nbasis = n_basis, norder = 4)

X_smooth_basis_log <- smooth.basis(argvals = as.numeric(dimnames(X_sub)[[2]]),
                                   y = log(aperm(X_sub, perm = c(2,1,3))), fdParobj = basis)
# Exponentiate the basis coefficients to obtain non-negative coefficients on the original scale.
# Downstream MMLE/MMCD/Shapley/PCA all operate on these exponentiated coefficients.
X_smooth_basis <- X_smooth_basis_log
X_smooth_basis$fd$coefs <- exp(X_smooth_basis$fd$coefs)
str(coef(X_smooth_basis))
coeffs <- coef(X_smooth_basis)
X_smoothed <- fd(coef = coef(X_smooth_basis), basisobj = basis)
# Sanity-check plot for the first country: log-scale fit (left) and original-scale fit (right) vs. data
par(mfrow = c(1,2))
plot.fd(fd(coef = (coef(X_smooth_basis_log)[,,1]), basisobj = basis), lwd = 2, lty = 2)
matlines(x = as.numeric(dimnames(X_sub)[[2]]), log(t(X_sub[,,1])), col = "black", lwd = 1, lty = 3)
plot.fd(fd(coef = (coef(X_smooth_basis)[,,1]), basisobj = basis), lwd = 2, lty = 2)
matlines(x = as.numeric(dimnames(X_sub)[[2]]), (t(X_sub[,,1])), col = "black", lwd = 1, lty = 3)
par(mfrow = c(1,1))

# Evaluate the smoothed curves on a fine age grid for plotting
eval_seq <- seq(from = 15, to = 45, by = 0.3)
X_smoothed_eval <- eval.fd(evalarg = eval_seq, fdobj = X_smoothed)
dimnames(X_smoothed_eval)[[1]] <- eval_seq
dimnames(X_smoothed_eval)[[2]] <- dimnames(X_sub)[[1]] 
dimnames(X_smoothed_eval)[[3]] <- dimnames(X_sub)[[3]]
X_smoothed_long <- bind_rows(
  apply(X_smoothed_eval, 2, function(x){
    x_long <- data.frame(x) %>% 
      rownames_to_column(var = "Year") %>% 
      pivot_longer(-Year) %>%
      mutate(Year = as.numeric(Year), name = factor(name))
    x_long
  }), .id = "Age")  

p_smooth <- ggplot(X_smoothed_long, aes(x = Year, y = value, group = name, color = name)) +
  geom_line() + 
  facet_wrap(~Age) + 
  scale_color_manual(values = color_vec) + 
  labs(x = "Age", y = "ASFR", color = "Country/Region") + 
  theme_bw()
p_smooth
ggsave(filename = "plots/p_smooth.pdf", plot = p_smooth, height = 4, width = 9)



# Matrix-variate location/scatter on the basis coefficients: classical MMLE and robust MMCD.
par_fd_mmle <- mmle(coeffs)
set.seed(1)
par_fd_mmcd <- mmcd(coeffs, nthreads = -1, nsamp = 5000, alpha = 0.5)

#############################################
# Test separability of the covariance (year vs. age) on raw and smoothed evaluations
#############################################

library(covsep)
Z <- eval.fd(as.numeric(dimnames(X_sub)[[2]]), fdobj = X_smoothed)

covsep::clt_test(Data = aperm(X_sub, c(3,1,2)), 1:12, 1:31) |> round(5)
covsep::clt_test(Data = aperm(Z, c(3,2,1)), 1:12, 1:31) |> round(5)

#############################################
# end test separability
#############################################

# Classical (MMLE) matrix Mahalanobis distance and the chi-square 0.99 cutoff
md_fd <- md.mfd(X_smooth_basis, par = par_fd_mmle)
quant_fd <- sqrt(qchisq(0.99, prod(dim(coeffs)[1:2])))

# Compare to two competing functional outlier detectors and print the flagged countries
tmp_depth <- mrfDepth::mfd(x = aperm(X_sub,c(2,3,1)), type = "projdepth", diagnostic = TRUE)
dimnames(X)[[3]][which(rowSums(tmp_depth$locOutlX)!=0)]
tmp_ms <- fdaoutlier::msplot(aperm(X_sub,c(3,2,1)), n_projections = 2000, plot = FALSE)
dimnames(X)[[3]][tmp_ms$outliers]

#######################################
# QQ plot of the robust squared Mahalanobis distances against chi^2 with df = 12 * n_basis
#######################################
qq_all <- ggplot(data = NULL, aes(sample = par_fd_mmcd$md)) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  geom_qq(distribution = stats::qchisq, dparams = list("df" = 12 * n_basis)) +
  labs(
    # title = "QQ Plot: Robust Mahalanobis Distance vs Chi-Square",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_classic()
qq_all
ggsave(filename = "plots/qq_plot.pdf", qq_all, height = 5, width = 5)

# Kolmogorov-Smirnov goodness-of-fit against chi^2 for the full sample and the reweighted h-subset
ks.test(par_fd_mmcd$md, "pchisq", df = 12 * n_basis)
ks.test(par_fd_mmcd$md[par_fd_mmcd$h_subset_reweighted], "pchisq", df = 12 * n_basis)
#######################################
#######################################

# Distance-distance plot (classical fMMD vs. robust fMMD) with outliers labelled
library(ggrepel)
outlier_ind_fd <- sqrt(par_fd_mmcd$md) > quant_fd
p_dd_label <- ggplot(data = NULL) +
  geom_hline(yintercept = quant_fd, lty = 2) + 
  geom_vline(xintercept = quant_fd, lty = 2) + 
  geom_point(aes(x = sqrt(md_fd), y = sqrt(par_fd_mmcd$md), label = dimnames(coeffs)[[3]], 
                 color = ifelse(sqrt(par_fd_mmcd$md) > quant_fd, "outlier","regular"))) + 
  geom_label_repel(aes(x = sqrt(md_fd)[outlier_ind_fd], y = sqrt(par_fd_mmcd$md)[outlier_ind_fd], 
                       label = dimnames(coeffs)[[3]][outlier_ind_fd]), nudge_y = 5, max.overlaps = 40) + 
  theme_classic() + 
  scale_color_manual(values = c("#B2182B", "black"), guide = "none") + 
  labs(x = "fMMD", y = "Robust fMMD") 
p_dd_label
ggsave(filename = "plots/p_dd_label.pdf", plot = p_dd_label, height = 4, width = 4)

#################################################################################
##### Cell-wise Shapley decomposition of the MMCD distance (12 year groups x 10 age intervals)
#################################################################################

number_of_pc_for_kernel <- NULL
shv_cell2 <- shapley.mfd(fdobj = X_smooth_basis, par = par_fd_mmcd,
                         time_intervals = 10, sum_coordinates = FALSE, n_pc_kernel = number_of_pc_for_kernel)$shv
# Interval labels look like "15 to 18", "18 to 21", ...; keep the left endpoint as a numeric x-axis value
dimnames(shv_cell2)[[2]] <- sub(" .*", "", dimnames(shv_cell2)[[2]])
shv_cell_long2 <- bind_rows(
  apply(shv_cell2, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "Year") %>%
      pivot_longer(-Year, names_to = "Age", values_to = "shv") %>%
      mutate(shv = shv,
             Age = as.numeric(Age))
    x_long
  }), .id = "Country")

# Evaluate the smoothed curves and the robust-centred curves on the fine age grid for joining with Shapley values
X_eval2 <- eval.fd(evalarg = eval_seq, fdobj = X_smoothed)
dimnames(X_eval2)[[1]] <- eval_seq
dimnames(X_eval2)[[2]] <- dimnames(coeffs)[[2]]
dimnames(X_eval2)[[3]] <- dimnames(coeffs)[[3]]
X_eval2_long <-  bind_rows(
  apply(X_eval2, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "Age") %>%
      pivot_longer(-Age, names_to = "Year", values_to = "asfr") %>% 
      mutate(Age = as.numeric(Age))
    x_long
  }), .id = "Country")

X_center_smoothed_eval <- eval.fd(evalarg = eval_seq, fd(sweep(x = coeffs, MARGIN = 1:2, STATS = par_fd_mmcd$mu, FUN = "-"), basis))
dimnames(X_center_smoothed_eval) <- dimnames(X_eval2)
X_center_smoothed_long <- bind_rows(
  apply(X_center_smoothed_eval, 2, function(x){
    x_long <- data.frame(x) %>%
      rownames_to_column(var = "Age") %>%
      pivot_longer(-Age, names_to = "Country", values_to = "asfr_centered") %>%
      mutate(Age = as.numeric(Age))
    x_long
  }), .id = "Year")%>% mutate(sign = sign(asfr_centered))

# Join curves, Shapley values and centred curves and carry each interval value forward in age
X_shapley_and_eval <- left_join(left_join(X_eval2_long, shv_cell_long2 %>% mutate(Age_group = Age)),X_center_smoothed_long) %>%
  group_by(Country, Year) %>%
  fill(shv, Age_group, .direction = "down")


# Build age-group labels of the form "15 to 18" for the row-aggregated heatmap
age_seq <- seq(from = min(X_shapley_and_eval$Age), to = max(X_shapley_and_eval$Age), length.out = 11)
age_seq_names <- paste0(age_seq[-length(age_seq)], " to ", age_seq[-1])
# Aggregate Shapley values per (Age_group, Country): sum the contributions, zero out non-outlier countries,
# then clip to [-0.2, 0.2] (and threshold below 0.05) so the colour scale stays interpretable.
X_shapley_row_and_eval <- X_shapley_and_eval %>%
  group_by(Age_group, Country) %>% 
  reframe(shv = sum(shv), sign = sign(mean(asfr_centered))) %>% 
  group_by(Country) %>% 
  mutate(shv_prop = ifelse(Country %in% dimnames(X_sub)[[3]][par_fd_mmcd$h_subset_reweighted], 0, shv/sum(shv)),
         shv_prop = ifelse(shv_prop > 0.05, ifelse(shv_prop > 0.2, 0.2*sign, abs(shv_prop)*sign),0)) %>% 
  rename(Age = Age_group) %>%
  mutate(Age = factor(Age, labels = age_seq_names))


p_shapley_row <- ggplot(data = X_shapley_row_and_eval, aes(x = Age, y = Country, fill = shv_prop)) +
  geom_tile(color = "gray") +
  scale_fill_gradient2(low = "#2166AC", mid = "transparent", high = "#B2182B", midpoint = 0, na.value = "transparent", guide = "none") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  scale_x_discrete(expand = c(0,0))+
  labs(x = "Age", y = "Country/Region")

# Same aggregation but per (Year, Country) for the column-aggregated heatmap
X_shapley_col_and_eval <- X_shapley_and_eval %>%
  group_by(Year, Country) %>%
  reframe(shv = sum(shv), sign = sign(mean(asfr_centered))) %>%
  group_by(Country) %>% 
  mutate(shv_prop = ifelse(Country %in% dimnames(X_sub)[[3]][par_fd_mmcd$h_subset_reweighted], 0, shv/sum(shv)),
         shv_prop = ifelse(shv_prop > 0.05, ifelse(shv_prop > 0.2, 0.2*sign, abs(shv_prop)*sign),0))

p_shapley_col <- ggplot(data = X_shapley_col_and_eval, aes(x = Year, y = Country, fill = shv_prop)) +
  geom_tile(color = "gray") +
  scale_fill_gradient2(low = "#2166AC", mid = "transparent", high = "#B2182B", midpoint = 0, na.value = "transparent", guide = "none") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  scale_x_discrete(expand = c(0,0))+
  labs(x = "Year", y = "Country/Region")

p_shapley_row_and_col <- grid.arrange(p_shapley_row, p_shapley_col, nrow = 1)

ggsave(filename = "plots/p_shapley_row_and_col.pdf", p_shapley_row_and_col, height = 4, width = 9)

# Per-country age-by-year raster of ASFR overlaid with cell-wise Shapley contributions.
library(ggnewscale)
p_raster1 <- X_shapley_and_eval %>%
  group_by(Country) %>%
  mutate(shv_prop = shv/sum(shv)) %>%
  filter(Country %in% c("BEL", "NOR", "POL", "USA")) %>%
  mutate(shv_prop = if_else(shv_prop > 0.01, if_else(shv_prop > 0.02, 0.02*sign, shv_prop*sign), 0)) %>%
  ggplot() +
  geom_raster(aes(x = Age, y = Year, fill = asfr), interpolate = TRUE) + 
  scale_fill_gradient(low = "white", high = "black", name = "ASFR") + 
  theme_classic() + 
  scale_x_continuous(expand = c(0,0)) + 
  facet_wrap(~Country)

p_raster2 <- p_raster1 + 
  new_scale("fill") +
  new_scale("alpha") +
  geom_raster(aes(x = Age, y = Year, fill = shv_prop))+
  scale_fill_gradient2(low = alpha("#2166AC", 0.9), mid = "transparent", high = alpha("#B2182B",0.9), midpoint = 0, guide=FALSE)

p_raster2
ggsave(filename = "plots/p_shapley_cell.pdf", p_raster2, height = 6, width = 9)

####################################################################
####################################################################
##### Functional PCA on the kernel covariance (MMLE vs. MMCD)
####################################################################
####################################################################

n_pc_kernel <- 3
n_pc_var <- NULL

# Classical PCA (uses MMLE-based mean and covariance internally) and robust PCA driven by MMCD
X_pca <- pca.mfd(fdobj = X_smoothed, n_pc_kernel = n_pc_kernel, n_pc_var = n_pc_var)
X_pca_rob <- pca.mfd(fdobj = X_smoothed, n_pc_kernel = n_pc_kernel, n_pc_var = n_pc_var, par = par_fd_mmcd)

# Helper to convert an fd object into a long-format data frame for ggplot
fd_to_ggplot <- function(evalarg = NULL, fdobj, colnames = NULL){
  if(is.null(evalarg)){
    evalarg = seq(from = fdobj$basis$rangeval[1], to = fdobj$basis$rangeval[2], length = 100)
  }
  data_wide <- data.frame(eval.fd(evalarg = evalarg, fdobj = fdobj))
  if(!is.null(colnames)){
    colnames(data_wide) <- colnames
  }
  rownames(data_wide) <- evalarg
  data_long <- data_wide %>% 
    rownames_to_column(var = "x") %>% 
    pivot_longer(-x) %>%
    mutate(x = as.numeric(x))
  data_long
}


# Prepend the (age-averaged) mean curve to the PC coefficients so it can be plotted alongside the PCs
X_pca_kernel_with_mean <- X_pca$pc_kernel
X_pca_kernel_with_mean$coefs <- cbind("mean" = rowMeans(X_pca$par$mu), X_pca_kernel_with_mean$coefs)

X_pca_rob_kernel_with_mean <- X_pca_rob$pc_kernel
X_pca_rob_kernel_with_mean$coefs <- cbind("mean" = rowMeans(X_pca_rob$par$mu), X_pca_rob_kernel_with_mean$coefs)

pc_labels <- c("Mean" , paste0("PC ",formatC(1:n_pc_kernel, width = 2, format = "d", flag = "0"), 
                               " \n(MMLE: ", round(diff(c(0,X_pca$explained_variance$kernel[1:n_pc_kernel])),3),
                               " | MMCD: ", round(diff(c(0,X_pca_rob$explained_variance$kernel[1:n_pc_kernel])),3),
                               ")"))

# PC signs are arbitrary; the flips below align them for visual interpretation.
# Reproducibility is provided by the set.seed(1) before mmcd above.
X_pca_kernel_with_mean$coefs[,2] <- -X_pca_kernel_with_mean$coefs[,2]
X_pca_kernel_with_mean$coefs[,3] <- -X_pca_kernel_with_mean$coefs[,3]
X_pca_kernel_with_mean$coefs[,4] <- -X_pca_kernel_with_mean$coefs[,4]
pc_long <- fd_to_ggplot(evalarg = NULL,
                        fdobj = X_pca_kernel_with_mean,
                        colnames = pc_labels)

X_pca_rob_kernel_with_mean$coefs[,2] <- -X_pca_rob_kernel_with_mean$coefs[,2]
X_pca_rob_kernel_with_mean$coefs[,3] <- -X_pca_rob_kernel_with_mean$coefs[,3]
pc_long_rob <- fd_to_ggplot(evalarg = NULL,
                            fdobj = X_pca_rob_kernel_with_mean,
                            colnames = pc_labels)

# Combine classical and robust PCs into a single tidy frame for the side-by-side plot
plt_data_pc <- rbind(pc_long %>% mutate(method = "Classic"), pc_long_rob %>% mutate(method = "Robust")) %>%
  mutate(method = factor(method, levels = c("Classic", "Robust")))

p_pc_MMLE_vs_MMCD <- ggplot() +
  geom_line(data = plt_data_pc,
            aes(x = x, y = value, color = method, linetype = method), linewidth = 1) + 
  facet_wrap(~name, scales = "free_y", nrow = 2) + 
  labs(x = "Age", y = "ASFR") + 
  theme_classic() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  scale_color_manual(values = c("#0F65FF", "#FFA90F"))
p_pc_MMLE_vs_MMCD
ggsave(filename = "plots/p_pc_MMLE_vs_MMCD.pdf", plot = p_pc_MMLE_vs_MMCD, height = 4.6, width = 4.6)
