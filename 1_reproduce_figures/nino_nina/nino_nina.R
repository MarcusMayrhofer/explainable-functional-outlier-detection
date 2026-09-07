library(fda)
library(tidyverse)
library(plot.matrix)
library(abind)
library(robustmatrix)
library(gridExtra)
library(Polychrome)
library(zoo)
library(ggh4x)
library(ggrepel)
library(geomtextpath)
library(mrfDepth)
library(fdaoutlier)
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

# To get english date labels in ggplot
Sys.setlocale("LC_ALL", "English")

# Convert a numeric month-of-period index (1..12 starting in June) to a proper date in 2023/2024
# so the x-axis of the plots can use ggplot's date scale.
convert_date <- function(numeric_date, year = 2023) {
  numeric_date2 <- (numeric_date + 5)%%13
  numeric_date2 <- c(numeric_date2[numeric_date2 >= 6], numeric_date2[numeric_date2 < 6] + 1)
  month <- floor(numeric_date2)
  fractional_part <- numeric_date2 - month + 0.001
  years <- ifelse(month %in% 6:12, year, year + 1)
  day <- ceiling(fractional_part * days_in_month(lubridate::ym(paste0(years, "-", month))))
  lubridate::ymd(paste0(year = years, "-", month, "-", day))
}

# Monthly ERSSTv5 SST values and 3-month anomalies for the four Niño regions.
# source: https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii
data <- read.csv("data/data_ERSSTv5.txt", sep = "", skip = 0, head = TRUE)

# For each value, return the highest threshold whose absolute value is exceeded for n consecutive months.
# The sign indicates the direction (positive = El Niño, negative = La Niña).
check_consecutive <- function(values, thresholds, n = 5) {
  labels <- rep(0, length(values))
  for (i in seq_along(thresholds)) {
    condition <- rollapply(values, width = n, FUN = function(x) all(x > thresholds[i]), align = "right", fill = NA)
    condition2 <- rollapply(-values, width = n, FUN = function(x) all(x > thresholds[i]), align = "right", fill = NA)
    labels[condition] <- thresholds[i]
    labels[condition2] <- -thresholds[i]
  }
  labels
}

# Tag each month with: its "El Niño period" label (June year_t to May year_{t+1}),
# a 3-month rolling average of ANOM.3, and an El Niño / La Niña classification (threshold * sign).
# Incomplete leading and trailing periods are dropped.
data_nino_nina <- data %>%
  mutate(
    period = case_when(
      MON >= 6 ~ paste0(YR, ":", YR + 1),
      MON < 6 ~ paste0(YR - 1, ":", YR)
    ),
    date = lubridate::ym(paste0(YR, "-", MON)),
    avg_anom = rollapply(ANOM.3, width = 3, FUN = function(x) round(mean(x),2), fill = NA, align = "right"),
    type = check_consecutive(avg_anom, c(0.5,1,1.5,2))
  )%>%
  filter(!(period %in% c("1949:1950", "2025:2026"))) %>%
  mutate(month_name = factor(lubridate::month(date, label = TRUE, locale = "UCT"),
                             levels = month.abb[c(6:12,1:5)]))

# Per-period summary: number of El Niño / La Niña months, an overall class, and the strongest threshold reached.
data_overview <- data_nino_nina %>%
  group_by(period) %>%
  summarize(nino = sum(type > 0, na.rm = TRUE),
            nina = sum(type < 0, na.rm = TRUE),
            strength = type[which.max(abs(type))]) %>%
  mutate(class = ifelse(nino == 0 & nina == 0, "none", ifelse(nino > 0, "nino", "nina"))) %>%
  relocate(class, .before = strength)

# Variant of data_overview where strength is the largest 3-month anomaly (continuous), used for the
# method-comparison correlations further below.
data_overview2 <- data_nino_nina %>%
  group_by(period) %>%
  summarize(nino = sum(type > 0, na.rm = TRUE),
            nina = sum(type < 0, na.rm = TRUE),
            strength = avg_anom[which.max(abs(avg_anom))]) %>%
  mutate(class = ifelse(nino == 0 & nina == 0, "none", ifelse(nino > 0, "nino", "nina"))) %>%
  relocate(class, .before = strength)

# Wide-format SST values per (period, month, station). Drops anomaly columns; only raw SSTs go into the array.
data_wide <- data %>%
  dplyr::select(-c("ANOM", "ANOM.1", "ANOM.2", "ANOM.3")) %>%
  mutate(
    period = case_when(
      MON >= 6 ~ paste0(YR, ":", YR + 1),
      MON < 6 ~ paste0(YR - 1, ":", YR)
    )
  )%>% 
  filter(!(period %in% c("1949:1950", "2025:2026")))

data_wide_split <- data_wide %>% 
  group_by(period) %>% 
  group_split() 

data_list <- lapply(data_wide_split, function(x) as.matrix(x[,-c(1,2,7)]))

#######################################################################################################################################

# 3D SST array with dimensions (4 Niño regions, 12 months June-May, n periods)
X <- aperm(abind(data_list, along = 3), c(2,1,3))
dimnames(X)[[2]] <- c(6:12,1:5)
dimnames(X)[[3]] <- unique(data_wide$period)
n <- dim(X)[[3]]

############################################################################################
# Competing functional outlier-detection methods (mrfDepth and fdaoutlier)
# Scores are used for the method-comparison correlations below
############################################################################################


ms_outliers <- fdaoutlier::msplot(dts = aperm(X,c(3,2,1)), plot = FALSE)
dimnames(X)[[3]][ms_outliers$outliers]


# Multivariate functional depth and outlyingness measures
Result <- mrfDepth::mfd(x = aperm(X,c(2,3,1)), diagnostic = TRUE)
Plot <- mrfDepth::fHeatmap(rowValues = Result$MFDdepthZ,
                 cellValues = Result$crossdepthZ,
                 type = "depth",
                 legend.title = "HD")

f_out_fAO <- mrfDepth::fOutl(x = aperm(X,c(2,3,1)), type = "fAO", diagnostic = TRUE)
f_out_fSDO <- mrfDepth::fOutl(x = aperm(X,c(2,3,1)), type = "fSDO", diagnostic = TRUE)
f_out_fDO <- mrfDepth::fOutl(x = aperm(X,c(2,3,1)), type = "fDO", diagnostic = TRUE)
dout <- fdaoutlier::dir_out(aperm(X,c(3,2,1)))
grid.arrange(
  mrfDepth::fom(f_out_fAO),
  mrfDepth::fom(f_out_fSDO),
  mrfDepth::fom(f_out_fDO), nrow = 2#,
)

fdaoutlier::functional_boxplot(t(X[4,,]), depth_method = "mbd")
boxplot.fd((X[4,,]))

dimnames(X)[[3]][which(rowSums(f_out_fAO$locOutlX)!=0)]
dimnames(X)[[3]][which(rowSums(f_out_fSDO$locOutlX)!=0)]
dimnames(X)[[3]][which(rowSums(f_out_fDO$locOutlX)!=0)]
dimnames(X)[[3]][ms_outliers$outliers]

##################################
# Run Multivariate Functional Isolation Forest with three different dictionaries
##################################

# Convert the data and time grid into numpy arrays expected by the Python implementation
X_py <- np$array(
  aperm(X,c(3,1,2)),
  dtype = "float64"
)

time_py <- np$array(
  1:12,
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

############################################################################################
## FDA - B-spline smoothing of the SST curves
############################################################################################

n_basis <- 6
basis <- create.bspline.basis(rangeval = range(as.numeric(dimnames(X)[[2]])), nbasis = n_basis, norder = 4)

fdParobj = fdPar(fdobj=basis, Lfdobj=0, lambda=0.0)
X_smooth_basis <- smooth.basis(argvals = 1:12, y = aperm(X, perm = c(2,1,3)), fdParobj = fdParobj)
coeffs <- coef(X_smooth_basis)
X_smoothed <- fd(coef = coef(X_smooth_basis), basisobj = basis)
# Sanity check: smoothed fits vs. raw points for the first three periods
par(mfrow = c(1,3))
plot.fd(fd(coef = (coef(X_smooth_basis)[,,1]), basisobj = basis), lwd = 2, lty = 2)
matlines(t(X[,,1]))
plot.fd(fd(coef = (coef(X_smooth_basis)[,,2]), basisobj = basis), lwd = 2, lty = 2)
matlines(t(X[,,2]))
plot.fd(fd(coef = (coef(X_smooth_basis)[,,3]), basisobj = basis), lwd = 2, lty = 2)
matlines(t(X[,,3]))
par(mfrow = c(1,1))

# Matrix-variate location/scatter on the basis coefficients: classical MMLE and robust MMCD.
# Seed makes the MMCD random subsampling reproducible.
par_fd_mmle <- mmle(coeffs)
set.seed(1)
par_fd_mmcd <- mmcd(coeffs, nthreads = 1, nsamp = 5000, alpha = 0.5, scale_consistency = "mmd_med")

# Classical matrix Mahalanobis distance and the chi-square 0.99 cutoff (used by the DD-plot below)
md_fd <- mmd(coeffs, par_fd_mmle$mu, par_fd_mmle$cov_row, par_fd_mmle$cov_col)

quant_fd <- sqrt(qchisq(0.99, prod(dim(coeffs)[1:2])))

##################################################################################################
# Univariate analysis on the Niño 3.4 region alone (the standard region used to classify ENSO events)
##################################################################################################
coeffs_uni <- t(coeffs[,4,])
center_fd_uni <- colMeans(coeffs_uni)
cov_fd_uni <- cov(coeffs_uni)
# Classical and robust (MCD) Mahalanobis distances + chi-square cutoff for the univariate case
set.seed(1)
par_fd_uni_mcd <- robustbase::covMcd(x = coeffs_uni, nsamp = 5000)
md_fd_uni <- mahalanobis(coeffs_uni, center_fd_uni, cov_fd_uni)
md_rob_fd_uni <- mahalanobis(coeffs_uni, par_fd_uni_mcd$center, par_fd_uni_mcd$cov)
quant_fd_uni <- sqrt(qchisq(0.99, ncol(coeffs_uni)))

# Spearman correlations of each outlier-score with |strongest 3-month anomaly| of the period.
# Used for the method-comparison table in the manuscript; nothing here is saved to disk.
rbind(
  data.frame(method = "MFIF self",cor = cor(abs(pull(data_overview2, strength)), scores_Self_py, method = "spearman")),
  data.frame(method = "MFIF gauss",cor = cor(abs(pull(data_overview2, strength)), scores_gaussian_wavelets_py, method = "spearman")),
  data.frame(method = "MFIF brown",cor = cor(abs(pull(data_overview2, strength)), scores_Brownian_py, method = "spearman")),
  data.frame(method = "fAO",cor = cor(abs(pull(data_overview2, strength)), f_out_fAO$fOutlyingnessX, method = "spearman")),
  data.frame(method = "fSDO",cor = cor(abs(pull(data_overview2, strength)), f_out_fSDO$fOutlyingnessX, method = "spearman")),
  data.frame(method = "fDO",cor = cor(abs(pull(data_overview2, strength)), f_out_fDO$fOutlyingnessX, method = "spearman")),
  data.frame(method = "dout",cor = cor(abs(pull(data_overview2, strength)), dout$distance, method = "spearman")),
  data.frame(method = "md",cor = cor(abs(pull(data_overview2, strength)), md_rob_fd_uni, method = "spearman")),
  data.frame(method = "mmd",cor = cor(abs(pull(data_overview2, strength)), par_fd_mmcd$md, method = "spearman"))
) %>% mutate(cor = round(cor,2))


##################################################################################################
# Test separability of the covariance (region vs. month) on raw and smoothed evaluations
##################################################################################################

library(covsep)
Z <- eval.fd(1:12, fdobj = X_smoothed)
covsep::clt_test(Data = aperm(X, c(3,1,2)), 1:4, 1:12)
covsep::clt_test(Data = aperm(Z, c(3,2,1)), 1:4, 1:12)

##################################################################################################
# QQ plot of the robust squared Mahalanobis distances against chi^2 with df = 4 * n_basis
##################################################################################################

qq_all <- ggplot(data = NULL, aes(sample = par_fd_mmcd$md)) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") + 
  geom_qq(distribution = stats::qchisq, dparams = list("df" = 4 * n_basis)) +
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_classic()
qq_all
ggsave(filename = "plots/plt_qq_all.pdf", plot = qq_all, height = 2.5, width = 2.5)


# Kolmogorov-Smirnov goodness-of-fit against chi^2 for the full sample and the reweighted h-subset
ks.test(par_fd_mmcd$md, "pchisq", df = 4*n_basis)
ks.test(par_fd_mmcd$md[par_fd_mmcd$h_subset_reweighted], "pchisq", df = 4*n_basis)

##################################################################################################
# Univariate vs multivariate robust distances (Niño 3.4 alone vs. all four regions jointly)
##################################################################################################
plt_data_uni_multi <- data_overview %>%
  mutate("fMD_rob" = sqrt(md_rob_fd_uni), 
         "fMMD_rob" = sqrt(par_fd_mmcd$md),
         "uni_out" = (sqrt(md_rob_fd_uni) > quant_fd_uni),
         "multi_out" = (sqrt(par_fd_mmcd$md) > quant_fd),
         "duration" = nino+nina)

plt_uni_vs_multi <- ggplot(plt_data_uni_multi, aes(x = fMD_rob, y = fMMD_rob, label = period, 
                               color = sign(strength), size = duration)) + 
  geom_vline(xintercept = quant_fd_uni) + 
  geom_hline(yintercept = quant_fd) +
  geom_point() +
  geom_point(data = plt_data_uni_multi %>% filter(uni_out | multi_out)) +
  geom_label_repel(data = plt_data_uni_multi %>% filter(uni_out | multi_out), 
                   max.overlaps = 30, force = 5, nudge_y = 0.3, color = "black", size = 3) + 
  scale_color_gradient2(low = "#2166AC", mid = "gray", high = "#B2182B", midpoint = 0, guide = "none") +
  labs(x = "Robust fMD of region 'Niño 3.4'", y = "Robust fMMD",
       color = "Highest\nSST anomaly\nin period", size = "Number of 'El Niño' and\n'La Niña' months in period") + 
  theme_classic() + 
  theme(legend.position = "bottom") + 
  scale_size_continuous(range = c(2,5))
plt_uni_vs_multi
ggsave(filename = "plots/plt_uni_vs_multi.pdf", plot = plt_uni_vs_multi, height = 5, width = 5)

# Evaluate the smoothed curves on a fine month grid (mapped to calendar dates) for the overview plot
eval_seq <- seq(from = 1, to = 12, length = 100)
eval_seq_date <- convert_date(eval_seq)
X_eval <- eval.fd(evalarg = eval_seq, fdobj = X_smoothed) |> aperm(c(2,1,3))
dimnames(X_eval) <- list(dimnames(X)[[1]],
                         as.character(eval_seq_date),
                         #paste0(month(eval_seq_date),"-",day(eval_seq_date)),
                         dimnames(X)[[3]])
X_eval_bind <- bind_rows(
  apply(X_eval, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "station") %>%
      pivot_longer(-station, names_to = "Month", values_to = "value")
    x_long
  }), .id = "period") %>% 
  left_join(data_overview[,c(1,4)])

# Flag periods not in the MMCD reweighted h-subset as outliers; relabel stations for plotting
plt_overview_data <- X_eval_bind %>%
  mutate(outlier = if_else(period %in% dimnames(X)[[3]][par_fd_mmcd$h_subset_reweighted],
                           "Regular",
                           "Outlier"),
         station = factor(station, labels = c("Niño 1+2", "Niño 3", "Niño 3.4", "Niño 4")))

# Robust mean curve (average over the regular periods only)
rob_mean <- plt_overview_data %>%
  filter(outlier == "Regular") %>%
  group_by(station, Month) %>%
  summarise(value = mean(value))

# Per-region overview: regular periods in gray, outliers coloured by El Niño / La Niña / neutral, with the robust mean overlaid
plt_overview_smooth2 <- ggplot() +
  geom_line(data = plt_overview_data %>% filter(outlier == "Regular"), 
            mapping = aes(x = ymd(Month), y = value, group = period), 
            linewidth = 0.3, color = "gray") + 
  geom_line(data = plt_overview_data %>% filter(outlier != "Regular"), 
            mapping = aes(x = ymd(Month), y = value, group = period, color = class), 
            linewidth = 0.3) + 
  geom_line(data = rob_mean, aes(x = ymd(Month), y = value), linewidth = 0.7) + 
  geom_labelpath(data = plt_overview_data %>%
                  filter(outlier != "Regular"),
                mapping = aes(x = ymd(Month), y = value, color = class, group = period, 
                              label = substr(period, 3,4), 
                              hjust = period),
                linewidth = 0.3, size = 1.6, label.padding = 0.02,
                straight = TRUE, show.legend = FALSE) +
  facet_wrap(station~., scales = "free") + 
  scale_color_manual(values = c("#2166AC", "#B2182B", "black"), 
                     labels = c("La Niña", "El Niño", "Neutral\nState")) + 
  labs(x = "Month", 
       y = "Sea surface temperatue", 
       color = "    Type")+
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand= c(0,1)) +
  theme_classic() + 
  theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plt_overview_smooth2
ggsave(filename = "plots/plt_overview_smooth2.pdf", plot = plt_overview_smooth2, height = 5, width = 7)

#################################################################################
##### Cell-wise Shapley decomposition of the MMCD distance (4 regions x 100 time intervals)
#################################################################################
number_of_pc_for_kernel <- NULL
shv_cell2 <- shapley.mfd(fdobj = X_smooth_basis, par = par_fd_mmcd,
                         time_intervals = 100, sum_coordinates = FALSE, n_pc_kernel = number_of_pc_for_kernel)$shv
# Interval labels look like "1 to 1.11", "1.11 to 1.22", ...; keep the left endpoint as a numeric x-axis value
dimnames(shv_cell2)[[2]] <- sub(" .*", "", dimnames(shv_cell2)[[2]])
shv_cell_long2 <- bind_rows(
  apply(shv_cell2, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "Station") %>%
      pivot_longer(-Station, names_to = "Month", values_to = "shv") %>%
      mutate(Month = as.numeric(Month))
    x_long
  }), .id = "Period")

# Evaluate the smoothed curves and the robust mean curve on the same grid as the Shapley intervals
mu_eval <- eval.fd(evalarg = as.numeric(dimnames(shv_cell2)[[2]]), fdobj = fd(coef = par_fd_mmcd$mu, basisobj = basis))
X_eval2 <- eval.fd(evalarg = as.numeric(dimnames(shv_cell2)[[2]]), fdobj = X_smoothed)
dimnames(X_eval2)[[1]] <- dimnames(shv_cell2)[[2]]
dimnames(X_eval2)[[2]] <- dimnames(coeffs)[[2]]
dimnames(X_eval2)[[3]] <- dimnames(coeffs)[[3]]
X_eval2_long <-  bind_rows(
  apply(X_eval2, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "Month") %>%
      pivot_longer(-Month, names_to = "Station", values_to = "SST") %>% 
      mutate(Month = as.numeric(Month))
    x_long
  }), .id = "Period")

# Centred curves (used to attach a sign to each cell's Shapley contribution)
X_eval2_center <- array(dim = dim(X_eval2), dimnames = dimnames(X_eval2))
X_eval2_center[] <- apply(X_eval2, MARGIN = 3, function(x) x - mu_eval)

X_eval2_center_long <-  bind_rows(
  apply(X_eval2_center, 3, function(x){
    x_long <- data.frame(x, check.names = FALSE) %>%
      rownames_to_column(var = "Month") %>%
      pivot_longer(-Month, names_to = "Station", values_to = "SST") %>% 
      mutate(Month = as.numeric(Month))
    x_long
  }), .id = "Period") %>% 
  mutate(sign = sign(SST)) %>% 
  dplyr::select(-SST)

# Join curves, Shapley values, sign of the centred curve, and the El Niño / La Niña / Neutral class per period
X_shapley_and_eval <- left_join(left_join(left_join(X_eval2_long, shv_cell_long2),X_eval2_center_long),
                                data_overview2 %>%
                                  rename(Period = period) %>%
                                  mutate(class = ifelse(class == "nino", "El Niño", ifelse(class == "nina", "La Niña", "Neutral"))) %>%
                                  select(Period, class)) %>%
  group_by(Period) %>%
  mutate(Month_raw = Month, Month = convert_date(Month))

# Faceted overlay: every outlier period, its curves coloured by signed Shapley contribution (clipped at the 90% quantile).
p_shapley_lines_outliers <- ggplot(X_shapley_and_eval %>%
         filter(!(Period %in% dimnames(X)[[3]][par_fd_mmcd$h_subset_reweighted])) %>%
         mutate(shv = ifelse(shv > quantile(shv, 0.90), quantile(shv, 0.90), shv),
                shv = ifelse(shv < 0, 0, shv), 
                shv_sign = shv*sign) ,
       aes(x = Month, y = SST, color = shv_sign, group = Period)) + 
  geom_line(lwd = 0.2, color = "black")  + 
  geom_line(lwd = 1.2) + 
  ggh4x::facet_grid2(class~Station, scales = "free", independent = "all") +  
  scale_color_gradient2(low = "#2166AC", mid = "transparent", high = "#B2182B", 
                        midpoint = 0, na.value = "transparent")  +    
  theme_classic() + 
  labs(x = "Month", 
       y = "Sea surface temperatue")+
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand= c(0,5)) +
  coord_cartesian(clip = 'off') +
  theme(panel.spacing = unit(0.5, "lines"), 
        plot.margin = unit(c(1,2,1,0.5), "lines"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p_shapley_lines_outliers
ggsave(filename = "plots/shapley_lines_outliers.pdf", plot = p_shapley_lines_outliers, height = 5, width = 9.5)



# Same Shapley overlay but restricted to two periods
p_shapley_lines <- ggplot(X_shapley_and_eval %>%
                            group_by(Period) %>% 
                            mutate(shv_prop = shv/sum(shv)) %>%  
                            ungroup() %>% 
                            group_by(Station, Month_raw) %>% 
                            mutate(SST_mean = mean(SST),
                                   SST_center = SST - SST_mean) %>%
                            ungroup() %>%
                            filter(Period %in% c("1956:1957", "2015:2016")) %>%
                            mutate(shv_prop = ifelse(shv_prop > quantile(shv_prop, 0.9), quantile(shv_prop, 0.9), shv_prop),
                                   shv_prop = sign*ifelse(shv_prop < 0, 0, shv_prop),
                                   label = ifelse(Month_raw == max(Month_raw), Station, NA)),
                          aes(x = ymd(Month), y = SST, color = shv_prop, group = Station, label = Station)) + 
  geom_label_repel(aes(label = label), nudge_x = 5, xlim = as.Date(c('2023-06-01','2024-08-01')), color = "black", size = 2,
                   arrow = arrow(angle = 30, length = unit(0.25, "lines"),
                                 ends = "last", type = "open"),
                   box.padding = unit(0.05, "lines")) +
  geom_line(size = 0.2, color = "black") + 
  geom_line(size = 1) +
  facet_wrap(~Period, axes = "all_y") +
  scale_color_gradient2(low = "#2166AC", mid = "transparent", high = "#B2182B", 
                        midpoint = 0, na.value = "transparent")  +    
  theme_classic() + 
  labs(x = "Month", 
       y = "Sea surface temperatue")+
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand= c(0,5)) +
  coord_cartesian(clip = 'off') +
  theme(panel.spacing = unit(4, "lines"), 
        plot.margin = unit(c(1,2,1,0.5), "lines"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p_shapley_lines
ggsave(filename = "plots/shapley_lines.pdf", plot = p_shapley_lines, height = 3, width = 7)
