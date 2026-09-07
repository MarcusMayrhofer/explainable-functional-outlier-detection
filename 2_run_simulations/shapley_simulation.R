library(robustmatrix)
library(foreach)
library(doParallel)
library(tidyverse)
library(doSNOW)
library(ggh4x)
library(pROC)

library(fda)
library(mrfDepth)
library(fdaoutlier)
library(abind)
library(rSPDE)
require(this.path)
setwd(this.path::this.dir())
source("../functions/multivariate_functional_pca_function.R")
source("../functions/multivariate_functional_data_functions.R")
source("../functions/helper_functions_simulation.R")
source("../functions/functions_for_simulations.R")
source("../functions/shapley_simulation_helper.R")


#MEAN FUNCTIONS
f1 <- function(t,lambda = 0.5) 30*t^(1 + 0.5*(1-lambda))*(1-t)^(1 + 0.5*(lambda))
f2 <- function(t, lambda = 1) 4*t + lambda*(-1)^rbinom(1, 1, 0.5)*(1.8-(0.02*pi)^(-0.5)*exp(-(t-runif(n = 1, min = 0.25, max = 0.75))^2/0.02))
f3 <- function(t, lambda = 1) 4*t + lambda*2*sin(8*(t+runif(n = 1, min = 0.25, max = 0.75))*pi)
mu1.0 = function(t, p) (replicate(n = p, expr = f1(t, lambda = 1)))
mu1.2 = function(t, p) (replicate(n = p, expr = f1(t, lambda = 0)))
mu2.0 = function(t, p) (replicate(n = p, expr = 4*t)) 
mu2.1 = function(t, p, lambda = 1) (replicate(n = p, expr = f2(t, lambda = lambda)))
mu3.1 = function(t, p) (replicate(n = p, expr = f3(t)))
mu4.1 = function(t, p) (replicate(n = p, expr = f3(t, lambda = 0.1)))

#COVARIANCE FUNCTIONS
K1.0 = function(s,t, rho = 0.3)  rho*exp((-abs(s-t))/rho) 
K2.0 = function(s,t, kappa = 5, nu = 0.5, sigma = 1) rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = sigma)
K2.1 = function(s,t, kappa = 10, nu = 0.2, sigma = 1) rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = sigma)



#########################################
from = 0 
to = 1
#########################################

n_run <- 100
eps_vec = c(0.1,0.3)
# p_vec = c(3,10,50) # other simulations not shown in paper
p_vec = c(10)
dim_vec <- c(1000)
outlier_shift_vec <- c(1,1.5,2)
cov_function_type = c("M", "OU")
n_basis_vec <- c(10,20,30)


npar <- n_run*length(p_vec)*length(eps_vec)*length(dim_vec)*length(outlier_shift_vec)*
  length(cov_function_type)*length(n_basis_vec)
cl <- makeCluster(120)
registerDoSNOW(cl)

# Progress Bar
pb <- txtProgressBar(max = npar, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tstart <- Sys.time()
res <-
  foreach(i1 = 1:n_run, .combine = "rbind") %:%
  foreach(i2 = 1:length(eps_vec), .combine = "rbind") %:%
  foreach(i4 = 1:length(dim_vec), .combine = "rbind") %:%
  foreach(i5 = 1:length(outlier_shift_vec), .combine = "rbind") %:%
  foreach(i7 = 1:length(cov_function_type), .combine = "rbind") %:%
  foreach(i8 = 1:length(n_basis_vec), .combine = "rbind") %:%
  foreach(i3 = 1:length(p_vec), .combine = "rbind", .packages = c("robustmatrix", "mrfDepth", "fdaoutlier", "fda", "dplyr", "tidyr", "tibble", "rSPDE"), .options.snow = opts) %dopar% {
    run <- i1
    set.seed(run)
    outlier_percentage = eps_vec[i2]
    p = p_vec[i3]
    n = dim_vec[i4]
    outlier_shift = outlier_shift_vec[i5] # controling strenght of shape and shift outliers
    cov_function = cov_function_type[i7]
    n_basis <- n_basis_vec[i8]

    q = 100
    quant <- qchisq(0.99, p*q)
    
    random_mu = FALSE
    mu_function <- mu1.0
    mu_function_outlier <- mu1.0
    
    Sigma <- cellWise::generateCorMat(d = p, corrType = "ALYZ")
    Sigma_outlier <- Sigma

    if(cov_function == "OU"){
      K <- K1.0
      K_outlier <- K1.0
    } else{
      K <- K2.0
      K_outlier <- K2.0
    }
    
    set.seed(run) # we set the seed again to be able to create the same, but shifted outliers if we only shift in some of the coordinates
    create_data <- msp(from = from, to = to, p = p, q = q, n = n, 
                       mu = mu_function, cov_row = Sigma, K = K)
    mu = create_data$mu
    cov_row = create_data$cov_row
    cov_col = create_data$cov_col
    eigen_cov_col <- eigen(cov_col)
    cov_col_projected = list(
      "proj_k5" = tcrossprod(eigen_cov_col$vectors[,1:5]),
      "proj_k10" = tcrossprod(eigen_cov_col$vectors[,1:10]),
      "proj_k20" = tcrossprod(eigen_cov_col$vectors[,1:20])
    )
    cov_col_projected$cov_k5 = cov_col_projected$proj_k5%*%cov_col%*%cov_col_projected$proj_k5
    cov_col_projected$cov_k10 = cov_col_projected$proj_k10%*%cov_col%*%cov_col_projected$proj_k10
    cov_col_projected$cov_k20 = cov_col_projected$proj_k20%*%cov_col%*%cov_col_projected$proj_k20
    cov_row_inv = solve(create_data$cov_row)
    cov_col_inv = solve(create_data$cov_col)
    n_outlier <- ceiling(n*outlier_percentage)
    X <- create_data$X
    X_clean <- X
    
    outlier_index <- 0
    if(outlier_percentage > 0){
      outlier_index <- sample(1:n, n_outlier)
      for(i in outlier_index){
        n_coord <- ceiling(runif(1,0,p-1))
        outlier_coord <- sample(1:p, n_coord)
        for(j in outlier_coord){
          sub_ind <- 1:(q)
          weight_mean <- sample(sub_ind, 1)
          weight_sd <- runif(1, min = q/25, q/5)
          prob <- dnorm(x = sub_ind, mean = weight_mean, sd = weight_sd)
          weight <- prob/max(prob)
          weight[weight < 0.5] <- 0
          outlying_time <- which(weight > 0)
          wave <- sin(outlying_time/(length(outlying_time)/ceiling(runif(1,1,15))))
          X[j,outlying_time,i] <- X[j,outlying_time,i] + runif(1,outlier_shift,outlier_shift+1)*weight[outlying_time]*wave
        }
      }
    }
    
    quant_fda <- qchisq(0.99, p*n_basis)
    basis <- create.bspline.basis(rangeval = c(from,to), nbasis = n_basis, norder = 4)
    X_smooth_basis <- smooth.basis(argvals = as.numeric(dimnames(X)[[2]]), y = aperm(X, perm = c(2,1,3)), fdParobj = basis)
    coeffs <- coef(X_smooth_basis)
    X_smoothed <- X_smooth_basis$fd
    X_smooth <- aperm(eval.fd(evalarg = create_data$eval_arg, fdobj = X_smoothed), c(2,1,3))
    
    info_raw = data.frame("n" = n, 
                          "eps" = outlier_percentage, 
                          "run" = run, 
                          "p" = p, 
                          "q" = q, 
                          "n_basis" = "raw", 
                          "cov_matrix" = "ALYZ", 
                          "cov_function" = cov_function, 
                          "outlier_shift" = outlier_shift)
    res_projdepth <-  outlier_summary_skeleton(method = "projdepth", info = info_raw)
    res_sprojdepth <-  outlier_summary_skeleton(method = "sprojdepth", info = info_raw)
    res_dprojdepth <-  outlier_summary_skeleton(method = "dprojdepth", info = info_raw)
    res_ms_plot <-  outlier_summary_skeleton(method = "ms_plot", info = info_raw)
    
    info_smooth = data.frame("n" = n, 
                             "eps" = outlier_percentage, 
                             "run" = run, 
                             "p" = p, 
                             "q" = q, 
                             "n_basis" = n_basis, 
                             "cov_matrix" = "ALYZ", 
                             "cov_function" = cov_function, 
                             "outlier_shift" = outlier_shift)
    res_mmcd_fda <-  outlier_summary_skeleton(method = "mmcd_fda", info = info_smooth)
      try({
        tstart_projdepth <- Sys.time()
        projdepth <- mrfDepth::mfd(x = aperm(X,c(2,3,1)), type = "projdepth", diagnostic = TRUE)
        t_projdepth <- Sys.time()-tstart_projdepth
        outliers_projdepth <- which(rowSums(projdepth$locOutlX)!=0)
        par_projdepth <- mmle(X[,,if(length(outliers_projdepth) == 0){TRUE} else{-outliers_projdepth}])
        res_projdepth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_projdepth,
                                         method = "projdepth", info = info_raw, time = t_projdepth, 
                                         mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                         outlier_score = as.numeric(projdepth$MFDdepthX),
                                         mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                         mu_function_est = par_projdepth$mu, 
                                         cov_matrix_est = par_projdepth$cov_row, 
                                         cov_function_est = par_projdepth$cov_col,
                                         eval_arg = NULL)
      })
      try({
        tstart_sprojdepth <- Sys.time()
        sprojdepth <- mrfDepth::mfd(x = aperm(X,c(2,3,1)), type = "sprojdepth", diagnostic = TRUE)
        t_sprojdepth <- Sys.time()-tstart_sprojdepth
        outliers_sprojdepth <- which(rowSums(sprojdepth$locOutlX)!=0)
        par_sprojdepth <- mmle(X[,,if(length(outliers_sprojdepth) == 0){TRUE} else{-outliers_sprojdepth}])
        res_sprojdepth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_sprojdepth,
                                          method = "sprojdepth", info = info_raw, time = t_sprojdepth, 
                                          mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                          outlier_score = as.numeric(sprojdepth$MFDdepthX),
                                          mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                          mu_function_est = par_sprojdepth$mu, 
                                          cov_matrix_est = par_sprojdepth$cov_row, 
                                          cov_function_est = par_sprojdepth$cov_col,
                                          eval_arg = NULL)
      })
      try({
        tstart_dprojdepth <- Sys.time()
        dprojdepth <- mrfDepth::mfd(x = aperm(X,c(2,3,1)), type = "dprojdepth", diagnostic = TRUE)
        t_dprojdepth <- Sys.time()-tstart_dprojdepth
        outliers_dprojdepth <- which(rowSums(dprojdepth$locOutlX)!=0)
        par_dprojdepth <- mmle(X[,,if(length(outliers_dprojdepth) == 0){TRUE} else{-outliers_dprojdepth}])
        res_dprojdepth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_dprojdepth,
                                          method = "dprojdepth", info = info_raw, time = t_dprojdepth, 
                                          mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                          outlier_score = as.numeric(dprojdepth$MFDdepthX),
                                          mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                          mu_function_est = par_dprojdepth$mu, 
                                          cov_matrix_est = par_dprojdepth$cov_row, 
                                          cov_function_est = par_dprojdepth$cov_col,
                                          eval_arg = NULL)
      })
      try({
        tstart_ms_plot <- Sys.time()
        ms_plot <- msplot_adv(aperm(X,c(3,2,1)), plot = FALSE)
        t_ms_plot <- Sys.time()-tstart_ms_plot
        outliers_ms_plot <- ms_plot$outliers
        par_ms_plot <- mmle(X[,,if(length(outliers_ms_plot) == 0){TRUE} else{-outliers_ms_plot}])
        res_ms_plot <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_ms_plot,
                                       method = "ms_plot", info = info_raw, time = t_ms_plot,
                                       mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                       outlier_score = ms_plot$dir_out_distance,
                                       mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                       mu_function_est = par_ms_plot$mu,
                                       cov_matrix_est = par_ms_plot$cov_row,
                                       cov_function_est = par_ms_plot$cov_col,
                                       eval_arg = NULL)
      })
    ###################################################################################
    ###################################################################################
    #smooth and coef methods
    ###################################################################################
    ###################################################################################
    try({
      tstart_mmcd_fda <- Sys.time()
      par_mmcd_fda <- mmcd(coeffs, scale_consistency = "mmd_med", outlier_quant = 0.99)
      t_mmcd_fda <- Sys.time()-tstart_mmcd_fda
      pca_X_mmcd <- pca.mfd(fdobj = X_smoothed, par = par_mmcd_fda)
      mmd_mmcd_fda <- mmd(coeffs, par_mmcd_fda$mu, par_mmcd_fda$cov_row_inv, par_mmcd_fda$cov_col_inv, inverted = TRUE)
      outliers_mmcd_fda <- which(mmd_mmcd_fda > quant_fda)
      
      res_mmcd_fda <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_mmcd_fda,
                                      method = "mmcd_fda", info = info_smooth, time = t_mmcd_fda, 
                                      mutlivariate_score = TRUE, functional = TRUE, reference_method = res_mmle_clean,
                                      outlier_score = mmd_mmcd_fda,
                                      mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                      mu_function_est = pca_X_mmcd$mu_function, 
                                      cov_matrix_est = pca_X_mmcd$par$cov_col, 
                                      cov_function_est = pca_X_mmcd$cov_function,
                                      eval_arg = create_data$eval_arg)
    })
    
    
    AUC_all <- data.frame(out_numb = NA,
                          method = "dummy",
                          type = NA,
                          absolute_score = NA,
                          n_interval = NA,
                          AUC = NA)
    
    try({
      outliers_detected <- outlier_index
      X_diff <- abs(X - X_clean)
      
      
      shap_mmcd_fda <- shapley.mfd(fdobj = X_smoothed, par = par_mmcd_fda, time_intervals = 100)
      shap_actual <-  matrixShapley(X, mu, cov_row_inv, cov_col_inv, inverted = TRUE)
      abs_centered_data <- abs(sweep(X, 1:2, mu, "-"))
      
      AUC_actual <- data.frame()
      AUC_mmcd_fda <- data.frame()
      AUC_abs_centered <- data.frame()
      AUC_projdepth <- data.frame()
      AUC_sprojdepth <- data.frame()
      AUC_dprojdepth <- data.frame()
      AUC_ms <- data.frame()
      n_interval_vec <- c(1,5,10,20)
      for(i in seq_along(outliers_detected)){
        sel_id <- outliers_detected[i]
        AUC_mmcd_fda <- rbind(
          AUC_mmcd_fda,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = shap_mmcd_fda$shv[,,sel_id], X_diff = X_diff[,,sel_id], n_interval = n_interval_vec,
                                                  method = "mmcd_fda", types = c("time_coord", "time", "coord")))
        )
        AUC_actual <- rbind(
          AUC_actual,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = shap_actual[,,sel_id], X_diff = X_diff[,,sel_id], n_interval = n_interval_vec,
                                                  method = "actual", types = c("time_coord", "time", "coord")))
        )
        AUC_abs_centered <- rbind(
          AUC_abs_centered,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = abs_centered_data[,,sel_id], X_diff = X_diff[,,sel_id], n_interval = n_interval_vec,
                                                  method = "abs_centered", types = c("time_coord", "time", "coord")))
        )
        AUC_projdepth <- rbind(
          AUC_projdepth,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = t(matrix(projdepth$crossdepthX[sel_id,])), X_diff = colSums(X_diff[,,sel_id]),
                                                  n_interval = n_interval_vec, method = "projdepth", types = c("time")))
        )
        AUC_sprojdepth <- rbind(
          AUC_sprojdepth,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = t(matrix(sprojdepth$crossdepthX[sel_id,])), X_diff = colSums(X_diff[,,sel_id]),
                                                  n_interval = n_interval_vec, method = "sprojdepth", types = c("time")))
        )
        AUC_dprojdepth <- rbind(
          AUC_dprojdepth,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = t(matrix(dprojdepth$crossdepthX[sel_id,])), X_diff = colSums(X_diff[,,sel_id]),
                                                  n_interval = n_interval_vec, method = "dprojdepth", types = c("time")))
        )
        AUC_ms <- rbind(
          AUC_ms,
          cbind("out_numb" = i, AUC_func_advanced(outlier_score = t(matrix(ms_plot$mean_outlyingness[sel_id,])), X_diff = rowSums(X_diff[,,sel_id]), 
                                                  n_interval = n_interval_vec, method = "ms", types = c("coord")))
        )
      }
      
      AUC_all <- rbind(AUC_mmcd_fda,
                       AUC_actual,
                       AUC_abs_centered,
                       AUC_projdepth,
                       AUC_sprojdepth,
                       AUC_dprojdepth,
                       AUC_ms)
    })
    cbind(info_smooth,AUC_all)
  }
tend <- Sys.time()
tend-tstart
close(pb)
stopCluster(cl)
save("res",file = "raw_results/shapley_simulation_res.RData")



