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
eps_vec = c(0.05,0.1,0.2,0.3)
p_vec = c(3)
dim_vec <- c(300,1000)
#shape and shift outliers
  # eigen_shift_vec <- NA
  # eigen_shift_vec <- c(16,20,24) #use this for eigen_id = 1, p = 3
  # eigen_shift_vec <- c(12,15,18) #use this for eigen_id = 1, p = 10
  # eigen_shift_vec <- c(8,10,12) #use this for eigen_id = 1, p = 50
  eigen_shift_vec <- c(5,8,11) #use this for eigen_id = 10, p = 3
  # eigen_shift_vec <- c(2,3,4) #use this for eigen_id = 10, p = 10
  # eigen_shift_vec <- c(1.5,2,2.5) #use this for eigen_id = 10, p = 50
  # eigen_id_vec <- NA
  # eigen_id_vec <- c(1,10)
  eigen_id_vec <- c(10)
cov_function_type = c("M", "OU")
n_basis_vec <- c(10,20,30)
eps_coord_vec <- c(0.1,0.5,1)
#isolated outleirs
  lambda_vec <- NA
  # lambda_vec <- c(2,5,10)/10
#covariance outliers
  kappa_vec <- NA
  # kappa_vec <- c(7, 10, 15)
  nu_vec <- NA
  # nu_vec <- c(0.1,0.2,0.5)


npar <- n_run*length(p_vec)*length(eps_vec)*length(dim_vec)*length(eigen_shift_vec)*
  length(cov_function_type)*length(eigen_id_vec)*length(n_basis_vec)*
  length(eps_coord_vec)*length(lambda_vec)*length(kappa_vec)*length(nu_vec)

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
  foreach(i5 = 1:length(eigen_shift_vec), .combine = "rbind") %:%
  foreach(i6 = 1:length(eigen_id_vec), .combine = "rbind") %:%
  foreach(i7 = 1:length(cov_function_type), .combine = "rbind") %:%
  foreach(i8 = 1:length(n_basis_vec), .combine = "rbind") %:%
  foreach(i9 = 1:length(eps_coord_vec), .combine = "rbind") %:%
  foreach(i10 = 1:length(lambda_vec), .combine = "rbind") %:%
  foreach(i11 = 1:length(kappa_vec), .combine = "rbind") %:%
  foreach(i12 = 1:length(nu_vec), .combine = "rbind") %:%
  foreach(i3 = 1:length(p_vec), .combine = "rbind", .packages = c("robustmatrix", "mrfDepth", "fdaoutlier", "fda", "dplyr", "tidyr", "tibble", "rSPDE"), .options.snow = opts) %dopar% {
    run <- i1
    set.seed(run)
    outlier_percentage = eps_vec[i2]
    p = p_vec[i3]
    n = dim_vec[i4]
    eigen_shift = eigen_shift_vec[i5] # controling strenght of shape and shift outliers
    eigen_id = eigen_id_vec[i6] # controlling shift or shape outliers
    cov_function = cov_function_type[i7]
    n_basis <- n_basis_vec[i8]
    outlier_coord_percentage = eps_coord_vec[i9] # controls in which coordinates the outliers are
    lambda = lambda_vec[i10] # controling the isolated outliers
    kappa = kappa_vec[i11] # controlls covariance induced outliers
    nu = nu_vec[i12] # controlls covariance induced outliers
    
    q = 100
    # q = 50
    quant <- qchisq(0.99, p*q)
    
    random_mu = FALSE
    if(!is.na(lambda_vec)){
      random_mu = TRUE
      mu_function <- mu2.0
      mu_function_outlier <- function(t, p ,lambda_sel = lambda) mu2.1(t, p, lambda_sel)
    } else{
      mu_function <- mu1.0
      mu_function_outlier <- mu1.0
    }
    
    Sigma <- cellWise::generateCorMat(d = p, corrType = "ALYZ")
    Sigma_outlier <- Sigma
    
    if(!is.na(kappa_vec) | !is.na(nu_vec)){
      K <- K2.0
      K_outlier <- function(s, t, kappa_sel = kappa, nu_sel = nu, sigma = 1) K2.1(s, t, kappa_sel, nu_sel, sigma)
    } else{
      if(cov_function == "OU"){
        K <- K1.0
        K_outlier <- K1.0
      } else{
        K <- K2.0
        K_outlier <- K2.0
      }
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
      #here we set the seed again, we create the same process as for the clean data but this time shifted by gamma or with a different covariance/lambda
      set.seed(run)
      if(!is.na(eigen_id) | !is.na(eigen_shift)){
        create_data_outlier <- msp(from = from, to = to, p = p, q = q, n = n,
                                   mu = mu_function_outlier, cov_row = Sigma_outlier, 
                                   K = K_outlier, mu_rand = TRUE, eigen_id = eigen_id, gamma = eigen_shift)
      } else{
        create_data_outlier <- msp(from = from, to = to, p = p, q = q, n = n,
                                   mu = mu_function_outlier, cov_row = Sigma_outlier, 
                                   K = K_outlier, mu_rand = TRUE)
      }
      X_outlier <- create_data_outlier$X
      outlier_index <- sample(1:n, n_outlier)
      outlier_coord <- ceiling(outlier_coord_percentage*p)
      X[1:outlier_coord,,outlier_index] <- X_outlier[1:outlier_coord,,outlier_index]
    }
    
    quant_fda <- qchisq(0.99, p*n_basis)
    basis <- create.bspline.basis(rangeval = c(from,to), nbasis = n_basis, norder = 4)
    X_smooth_basis <- smooth.basis(argvals = as.numeric(dimnames(X)[[2]]), y = aperm(X, perm = c(2,1,3)), fdParobj = basis)
    coeffs <- coef(X_smooth_basis)
    X_smoothed <- X_smooth_basis$fd
    X_smooth <- aperm(eval.fd(evalarg = create_data$eval_arg, fdobj = X_smoothed), c(2,1,3))
    
    info_raw = data.frame("n" = n, 
                          "eps" = outlier_percentage, 
                          "eps_coord_vec" = outlier_coord_percentage, 
                          "run" = run, 
                          "p" = p, 
                          "q" = q, 
                          "n_basis" = "raw", 
                          "cov_matrix" = "ALYZ", 
                          "cov_function" = cov_function, 
                          "eigen_id" = eigen_id, 
                          "eigen_shift" = eigen_shift,
                          "nu" = nu, 
                          "kappa" = kappa,
                          "lambda" = lambda)
    res_actual <- outlier_summary_skeleton(method = "actual", info = info_raw)
    res_mmle_clean <-  outlier_summary_skeleton(method = "mmle_clean", info = info_raw)
    res_mmle <-  outlier_summary_skeleton(method = "mmle", info = info_raw)
    res_mmcd <-  outlier_summary_skeleton(method = "mmcd", info = info_raw)
    res_projdepth <-  outlier_summary_skeleton(method = "projdepth", info = info_raw)
    res_sprojdepth <-  outlier_summary_skeleton(method = "sprojdepth", info = info_raw)
    res_dprojdepth <-  outlier_summary_skeleton(method = "dprojdepth", info = info_raw)
    res_ms_plot <-  outlier_summary_skeleton(method = "ms_plot", info = info_raw)
    
    info_smooth = data.frame("n" = n, 
                             "eps" = outlier_percentage, 
                             "eps_coord_vec" = outlier_coord_percentage, 
                             "run" = run, 
                             "p" = p, 
                             "q" = q, 
                             "n_basis" = n_basis, 
                             "cov_matrix" = "ALYZ", 
                             "cov_function" = cov_function, 
                             "eigen_id" = eigen_id, 
                             "eigen_shift" = eigen_shift,
                             "nu" = nu, 
                             "kappa" = kappa,
                             "lambda" = lambda)
    res_mmle_fda <-  outlier_summary_skeleton(method = "mmle_fda", info = info_smooth)
    res_mmcd_fda <-  outlier_summary_skeleton(method = "mmcd_fda", info = info_smooth)
    res_projdepth_smooth <-  outlier_summary_skeleton(method = "projdepth_smooth", info = info_smooth)
    res_projdepth_coef <- outlier_summary_skeleton(method = "projdepth_coef", info = info_smooth)
    res_sprojdepth_smooth <-  outlier_summary_skeleton(method = "sprojdepth_smooth", info = info_smooth)
    res_sprojdepth_coef <- outlier_summary_skeleton(method = "sprojdepth_coef", info = info_smooth)
    res_dprojdepth_smooth <-  outlier_summary_skeleton(method = "dprojdepth_smooth", info = info_smooth)
    res_dprojdepth_coef <- outlier_summary_skeleton(method = "dprojdepth_coef", info = info_smooth)
    res_ms_plot_smooth <-  outlier_summary_skeleton(method = "ms_plot_smooth", info = info_smooth)
    res_ms_plot_coef <- outlier_summary_skeleton(method = "ms_plot_coef", info = info_smooth)
    try({
      tstart_mmle_clean <- Sys.time()
      par_mmle_clean <- mmle(X_clean)
      t_mmle_clean <- Sys.time()-tstart_mmle_clean
      mmd_mmle_clean <- mmd(X,par_mmle_clean$mu, par_mmle_clean$cov_row_inv, par_mmle_clean$cov_col_inv, inverted = TRUE)
      outliers_mmle_clean <- which(mmd_mmle_clean > quant)
      
      res_mmle_clean <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_mmle_clean,
                                        method = "mmle_clean", info = info_raw, time = t_mmle_clean, 
                                        mutlivariate_score = TRUE, functional = FALSE, reference_method = NULL,
                                        outlier_score = mmd_mmle_clean,
                                        mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                        mu_function_est = par_mmle_clean$mu, 
                                        cov_matrix_est = par_mmle_clean$cov_row, 
                                        cov_function_est = par_mmle_clean$cov_col,
                                        eval_arg = NULL)
    })
    if(n_basis == n_basis_vec[1]){
      try({
        tstart_mmle <- Sys.time()
        par_mmle <- mmle(X)
        t_mmle <- Sys.time()-tstart_mmle
        mmd_mmle <- mmd(X,par_mmle$mu, par_mmle$cov_row_inv, par_mmle$cov_col_inv, inverted = TRUE)
        outliers_mmle <- which(mmd_mmle > quant)
        
        res_mmle <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_mmle,
                                    method = "mmle", info = info_raw, time = t_mmle, 
                                    mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                    outlier_score = mmd_mmle,
                                    mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected, 
                                    mu_function_est = par_mmle$mu, 
                                    cov_matrix_est = par_mmle$cov_row, 
                                    cov_function_est = par_mmle$cov_col,
                                    eval_arg = NULL)
      })
      try({
        tstart_mmcd <- Sys.time()
        par_mmcd <- mmcd(X = X, scale_consistency = "mmd_med", outlier_quant = 0.99)
        t_mmcd <- Sys.time()-tstart_mmcd
        mmd_mmcd <- mmd(X,par_mmcd$mu, par_mmcd$cov_row_inv, par_mmcd$cov_col_inv, inverted = TRUE)
        outliers_mmcd <- which(mmd_mmcd > quant)
        
        res_mmcd <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_mmcd,
                                    method = "mmcd", info = info_raw, time = t_mmcd,
                                    mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                    outlier_score = mmd_mmcd,
                                    mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                    mu_function_est = par_mmcd$mu,
                                    cov_matrix_est = par_mmcd$cov_row,
                                    cov_function_est = par_mmcd$cov_col,
                                    eval_arg = NULL)
      })
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
      try({
        mmd_actual <- mmd(X, mu, cov_row_inv, cov_col_inv, inverted = TRUE)
        outliers_actual <- which(mmd_actual > quant)
        res_actual <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_actual,
                                      method = "actual", info = info_raw, time = 0, 
                                      mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                      outlier_score = mmd_actual,
                                      mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                      mu_function_est = mu, 
                                      cov_matrix_est = cov_row, 
                                      cov_function_est = cov_col,
                                      eval_arg = NULL)
      })
    }
    ###################################################################################
    ###################################################################################
    #smooth and coef methods
    ###################################################################################
    ###################################################################################
    try({
      tstart_mmle_fda <- Sys.time()
      par_mmle_fda <- mmle(coeffs)
      t_mmle_fda <- Sys.time()-tstart_mmle_fda
      pca_X_mmle <- pca.mfd(fdobj = X_smoothed, par = par_mmle_fda)
      mmd_mmle_fda <- mmd(coeffs, pca_X_mmle$par$mu, pca_X_mmle$par$cov_row_inv, pca_X_mmle$par$cov_col_inv, inverted = TRUE)
      outliers_mmle_fda <- which(mmd_mmle_fda > quant_fda)
      
      res_mmle_fda <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_mmle_fda,
                                      method = "mmle_fda", info = info_smooth, time = t_mmle_fda, 
                                      mutlivariate_score = TRUE, functional = TRUE, reference_method = res_mmle_clean,
                                      outlier_score = mmd_mmle_fda,
                                      mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                      mu_function_est = pca_X_mmle$mu_function, 
                                      cov_matrix_est = pca_X_mmle$par$cov_col, 
                                      cov_function_est = pca_X_mmle$cov_function,
                                      eval_arg = create_data$eval_arg)
    })
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
    try({
      tstart_projdepth_smooth <- Sys.time()
      projdepth_smooth <- mrfDepth::mfd(x = aperm(X_smooth,c(2,3,1)), type = "projdepth", diagnostic = TRUE)
      t_projdepth_smooth <- Sys.time()-tstart_projdepth_smooth
      outliers_projdepth_smooth <- which(rowSums(projdepth_smooth$locOutlX)!=0)
      par_projdepth_smooth <- mmle(X[,,if(length(outliers_projdepth_smooth) == 0){TRUE} else{-outliers_projdepth_smooth}])
      res_projdepth_smooth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_projdepth_smooth,
                                              method = "projdepth_smooth", info = info_smooth, time = t_projdepth_smooth, 
                                              mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                              outlier_score = as.numeric(projdepth_smooth$MFDdepthX),
                                              mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                              mu_function_est = par_projdepth_smooth$mu, 
                                              cov_matrix_est = par_projdepth_smooth$cov_row, 
                                              cov_function_est = par_projdepth_smooth$cov_col,
                                              eval_arg = NULL)
    })
    try({
      tstart_projdepth_coef <- Sys.time()
      projdepth_coef <- mrfDepth::mfd(x = aperm(coeffs,c(1,3,2)), type = "projdepth", diagnostic = TRUE)
      t_projdepth_coef <- Sys.time()-tstart_projdepth_coef
      outliers_projdepth_coef <- which(rowSums(projdepth_coef$locOutlX)!=0)
      par_projdepth_coef <- mmle(X[,,if(length(outliers_projdepth_coef) == 0){TRUE} else{-outliers_projdepth_coef}])
      res_projdepth_coef <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_projdepth_coef,
                                            method = "projdepth_coef", info = info_smooth, time = t_projdepth_coef, 
                                            mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                            outlier_score = as.numeric(projdepth_coef$MFDdepthX),
                                            mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                            mu_function_est = par_projdepth_coef$mu, 
                                            cov_matrix_est = par_projdepth_coef$cov_row, 
                                            cov_function_est = par_projdepth_coef$cov_col,
                                            eval_arg = NULL)
    })
    try({
      tstart_sprojdepth_smooth <- Sys.time()
      sprojdepth_smooth <- mrfDepth::mfd(x = aperm(X_smooth,c(2,3,1)), type = "sprojdepth", diagnostic = TRUE)
      t_sprojdepth_smooth <- Sys.time()-tstart_sprojdepth_smooth
      outliers_sprojdepth_smooth <- which(rowSums(sprojdepth_smooth$locOutlX)!=0)
      par_sprojdepth_smooth <- mmle(X[,,if(length(outliers_sprojdepth_smooth) == 0){TRUE} else{-outliers_sprojdepth_smooth}])
      res_sprojdepth_smooth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_sprojdepth_smooth,
                                               method = "sprojdepth_smooth", info = info_smooth, time = t_sprojdepth_smooth, 
                                               mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                               outlier_score = as.numeric(sprojdepth_smooth$MFDdepthX),
                                               mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                               mu_function_est = par_sprojdepth_smooth$mu, 
                                               cov_matrix_est = par_sprojdepth_smooth$cov_row, 
                                               cov_function_est = par_sprojdepth_smooth$cov_col,
                                               eval_arg = NULL)
    })
    try({
      tstart_sprojdepth_coef <- Sys.time()
      sprojdepth_coef <- mrfDepth::mfd(x = aperm(coeffs,c(1,3,2)), type = "sprojdepth", diagnostic = TRUE)
      t_sprojdepth_coef <- Sys.time()-tstart_sprojdepth_coef
      outliers_sprojdepth_coef <- which(rowSums(sprojdepth_coef$locOutlX)!=0)
      par_sprojdepth_coef <- mmle(X[,,if(length(outliers_sprojdepth_coef) == 0){TRUE} else{-outliers_sprojdepth_coef}])
      res_sprojdepth_coef <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_sprojdepth_coef,
                                             method = "sprojdepth_coef", info = info_smooth, time = t_sprojdepth_coef, 
                                             mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                             outlier_score = as.numeric(sprojdepth_coef$MFDdepthX),
                                             mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                             mu_function_est = par_sprojdepth_coef$mu, 
                                             cov_matrix_est = par_sprojdepth_coef$cov_row, 
                                             cov_function_est = par_sprojdepth_coef$cov_col,
                                             eval_arg = NULL)
    })
    try({
      tstart_dprojdepth_smooth <- Sys.time()
      dprojdepth_smooth <- mrfDepth::mfd(x = aperm(X_smooth,c(2,3,1)), type = "dprojdepth", diagnostic = TRUE)
      t_dprojdepth_smooth <- Sys.time()-tstart_dprojdepth_smooth
      outliers_dprojdepth_smooth <- which(rowSums(dprojdepth_smooth$locOutlX)!=0)
      par_dprojdepth_smooth <- mmle(X[,,if(length(outliers_dprojdepth_smooth) == 0){TRUE} else{-outliers_dprojdepth_smooth}])
      res_dprojdepth_smooth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_dprojdepth_smooth,
                                               method = "dprojdepth_smooth", info = info_smooth, time = t_dprojdepth_smooth, 
                                               mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                               outlier_score = as.numeric(dprojdepth_smooth$MFDdepthX),
                                               mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                               mu_function_est = par_dprojdepth_smooth$mu, 
                                               cov_matrix_est = par_dprojdepth_smooth$cov_row, 
                                               cov_function_est = par_dprojdepth_smooth$cov_col,
                                               eval_arg = NULL)
    })
    try({
      tstart_dprojdepth_coef <- Sys.time()
      dprojdepth_coef <- mrfDepth::mfd(x = aperm(coeffs,c(1,3,2)), type = "dprojdepth", diagnostic = TRUE)
      t_dprojdepth_coef <- Sys.time()-tstart_dprojdepth_coef
      outliers_dprojdepth_coef <- which(rowSums(dprojdepth_coef$locOutlX)!=0)
      par_dprojdepth_coef <- mmle(X[,,if(length(outliers_dprojdepth_coef) == 0){TRUE} else{-outliers_dprojdepth_coef}])
      res_dprojdepth_coef <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_dprojdepth_coef,
                                             method = "dprojdepth_coef", info = info_smooth, time = t_dprojdepth_coef, 
                                             mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                             outlier_score = as.numeric(dprojdepth_coef$MFDdepthX),
                                             mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                             mu_function_est = par_dprojdepth_coef$mu, 
                                             cov_matrix_est = par_dprojdepth_coef$cov_row, 
                                             cov_function_est = par_dprojdepth_coef$cov_col,
                                             eval_arg = NULL)
    })
    try({
      tstart_ms_plot_smooth <- Sys.time()
      ms_plot_smooth <- msplot_adv(aperm(X_smooth,c(3,2,1)), plot = FALSE)
      t_ms_plot_smooth <- Sys.time()-tstart_ms_plot_smooth
      outliers_ms_plot_smooth <- ms_plot_smooth$outliers
      par_ms_plot_smooth <- mmle(X[,,if(length(outliers_ms_plot_smooth) == 0){TRUE} else{-outliers_ms_plot_smooth}])
      res_ms_plot_smooth <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_ms_plot_smooth,
                                            method = "ms_plot_smooth", info = info_smooth, time = t_ms_plot_smooth, 
                                            mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                            outlier_score = ms_plot_smooth$dir_out_distance,
                                            mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                            mu_function_est = par_ms_plot_smooth$mu, 
                                            cov_matrix_est = par_ms_plot_smooth$cov_row, 
                                            cov_function_est = par_ms_plot_smooth$cov_col,
                                            eval_arg = NULL)
    })
    try({
      tstart_ms_plot_coef <- Sys.time()
      ms_plot_coef <- msplot_adv(aperm(coeffs,c(3,1,2)), plot = FALSE)
      t_ms_plot_coef <- Sys.time()-tstart_ms_plot_coef
      outliers_ms_plot_coef <- ms_plot_coef$outliers
      par_ms_plot_coef <- mmle(X[,,if(length(outliers_ms_plot_coef) == 0){TRUE} else{-outliers_ms_plot_coef}])
      res_ms_plot_coef <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_ms_plot_coef,
                                          method = "ms_plot_coef", info = info_smooth, time = t_ms_plot_coef, 
                                          mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                          outlier_score = ms_plot_coef$dir_out_distance,
                                          mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                          mu_function_est = par_ms_plot_coef$mu, 
                                          cov_matrix_est = par_ms_plot_coef$cov_row, 
                                          cov_function_est = par_ms_plot_coef$cov_col,
                                          eval_arg = NULL)
    })
    rbind(res_mmle, res_mmcd, res_mmle_fda, res_mmcd_fda, 
          res_projdepth, res_projdepth_smooth, res_projdepth_coef,
          res_sprojdepth, res_sprojdepth_smooth, res_sprojdepth_coef,
          res_dprojdepth, res_dprojdepth_smooth, res_dprojdepth_coef,
          res_ms_plot, res_ms_plot_smooth, res_ms_plot_coef,
          res_mmle_clean, res_actual)
  }
tend <- Sys.time()
tend-tstart
close(pb)
stopCluster(cl)
save("res",file = "raw_results/simulation_res_p3_shape.RData")
