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
source("../functions/functions_for_simulations_non_separable.R")
source("../functions/rank_1_approx.R")



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
K2.1 = function(s,t, kappa = 10, nu = 2.5, sigma = 1.2) rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = sigma)
K2.2 = function(s,t, kappa = 20, nu = 2, sigma = 1) (rSPDE::matern.covariance(h = abs(s-t), kappa = kappa*10, nu = nu/4, sigma = sigma*2) + rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = sigma))
K2.3 = function(s,t, kappa = 20, nu = 0.5, sigma = 1) {
  sigma_s <- 0.3 + 0.8 * (0.5 + 0.5 * sin(2 * pi * s))
  sigma_t <- 0.3 + 0.8 * (0.5 + 0.5 * sin(2 * pi * t))
  sigma_s * sigma_t * rSPDE::matern.covariance(h = abs(s-t), kappa = kappa, nu = nu, sigma = 1)
}

#########################################
from = 0 
to = 1
#########################################

n_run <- 100
eps_vec = c(0.05,0.1,0.2,0.3)
# p_vec = c(3,10,50)
p_vec = c(10)
dim_vec <- c(1500,5000)
#shape and shift outliers
# eigen_shift_vec <- NA
# eigen_shift_vec <- c(15,22.5,30) #use this for eigen_id = 1, p = 3
eigen_shift_vec <- c(15,30)#c(12,15,18) #use this for eigen_id = 1, p = 10
# eigen_shift_vec <- c(8,10,12) #use this for eigen_id = 1, p = 50
# eigen_shift_vec <- c(5,8,11) #use this for eigen_id = 10, p = 3
# eigen_shift_vec <- c(2,3,4) #use this for eigen_id = 10, p = 10
# eigen_shift_vec <- c(1.5,2,2.5) #use this for eigen_id = 10, p = 50
# eigen_id_vec <- NA
# eigen_id_vec <- c(1,10)
eigen_id_vec <- c(1)
cov_function_type = "M"
n_basis_vec <- c(10,20,30)
eps_coord_vec <- 1
#isolated outleirs
lambda_vec <- NA
# lambda_vec <- c(2,5,10)/10
#covariance outliers
kappa_vec <- NA
# kappa_vec <- c(7, 10, 15)
nu_vec <- NA
# nu_vec <- c(0.1,0.2,0.5)
kron_sum_vec <- c(1:4)


npar <- n_run*length(p_vec)*length(eps_vec)*length(dim_vec)*length(eigen_shift_vec)*
  length(cov_function_type)*length(eigen_id_vec)*length(n_basis_vec)*
  length(eps_coord_vec)*length(lambda_vec)*length(kappa_vec)*length(nu_vec)*length(kron_sum_vec)
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
  foreach(i13 = 1:length(kron_sum_vec), .combine = "rbind") %:%
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
    if(!any(is.na(lambda_vec))){
      random_mu = TRUE
      mu_function <- mu2.0
      mu_function_outlier <- function(t, p ,lambda_sel = lambda) mu2.1(t, p, lambda_sel)
    } else{
      mu_function <- mu1.0
      mu_function_outlier <- mu1.0
    }
    
    Sigma <- cellWise::generateCorMat(d = p, corrType = "ALYZ")
    Sigma_outlier <- Sigma
    
    if(!any(is.na(kappa_vec)) | !any(is.na(nu_vec))){
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
    Sigma1 <- cellWise::generateCorMat(d = p, corrType = "ALYZ")
    Sigma2 <- cellWise::generateCorMat(d = p, corrType = "A09")
    Sigma3 <- cellWise::generateCorMat(d = p, corrType = "ALYZ")
    Sigma4 <- matrix(0.5, nrow = p, ncol = p)
    diag(Sigma4) <- 1
    Sigma_outlier <- Sigma1

    K1 <- K2.0
    K2 <- K2.1
    K3 <- K2.2
    K4 <- K2.3
    
    set.seed(run) # we set the seed again to be able to create the same, but shifted outliers if we only shift in some of the coordinates
    create_data1 <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function, cov_row = Sigma1, K = K1)
    create_data2 <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function, cov_row = Sigma2, K = K2)
    create_data3 <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function, cov_row = Sigma3, K = K3)
    create_data4 <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function, cov_row = Sigma4, K = K4)
    
    n_outlier <- ceiling(n*outlier_percentage)

    X1 <- X1_clean <- create_data1$X
    X2 <- X2_clean <- create_data2$X
    X3 <- X3_clean <- create_data3$X
    X4 <- X4_clean <- create_data4$X
    
    mu1 <- create_data1$mu
    mu2 <- create_data2$mu
    mu3 <- create_data3$mu
    mu4 <- create_data4$mu
    
    outlier_index <- 0
    if(outlier_percentage > 0){
      #here we set the seed again, we create the same process as for the clean data but this time shifted by gamma or with a different covariance/lambda
      set.seed(run)
      create_data1_outlier <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function_outlier, cov_row = Sigma1, K = K1, 
                                  mu_rand = TRUE, eigen_id = eigen_id, gamma = eigen_shift)
      create_data2_outlier <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function_outlier, cov_row = Sigma2, K = K2, 
                                  mu_rand = TRUE, eigen_id = eigen_id + 1, gamma = eigen_shift)
      create_data3_outlier <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function_outlier, cov_row = Sigma3, K = K3, 
                                  mu_rand = TRUE, eigen_id = eigen_id + 2, gamma = eigen_shift)
      create_data4_outlier <- msp(from = from, to = to, p = p, q = q, n = n, mu = mu_function_outlier, cov_row = Sigma4, K = K4, 
                                  mu_rand = TRUE, eigen_id = eigen_id + 3, gamma = eigen_shift)
      X1_outlier <- create_data1_outlier$X
      X2_outlier <- create_data2_outlier$X
      X3_outlier <- create_data3_outlier$X
      X4_outlier <- create_data4_outlier$X
      outlier_index <- sample(1:n, n_outlier)
      outlier_coord <- ceiling(outlier_coord_percentage*p)
      X1[1:outlier_coord,,outlier_index] <- X1_outlier[1:outlier_coord,,outlier_index]
      X2[1:outlier_coord,,outlier_index] <- X2_outlier[1:outlier_coord,,outlier_index]
      X3[1:outlier_coord,,outlier_index] <- X3_outlier[1:outlier_coord,,outlier_index]
      X4[1:outlier_coord,,outlier_index] <- X4_outlier[1:outlier_coord,,outlier_index]
    }
    
    if(i13 == 1){
      X_clean <- X1_clean
      X <- X1
      kron_true <- kronecker(create_data1$cov_col, create_data1$cov_row)
      mu <- mu1
    } else if(i13 == 2){
      X_clean <- X1_clean + X2_clean
      X <- X1 + X2
      kron_true1 <- kronecker(create_data1$cov_col, create_data1$cov_row)
      kron_true2 <- kronecker(create_data2$cov_col, create_data2$cov_row)
      kron_true <- kron_true1 + kron_true2
      mu <- mu1 + mu2
    } else if(i13 == 3){
      X_clean <- X1_clean + X2_clean + X3_clean
      X <- X1 + X2 + X3
      kron_true1 <- kronecker(create_data1$cov_col, create_data1$cov_row)
      kron_true2 <- kronecker(create_data2$cov_col, create_data2$cov_row)
      kron_true3 <- kronecker(create_data3$cov_col, create_data3$cov_row)
      kron_true <- kron_true1 + kron_true2 + kron_true3
      mu <- mu1 + mu2 + mu3
    } else if(i13 == 4){
      X_clean <- X1_clean + X2_clean + X3_clean + X4_clean
      X <- X1 + X2 + X3 + X4
      kron_true1 <- kronecker(create_data1$cov_col, create_data1$cov_row)
      kron_true2 <- kronecker(create_data2$cov_col, create_data2$cov_row)
      kron_true3 <- kronecker(create_data3$cov_col, create_data3$cov_row)
      kron_true4 <- kronecker(create_data4$cov_col, create_data4$cov_row)
      kron_true <- kron_true1 + kron_true2 + kron_true3 + kron_true4
      mu <- mu1 + mu2 + mu3 + mu4
    } 
    kron_true_inv <- solve(kron_true)
    kron_approx_rank1 <- best_kron_rank1(kron_true, p, q)
    
    quant_fda <- qchisq(0.99, p*n_basis)
    basis <- create.bspline.basis(rangeval = c(from,to), nbasis = n_basis, norder = 4)
    X_smooth_basis <- smooth.basis(argvals = as.numeric(dimnames(X)[[2]]), y = aperm(X, perm = c(2,1,3)), fdParobj = basis)
    coeffs <- coef(X_smooth_basis)
    X_smoothed <- X_smooth_basis$fd
    X_smooth <- aperm(eval.fd(evalarg = create_data1$eval_arg, fdobj = X_smoothed), c(2,1,3))
    
    ############################################
    # specific for non-separable setting
    ############################################
    X_vec <- t(apply(X,3,as.numeric))
    X_vec_clean <- t(apply(X_clean,3,as.numeric))
    
    mu_sample <- colMeans(X_vec_clean)
    Sigma_sample <- cov(X_vec_clean)
    
    tstart_mmle_clean <- Sys.time()
    par_mmle_clean <- mmle(X_clean)
    t_mmle_clean <- Sys.time()-tstart_mmle_clean
    kron_mmle <- kronecker(par_mmle_clean$cov_col, par_mmle_clean$cov_row)
    
    kron_norm_sample_clean <- norm(Sigma_sample - kron_true, "F")/norm(kron_true)
    kron_norm_mmle <- norm(kron_mmle - kron_true, "F")/norm(kron_true)
    ############################################
    ############################################
    
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
                          "lambda" = lambda,
                          "kron_sum_components" = i13,
                          "rank1_index" = kron_approx_rank1$sep_index, 
                          "rank1_norm" = kron_approx_rank1$sep_norm, 
                          "kron_norm_mmle" = kron_norm_mmle)
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
                             "lambda" = lambda,
                             "kron_sum_components" = i13,
                             "rank1_index" = kron_approx_rank1$sep_index, 
                             "rank1_norm" = kron_approx_rank1$sep_norm, 
                             "kron_norm_mmle" = kron_norm_mmle)
    res_mmle_fda <-  outlier_summary_skeleton(method = "mmle_fda", info = info_smooth)
    res_mmcd_fda <-  outlier_summary_skeleton(method = "mmcd_fda", info = info_smooth)
    
    try({
      md_actual <- mahalanobis(X_vec, as.numeric(mu), cov = kron_true, inverted = FALSE)
      outliers_actual <- which(md_actual > quant)
      res_actual <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_actual,
                                    method = "actual", info = info_raw, time = 0, 
                                    mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                    outlier_score = md_actual,
                                    eval_arg = NULL)
    })
    
    try({
      md_sample_clean <- mahalanobis(X_vec, center = mu_sample, cov = Sigma_sample, inverted = FALSE)
      outliers_sample_clean <- which(md_sample_clean > quant)
      res_sample_clean <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_sample_clean,
                                          method = "sample_clean", info = info_raw, time = 0, 
                                          mutlivariate_score = TRUE, functional = FALSE, reference_md = md_actual,
                                          outlier_score = md_sample_clean,
                                          eval_arg = NULL)
    })
    try({
      mmd_mmle_clean <- mmd(X,par_mmle_clean$mu, par_mmle_clean$cov_row_inv, par_mmle_clean$cov_col_inv, inverted = TRUE)
      outliers_mmle_clean <- which(mmd_mmle_clean > quant)
      
      res_mmle_clean <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_mmle_clean,
                                        method = "mmle_clean", info = info_raw, time = t_mmle_clean, 
                                        mutlivariate_score = TRUE, functional = FALSE, reference_md = md_actual,
                                        outlier_score = mmd_mmle_clean,
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
                                    mutlivariate_score = TRUE, functional = FALSE, reference_md = md_actual,
                                    outlier_score = mmd_mmle,
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
                                    mutlivariate_score = TRUE, functional = FALSE, reference_md = md_actual,
                                    outlier_score = mmd_mmcd,
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
                                         mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                         outlier_score = as.numeric(projdepth$MFDdepthX),
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
                                          mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                          outlier_score = as.numeric(sprojdepth$MFDdepthX),
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
                                          mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                          outlier_score = as.numeric(dprojdepth$MFDdepthX),
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
                                       mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                       outlier_score = ms_plot$dir_out_distance,
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
                                      mutlivariate_score = TRUE, functional = TRUE, reference_md = md_actual,
                                      outlier_score = mmd_mmle_fda,
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
                                      mutlivariate_score = TRUE, functional = TRUE, reference_md = md_actual,
                                      outlier_score = mmd_mmcd_fda,
                                      eval_arg = create_data$eval_arg)
    })
    
    rbind(res_mmle, res_mmcd, res_mmle_fda, res_mmcd_fda, 
          res_projdepth, res_sprojdepth, res_dprojdepth, res_ms_plot, 
          res_mmle_clean, res_sample_clean, res_actual)
  }
tend <- Sys.time()
tend-tstart
close(pb)
stopCluster(cl)
save("res",file = "raw_results/simulation_non_separable_p10.RData")
