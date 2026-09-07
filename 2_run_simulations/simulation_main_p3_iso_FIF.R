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


# Setup of Functional Isolation Forest
##################################
require(reticulate)
mfif_path <- file.path(this.path::this.dir(), "..", "MFIF.py")
###################################



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
eigen_shift_vec <- NA
# eigen_shift_vec <- c(15,22.5,30) #use this for eigen_id = 1, p = 3
# eigen_shift_vec <- c(12,15,18) #use this for eigen_id = 1, p = 10
# eigen_shift_vec <- c(8,10,12) #use this for eigen_id = 1, p = 50
# eigen_shift_vec <- c(5,8,11) #use this for eigen_id = 10, p = 3
# eigen_shift_vec <- c(2,3,4) #use this for eigen_id = 10, p = 10
# eigen_shift_vec <- c(1.5,2,2.5) #use this for eigen_id = 10, p = 50
eigen_id_vec <- NA
# eigen_id_vec <- c(1,10)
# eigen_id_vec <- c(1)
cov_function_type = c("M", "OU")
# n_basis_vec <- c(10,20,30)
eps_coord_vec <- c(0.1,0.5,1)
#isolated outleirs
lambda_vec <- NA
lambda_vec <- c(2,5,8)/10
#covariance outliers
kappa_vec <- NA
# kappa_vec <- c(7, 10, 15)
nu_vec <- NA
# nu_vec <- c(0.1,0.2,0.5)


npar <- n_run*length(p_vec)*length(eps_vec)*length(dim_vec)*length(eigen_shift_vec)*
  length(cov_function_type)*length(eigen_id_vec)*
  length(eps_coord_vec)*length(lambda_vec)*length(kappa_vec)*length(nu_vec)

cl <- makeCluster(120)
registerDoSNOW(cl)

clusterExport(cl, "mfif_path")

clusterEvalQ(cl, {
  
  # ---- 1. Disable BLAS/OpenMP threading ----
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    MKL_DYNAMIC = "FALSE"
  )
  
  # ---- 2. Initialize Python ----
  library(reticulate)
  
  # ---- 3. Import NumPy ----
  np <- import("numpy")
  gc_py <- import("gc")
  
  # ---- 4. Load module ----
  source_python(mfif_path)
  
  NULL
})


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
      if(!any(is.na(eigen_id)) | !any(is.na(eigen_shift))){
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
    res_MFIF_Self <-  outlier_summary_skeleton(method = "MFIF_Self", info = info_raw)
    res_MFIF_gaussian_wavelets <-  outlier_summary_skeleton(method = "MFIF_gaussian_wavelets", info = info_raw)
    res_MFIF_Brownian <-  outlier_summary_skeleton(method = "MFIF_Brownian", info = info_raw)
    
    
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
    try({
      ##################################
      # Run FIF 
      ##################################
      
      X_py <- np$array(
        aperm(X,c(3,1,2)),
        dtype = "float64"
      )
      
      time_py <- np$array(
        1:q,
        dtype = "float64"
      )
      
      tstart_MFIF_Self <- Sys.time()
      MFIF_Self = MFIForest(X_py, 
                            time=time_py, 
                            ntrees = 100L, 
                            D = 'Self',
                            Dsize = 1000L, 
                            innerproduct='auto1', 
                            alpha = 0.5)
      scores_Self_py <- MFIF_Self$compute_paths()
      t_MFIF_Self <- Sys.time()-tstart_MFIF_Self
      scores_Self_R <- py_to_r(scores_Self_py)
      outliers_MFIF_Self <- order(scores_Self_R,decreasing = TRUE)[1:ceiling(n*outlier_percentage)]
      par_MFIF_Self <- mmle(X[,,if(length(outliers_MFIF_Self) == 0){TRUE} else{-outliers_MFIF_Self}])
      res_MFIF_Self <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_MFIF_Self,
                                       method = "MFIF_Self", info = info_raw, time = t_MFIF_Self, 
                                       mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                       outlier_score = scores_Self_R,
                                       mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                       mu_function_est = par_MFIF_Self$mu, 
                                       cov_matrix_est = par_MFIF_Self$cov_row, 
                                       cov_function_est = par_MFIF_Self$cov_col,
                                       eval_arg = NULL)
    })
    try({
      rm(X_py)
      rm(time_py)
      rm(MFIF_Self)
      rm(scores_Self_py)
      gc()
      gc_py$collect()
    })
    try({
      ##################################
      # Run FIF 
      ##################################
      
      X_py <- np$array(
        aperm(X,c(3,1,2)),
        dtype = "float64"
      )
      
      time_py <- np$array(
        1:q,
        dtype = "float64"
      )
      
      tstart_MFIF_gaussian_wavelets <- Sys.time()
      MFIF_gaussian_wavelets = MFIForest(X_py, 
                                         time=time_py, 
                                         ntrees = 100L, 
                                         D = 'gaussian_wavelets',
                                         Dsize = 1000L, 
                                         innerproduct='auto1', 
                                         alpha = 0.5)
      scores_gaussian_wavelets_py <- MFIF_gaussian_wavelets$compute_paths()
      t_MFIF_gaussian_wavelets <- Sys.time()-tstart_MFIF_gaussian_wavelets
      scores_gaussian_wavelets_R <- py_to_r(scores_gaussian_wavelets_py)
      outliers_MFIF_gaussian_wavelets <- order(scores_gaussian_wavelets_R,decreasing = TRUE)[1:ceiling(n*outlier_percentage)]
      par_MFIF_gaussian_wavelets <- mmle(X[,,if(length(outliers_MFIF_gaussian_wavelets) == 0){TRUE} else{-outliers_MFIF_gaussian_wavelets}])
      res_MFIF_gaussian_wavelets <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_MFIF_gaussian_wavelets,
                                                    method = "MFIF_gaussian_wavelets", info = info_raw, time = t_MFIF_gaussian_wavelets, 
                                                    mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                                    outlier_score = scores_gaussian_wavelets_R,
                                                    mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                                    mu_function_est = par_MFIF_gaussian_wavelets$mu, 
                                                    cov_matrix_est = par_MFIF_gaussian_wavelets$cov_row, 
                                                    cov_function_est = par_MFIF_gaussian_wavelets$cov_col,
                                                    eval_arg = NULL)
    })
    try({
      rm(X_py)
      rm(time_py)
      rm(MFIF_gaussian_wavelets)
      rm(scores_gaussian_wavelets_py)
      gc()
      gc_py$collect()
    })
    try({
      ##################################
      # Run FIF 
      ##################################
      
      X_py <- np$array(
        aperm(X,c(3,1,2)),
        dtype = "float64"
      )
      
      time_py <- np$array(
        1:q,
        dtype = "float64"
      )
      
      tstart_MFIF_Brownian <- Sys.time()
      MFIF_Brownian = MFIForest(X_py, 
                                time=time_py, 
                                ntrees = 100L, 
                                D = 'Brownian',
                                Dsize = 1000L, 
                                innerproduct='auto1', 
                                alpha = 0.5)
      scores_Brownian_py <- MFIF_Brownian$compute_paths()
      t_MFIF_Brownian <- Sys.time()-tstart_MFIF_Brownian
      scores_Brownian_R <- py_to_r(scores_Brownian_py)
      outliers_MFIF_Brownian <- order(scores_Brownian_R,decreasing = TRUE)[1:ceiling(n*outlier_percentage)]
      par_MFIF_Brownian <- mmle(X[,,if(length(outliers_MFIF_Brownian) == 0){TRUE} else{-outliers_MFIF_Brownian}])
      res_MFIF_Brownian <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_MFIF_Brownian,
                                           method = "MFIF_Brownian", info = info_raw, time = t_MFIF_Brownian, 
                                           mutlivariate_score = TRUE, functional = FALSE, reference_method = res_mmle_clean,
                                           outlier_score = scores_Brownian_R,
                                           mu = mu, cov_row = cov_row, cov_col = cov_col, cov_col_projected = cov_col_projected,
                                           mu_function_est = par_MFIF_Brownian$mu, 
                                           cov_matrix_est = par_MFIF_Brownian$cov_row, 
                                           cov_function_est = par_MFIF_Brownian$cov_col,
                                           eval_arg = NULL)
    })
    try({
      rm(X_py)
      rm(time_py)
      rm(MFIF_Brownian)
      rm(scores_Brownian_py)
      gc()
      gc_py$collect()
    })
    rbind(res_MFIF_Self,res_MFIF_gaussian_wavelets,res_MFIF_Brownian)
  }
tend <- Sys.time()
tend-tstart
close(pb)
stopCluster(cl)
save("res",file = "raw_results/simulation_res_p3_iso_FIF.RData")
