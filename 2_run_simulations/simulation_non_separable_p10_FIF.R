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
# source("../classification_and_clustering_functions.R")
source("../functions/multivariate_functional_data_functions.R")
source("../functions/helper_functions_simulation.R")
source("../functions/functions_for_simulations_non_separable.R")
source("../functions/rank_1_approx.R")


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
# mu1.1 = function(t, p) c(rep(f1(t, lambda = 1), p-1),f1(t, lambda = 0))
mu1.2 = function(t, p) (replicate(n = p, expr = f1(t, lambda = 0)))
mu2.0 = function(t, p) (replicate(n = p, expr = 4*t)) 
mu2.1 = function(t, p, lambda = 1) (replicate(n = p, expr = f2(t, lambda = lambda)))
# mu2.2 = function(t, p) c(rep(4*t, p-1),f2(t)) 
mu3.1 = function(t, p) (replicate(n = p, expr = f3(t)))
# mu3.2 = function(t, p) c(rep(4*t, p-1),f3(t)) 
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
# K <- K1.0
# K_outlier <- K1.0
from = 0 
to = 1
# Sigma <- Sigma1
# Sigma_outlier <- Sigma1
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
# n_basis_vec <- c(10,20,30)
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
  length(cov_function_type)*length(eigen_id_vec)*#length(n_basis_vec)*
  length(eps_coord_vec)*length(lambda_vec)*length(kappa_vec)*length(nu_vec)*length(kron_sum_vec)

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
  
  # ---- 3. Import NumPy ONCE ----
  np <- import("numpy")
  gc_py <- import("gc")
  
  # ---- 4. Load your module ----
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
  # foreach(i8 = 1:length(n_basis_vec), .combine = "rbind") %:%
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
    # n_basis <- n_basis_vec[i8]
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
    res_MFIF_Self <-  outlier_summary_skeleton(method = "MFIF_Self", info = info_raw)
    res_MFIF_gaussian_wavelets <-  outlier_summary_skeleton(method = "MFIF_gaussian_wavelets", info = info_raw)
    res_MFIF_Brownian <-  outlier_summary_skeleton(method = "MFIF_Brownian", info = info_raw)
    
    
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
                            # subsample_size = 16L,
                            # D = 'Brownian',
                            # D = 'Dyadic_indicator',
                            # D = 'SinusCosinus',
                            D = 'Self',
                            # D = 'gaussian_wavelets',
                            Dsize = 1000L,
                            innerproduct='auto1',
                            alpha = 0.5)
      scores_Self_py <- MFIF_Self$compute_paths()
      t_MFIF_Self <- Sys.time()-tstart_MFIF_Self
      scores_Self_R <- py_to_r(scores_Self_py)
      outliers_MFIF_Self <- order(scores_Self_R,decreasing = TRUE)[1:ceiling(n*outlier_percentage)]
      res_MFIF_Self <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_MFIF_Self,
                                       method = "MFIF_Self", info = info_raw, time = t_MFIF_Self,
                                       mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                       outlier_score = scores_Self_R,
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
                                         # subsample_size = 16L,
                                         D = 'gaussian_wavelets',
                                         Dsize = 1000L,
                                         innerproduct='auto1',
                                         alpha = 0.5)
      scores_gaussian_wavelets_py <- MFIF_gaussian_wavelets$compute_paths()
      t_MFIF_gaussian_wavelets <- Sys.time()-tstart_MFIF_gaussian_wavelets
      scores_gaussian_wavelets_R <- py_to_r(scores_gaussian_wavelets_py)
      outliers_MFIF_gaussian_wavelets <- order(scores_gaussian_wavelets_R,decreasing = TRUE)[1:ceiling(n*outlier_percentage)]
      res_MFIF_gaussian_wavelets <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_MFIF_gaussian_wavelets,
                                                    method = "MFIF_gaussian_wavelets", info = info_raw, time = t_MFIF_gaussian_wavelets,
                                                    mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                                    outlier_score = scores_gaussian_wavelets_R,
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
                                # subsample_size = 16L,
                                D = 'Brownian',
                                Dsize = 1000L,
                                innerproduct='auto1',
                                alpha = 0.5)
      scores_Brownian_py <- MFIF_Brownian$compute_paths()
      t_MFIF_Brownian <- Sys.time()-tstart_MFIF_Brownian
      scores_Brownian_R <- py_to_r(scores_Brownian_py)
      outliers_MFIF_Brownian <- order(scores_Brownian_R,decreasing = TRUE)[1:ceiling(n*outlier_percentage)]
      res_MFIF_Brownian <- outlier_summary(outlier_ind = outlier_index, outliers_flagged = outliers_MFIF_Brownian,
                                           method = "MFIF_Brownian", info = info_raw, time = t_MFIF_Brownian,
                                           mutlivariate_score = TRUE, functional = FALSE, reference_md = NULL,
                                           outlier_score = scores_Brownian_R,
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
    rbind(res_MFIF_Brownian, res_MFIF_Self, res_MFIF_gaussian_wavelets)
  }
tend <- Sys.time()
tend-tstart
close(pb)
stopCluster(cl)
save("res",file = "raw_results/simulation_non_separable_p10_FIF.RData")
