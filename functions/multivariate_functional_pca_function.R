pca.mfd <- function(fdobj, n_pc_kernel = NULL, n_pc_var = NULL, par = NULL, method = "mmle", explained_variance = TRUE){
  X_coef <- coef(fdobj)
  if(!is.null(fdobj$basis)){
    X_basis <- fdobj$basis
  } else if(!is.null(fdobj$fd$basis)){
    X_basis <- fdobj$fd$basis
  } else{
    stop("No basis object contained in fdobj")
  }
  
  if(is.null(n_pc_kernel)){
    n_pc_kernel <- dim(X_coef)[1]
  } else{
    n_pc_kernel <- min(n_pc_kernel, dim(X_coef)[1])
  }
  if(is.null(n_pc_var)){
    n_pc_var <- dim(X_coef)[2]
  } else{
    n_pc_var <- min(n_pc_var, dim(X_coef)[2])
  }
  
  if(is.null(par)){
    if(method == "mmcd"){
      par_X <- mmcd(X_coef)
    } else{
      par_X <- mmle(X_coef)
    }
  } else{
    par_X <- par
  }
  
  W <- inprod(X_basis, X_basis)
  W_eigen <- eigen(W)
  W_sqrt <- W_eigen$vectors%*%diag(sqrt(W_eigen$values))%*%t(W_eigen$vectors)
  W_sqrt_inv <- solve(W_sqrt)
  
  # We have to account for the inner products of the basis functions, W, when we compute the eigen decomposition.
  # In notation of the paper this is the eigendecomposition 
  # $\bm{W}^{1/2}\bm{\Sigma}^{\col}\bm{W}^{1/2}= \sum_{i = 1}^m \lambda_i\textbf{u}_i$, 
  # with $\textbf{u}_i = \bm{W}^{1/2}\textbf{b}_i$.
  eigen_kernel <- eigen(W_sqrt %*% par_X$cov_row %*% W_sqrt) 
  # When we want the eigen vectors (functions) of the kernel (cov. matrix) we need to account again for W.
  # We want the PC functions/vectors $\textbf{b}_i$ of the covariance kernel/matrix alone, not including the W matrix.
  pc_coef_kernel <- W_sqrt_inv %*% eigen_kernel$vectors
  pc_kernel <- fd(coef = pc_coef_kernel[,1:n_pc_kernel], basisobj = X_basis)
  
  eigen_var <- eigen(par_X$cov_col)
  
  explained_variance_list <- NULL
  if(explained_variance){
    #explained variance by PCs of kernel
    expv_kernel <- cumsum(eigen_kernel$values)/sum(eigen_kernel$values)
    #explained variance by PCs of covariance between variables
    expv_var <- cumsum(eigen_var$values)/sum(eigen_var$values)
    
    expv <- data.frame(matrix(apply(expand.grid(expv_kernel,expv_var), 1, prod), ncol = dim(X_coef)[2])) 
    colnames(expv) <- 1:dim(X_coef)[2]
    expv_long <- expv %>% 
      rownames_to_column(var = "n_pc_kernel") %>% 
      pivot_longer(-n_pc_kernel, names_to = "n_pc_var", values_to = "explained_variance") %>%
      mutate(n_pc_kernel = as.numeric(n_pc_kernel),
             n_pc_var = as.numeric(n_pc_var),
             n_pc_total = n_pc_var * n_pc_kernel)
    
    expv_long_tmp <- expv_long %>% 
      group_by(n_pc_total) %>% 
      mutate(explained_variance_best = max(explained_variance)) %>%
      arrange(n_pc_total) %>%
      filter(explained_variance == explained_variance_best) %>%
      dplyr::select(-explained_variance_best)
    
    index_expv <- NULL
    for(i in 1:nrow(expv_long_tmp)){
      if(all(expv_long_tmp$explained_variance[1:max((i-1),1)] <= expv_long_tmp$explained_variance[i])){
        index_expv <- c(index_expv,i)
      }
    }
    
    expv_best <- expv_long_tmp[index_expv,]
    explained_variance_list <- list("best" = expv_best,
                                    "kernel" = expv_kernel,
                                    "var" = expv_var,
                                    "total" = expv_long)
  }
  
  mu_function <- fda::fd(coef = par_X$mu, basisobj = basis)
  cov_matrix <- par_X$cov_col
  cov_function <- fda::bifd(par_X$cov_row, basis, basis)

  return(list("pc_kernel" = pc_kernel, 
              "pc_var" = eigen_var$vectors[,1:n_pc_var],
              "values_kernel" = eigen_kernel$values,
              "values_var" = eigen_var$values,
              "mu_function" = mu_function,
              "cov_matrix" = cov_matrix,
              "cov_function" = cov_function,
              "par" = par_X,
              "basis" = X_basis,
              "explained_variance" = explained_variance_list,
              "W" = W,
              "W_sqrt" = W_sqrt,
              "W_sqrt_inv" = W_sqrt_inv))
}
