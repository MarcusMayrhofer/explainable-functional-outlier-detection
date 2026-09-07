md.mfd <- function(fdobj, par = NULL, method = "mmle"){
  X_coef <- coef(fdobj)

  if(is.null(par)){
    if(method == "mmcd"){
      par_X <- mmcd(X_coef)
    } else{
      par_X <- mmle(X_coef)
    }
  } else{
    par_X <- par
  }
  
  md <- robustmatrix::mmd(X = X_coef, 
                          mu = par_X$mu, 
                          cov_row = par_X$cov_row_inv, 
                          cov_col = par_X$cov_col_inv,
                          inverted = TRUE)
  md
}

shapley.mfd <- function(fdobj, par = NULL, method = "mmle", 
                        time_intervals = NULL, sum_coordinates = FALSE,
                        n_pc_kernel = NULL){
  X_coef <- coef(fdobj)
  if(!is.null(fdobj$basis)){
    X_basis <- fdobj$basis
  } else if(!is.null(fdobj$fd$basis)){
    X_basis <- fdobj$fd$basis
  } else{
    stop("No basis object contained in fdobj")
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
  
  if(!is.null(n_pc_kernel)){
    W <- inprod(X_basis, X_basis)
    pca_X <- pca.mfd(fdobj = fdobj, n_pc_kernel = n_pc_kernel, par = par_X)
    cov_kernel_regularized <- pca_X$pc_kernel$coefs %*% diag(pca_X$values_kernel[1:n_pc_kernel]) %*% t(pca_X$pc_kernel$coefs)
    cov_inv_kernel_regularized <- W %*% pca_X$pc_kernel$coefs %*% diag(1/pca_X$values_kernel[1:n_pc_kernel]) %*% t(pca_X$pc_kernel$coefs) %*% W
    par_X$cov_row <- cov_kernel_regularized
    par_X$cov_row_inv <- cov_inv_kernel_regularized
  }
  
  if(is.null(time_intervals)){
    shv <- robustmatrix::matrixShapley(X = X_coef, 
                                       mu = par_X$mu, 
                                       cov_row = par_X$cov_row_inv, 
                                       cov_col = par_X$cov_col_inv,
                                       inverted = TRUE, 
                                       type = "col")
  } else{
    if(length(time_intervals) == 1){
      time_seq <- seq(from = X_basis$rangeval[1], to = X_basis$rangeval[2], length.out = time_intervals + 1)
    } else{
      time_seq <- time_intervals
      time_intervals <- length(time_intervals) - 1
    }
    W <- inprod(X_basis, X_basis)
    W_inv <- solve(W)
    shv <- array(NA, dim = c(dim(X_coef)[2], length(time_seq) - 1, dim(X_coef)[3]))
    #get the first non zero decimal entry and add 1, yielding the number of decimal places diplayed for the names
    n_digits_tmp <- as.numeric(unlist(strsplit(formatC(median(diff(time_seq)), format = "e"), "e"))[2])
    n_digits <- ifelse(n_digits_tmp > 0, 0, abs(n_digits_tmp) + 1)
    dimnames(shv)[[2]] <- paste0(round(time_seq[-length(time_seq)],n_digits), " to ", round(time_seq[-1],n_digits))
    dimnames(shv)[[1]] <- dimnames(X_coef)[[2]]
    dimnames(shv)[[3]] <- dimnames(X_coef)[[3]]
    for(t_iter in 1:time_intervals){
      W_seq <- inprod(X_basis, X_basis, rng = time_seq[t_iter:(t_iter+1)])
      shv[,t_iter,] <- robustmatrix::matrixShapley(X = X_coef, 
                                                   mu = par_X$mu, 
                                                   cov_row = par_X$cov_row_inv%*% W_inv  %*% W_seq, 
                                                   cov_col = par_X$cov_col_inv,
                                                   inverted = TRUE, 
                                                   type = "col")
    }
    if(sum_coordinates){
      shv <- t(apply(shv,2,colSums))
    }
  }
  return(list("shv" = shv, 
              "par" = par_X))
}
