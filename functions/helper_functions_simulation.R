msp <- function(from = 0, to = 1, p = 2, q = to, n = 10,
                mu = function(t) {c(0,0)},
                cov_row = diag(1, nrow = p, ncol = p),
                K = function(s, t) {min(s, t)}, eigen_id = 1, gamma = 0,
                type = "normal", df = NULL, mu_rand = FALSE) {
  p = nrow(cov_row)
  
  t <- seq(from = from, to = to, length.out = q)
  cov_col <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  cov_row_eigen <- eigen(cov_row)
  cov_col_eigen <- eigen(cov_col)
  
  cov_row_sqrt <- cov_row_eigen$vectors %*% diag(sqrt(cov_row_eigen$values)) %*% t(cov_row_eigen$vectors)
  cov_col_sqrt <- cov_col_eigen$vectors %*% diag(sqrt(cov_col_eigen$values)) %*% t(cov_col_eigen$vectors)
  
  mu_array <- NULL
  if(mu_rand){
    paths <- array(dim = c(p,q,n))
    mu_array <- array(dim = c(p,q,n))
    if(type == "t"){
      for(i in 1:n){
        mu_mat <- matrix(mu(t = t, p = p), nrow = p, ncol = q, byrow = TRUE)
        mu_array[,,i] <- mu_mat
        paths[,,i] <- MixMatrix::rmatrixt(n = 1, df = df, mean = mu_mat, L = cov_row, R = cov_col)
      }
    } else {
      for(i in 1:n){
        mu_mat <- matrix(mu(t = t, p = p), nrow = p, ncol = q, byrow = TRUE)
        mu_array[,,i] <- mu_mat
        paths[,,i] <- robustmatrix::rmatnorm(n = 1, mu = mu_mat, cov_row = cov_row_sqrt, cov_col = cov_col_sqrt, sqrt = TRUE)
      }
    }
    mu_mat <- apply(paths, 1:2, FUN = mean)
  } else{
    mu_mat <- matrix(mu(t = t, p = p), nrow = p, ncol = q, byrow = TRUE)
    if(type == "t"){
      paths <- MixMatrix::rmatrixt(n = n, df = df, mean = mu_mat, L = cov_row, R = cov_col)
    } else {
      paths <- robustmatrix::rmatnorm(n = n, mu = mu_mat, cov_row = cov_row_sqrt, cov_col = cov_col_sqrt, sqrt = TRUE)
    }
  }
  
  if(gamma != 0){
    eigen_vec <- cov_col_eigen$vectors[,eigen_id]
    shift_mat <- matrix(rep(eigen_vec,p), byrow = TRUE, nrow = p)
    paths <- sweep(paths, MARGIN = 1:2, STATS = shift_mat*gamma, FUN = "+")
  }
  
  
  if(n == 1){
    dimnames(paths) <- list(NULL, t)
    
  } else{
    dimnames(paths) <- list(NULL, t, NULL)
  }
  
  return(list("X" = paths, "mu" = mu_mat, "cov_row" = cov_row, "cov_col" = cov_col, "eval_arg" = t, "mu_array" = mu_array))
}

plot_cov <- function(cov_matrix, eval_arg = NULL, plot_cor = FALSE, guide = "none"){
  if(is.null(eval_arg)){
    eval_arg <- seq(from = 0, to = 1, length.out = nrow(cov_matrix))
  }
  if(plot_cor){
    cov_matrix <- cov2cor(cov_matrix)
  } else{
    cov_matrix <- cov_matrix/cov_matrix[1,1]
  }
  dimnames(cov_matrix) <- list(eval_arg, eval_arg)
  cov_matrix_long <- data.frame(cov_matrix, check.names = FALSE) %>% 
    rownames_to_column() %>% 
    pivot_longer(-rowname) %>% 
    mutate(s = as.numeric(rowname), t = as.numeric(name)) %>%
    dplyr::select(-c(rowname, name))
  
  plt_cov <- ggplot(cov_matrix_long, aes(x = s, y = t, fill = value)) + 
    geom_tile() + 
    theme_classic() + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_fixed()
  
  if(plot_cor){
    plt_cov <- plt_cov + 
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", guide = guide, midpoint = 0, limits = c(-1,1))
  } else{
    plt_cov <- plt_cov + 
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", guide = guide, midpoint = 0)
  }
  plt_cov
}


subset_if_outliers <- function(X, outliers) {
  if (length(outliers) > 0) X[,,-outliers] else X
}

#######################################################
# This are functions from the fdaoutlier package needed to adapt the msplot function below
#######################################################

croux_hesbroeck_asymptotic <- function(n, dimension){
  h <- floor((n+dimension+1)/2)
  alpha <- (n-h)/n
  q_alpha <- qchisq(1-alpha, dimension)
  c_alpha <- (1 - alpha)/pchisq(q_alpha, dimension + 2)
  c2 <- -pchisq(q_alpha, dimension+2)/2
  c3 <- -pchisq(q_alpha, dimension + 4)/2
  c4 <- 3*c3
  b1 <- c_alpha*(c3-c4)/(1-alpha)
  b2 <- 0.5 + c_alpha/(1-alpha)*(c3-q_alpha*(c2+(1-alpha)/2)/dimension)
  v1  <- (1-alpha)*(b1^2)*
    (alpha*(c_alpha*q_alpha/dimension-1)^2-1)-2*c3*c_alpha^2*
    (3*(b1-dimension*b2)^2+(dimension+2)*b2*(2*b1-dimension*b2))
  v2 <- n*(b1*(b1-dimension*b2)*(1-alpha))^2*c_alpha^2
  v <- v1/v2
  m_asy <- 2/(c_alpha^2*v)
  m <- m_asy*exp(0.725-0.00663*dimension-0.078*log(n))
  if (m < dimension){ #if m is >= dimension, then line 18 works, if not change m to m_asy
    m <- m_asy
  }
  a1 <- rchisq(10000,dimension + 2)
  a2 <- rchisq(10000,dimension, h/n)
  c <- sum(a1 < a2)/(10000*h/n)
  factors <- c * (m - dimension + 1)/(dimension * m)
  cutoff <- qf(0.993, dimension, m - dimension + 1)
  list(factor1 = factors, factor2 = cutoff)
}

hardin_factor_numeric <- function(n, dimension){
  if (dimension == 2){
    if(n < 1000){
      k  <- floor(n / 5) + 1
      factor1 = hardin_factor_numeric_dimen_2$factor1[k]
      factor2 = hardin_factor_numeric_dimen_2$factor2[k]
    }
    if (n >= 1000){
      asymp_result <- croux_hesbroeck_asymptotic(n = n, dimension = dimension)
      factor1  <- asymp_result$factor1
      factor2 = asymp_result$factor2
    }
    
  } else if (dimension == 3){
    if (n < 1000){
      k  <- floor(n / 5) + 1
      factor1  <- hardin_factor_numeric_dimen_3$factor1[k]
      factor2  <- hardin_factor_numeric_dimen_3$factor2[k]
    }
    if (n >= 1000){
      asymp_result  <- croux_hesbroeck_asymptotic(n = n, dimension = dimension)
      factor1 <- asymp_result$factor1
      factor2 <- asymp_result$factor2
    }
  } else if (dimension > 3){
    asymp_result  <- croux_hesbroeck_asymptotic(n = n, dimension = dimension)
    factor1 <- asymp_result$factor1
    factor2 <- asymp_result$factor2
  } else{
    stop("Argument \'dimension\' must be greater than or equal to 2.")
  }
  return(list(factor1 = factor1, factor2 = factor2))
}


#######################################################
# This is a version of the MS plot that also outputs the distance from the dir_out function which is used as the outlier cutoff.
#######################################################
msplot_adv <- function (dts, data_depth = c("random_projections"), n_projections = 200, 
                        seed = NULL, return_mvdir = TRUE, plot = TRUE, plot_title = "Magnitude Shape Plot", 
                        title_cex = 1.5, show_legend = T, ylabel = "VO", xlabel) 
{
  data_dim <- dim(dts)
  n <- data_dim[1]
  dir_result <- dir_out(dts, data_depth = data_depth, n_projections = n_projections, 
                        seed = seed)
  if (length(data_dim) == 2) {
    dist <- dir_result$distance
    rocke_factors <- hardin_factor_numeric(n, 2)
    rocke_factor1 <- rocke_factors$factor1
    rocke_cutoff <- rocke_factors$factor2
    cutoff_value <- rocke_cutoff/rocke_factor1
    outliers_index <- which(dist > cutoff_value)
    median_curve <- which.min(dist)
    if (plot) {
      myx <- dir_result$mean_outlyingness
      myy <- dir_result$var_outlyingness
      if (missing(xlabel)) 
        xlabel <- "MO"
    }
  }
  else if (length(data_dim) == 3) {
    d <- data_dim[3]
    rocke_factors <- hardin_factor_numeric(n = n, dimension = d + 
                                             1)
    rocke_factor1 <- rocke_factors$factor1
    rocke_cutoff <- rocke_factors$factor2
    cutoff_value <- rocke_cutoff/rocke_factor1
    outliers_index <- which(dir_result$distance > cutoff_value)
    median_curve <- which.min(dir_result$distance)
    if (plot) {
      myx <- sqrt(rowSums(dir_result$mean_outlyingness^2, 
                          na.rm = T))
      myy <- dir_result$var_outlyingness
      if (missing(xlabel)) 
        xlabel <- "||MO||"
    }
  }
  if (plot) {
    plot(myx, myy, type = "n", xlab = xlabel, ylab = ylabel, 
         xlim = range(myx) + c(-sd(myx), 1.5 * sd(myx)), ylim = range(myy) + 
           c(-0.2 * sd(myy), 1 * sd(myy)), axes = F, col.lab = "gray20")
    axis(1, col = "white", col.ticks = "grey61", lwd.ticks = 0.5, 
         tck = -0.025, cex.axis = 0.9, col.axis = "gray30")
    axis(2, col = "white", col.ticks = "grey61", lwd.ticks = 0.5, 
         tck = -0.025, cex.axis = 0.9, col.axis = "gray30")
    grid(col = "grey75", lwd = 0.3)
    box(col = "grey51")
    if (length(outliers_index > 0)) {
      points(myx[-outliers_index], myy[-outliers_index], 
             bg = "gray60", pch = 21)
      points(myx[outliers_index], myy[outliers_index], 
             pch = 3)
    }
    else {
      points(myx, myy, bg = "gray60", pch = 21)
    }
    mtext(plot_title, 3, adj = 0.5, line = 1, cex = title_cex, 
          col = "gray20")
    if (show_legend) {
      legend("topright", legend = c("normal", "outlier"), 
             pch = c(21, 3), cex = 1, pt.bg = "gray60", col = "gray0", 
             text.col = "gray30", bty = "n", box.lwd = 0.1, 
             xjust = 0, inset = 0.01)
    }
  }
  if (return_mvdir) {
    return(list(outliers = outliers_index, median_curve = median_curve, 
                mean_outlyingness = dir_result$mean_outlyingness, 
                var_outlyingness = dir_result$var_outlyingness,
                dir_out_distance = dir_result$distance))
  }
  else {
    return(list(outliers = outliers_index, median_curve = median_curve, dir_out_distance = dir_result$distance))
  }
}
