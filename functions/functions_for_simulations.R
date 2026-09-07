frob_error <- function(cov_row, cov_col, cov_row_est, cov_col_est){
  cov_kron <- kronecker(cov_col, cov_row)
  cov_kron_est <- kronecker(cov_col_est, cov_row_est)
  norm(cov_kron - cov_kron_est, type = "F") / norm(cov_kron, type = "F")
}

frob_error_vec <- function(cov_row, cov_col, cov_est){
  cov_kron <- kronecker(cov_col, cov_row)
  norm(cov_kron - cov_est, type = "F") / norm(cov_kron, type = "F")
}

outlier_summary <- function(outlier_ind, outliers_flagged, method = "default", info = NULL, time = NA, 
                            mutlivariate_score = FALSE, functional = TRUE, reference_method = NULL,
                            outlier_score = NULL,
                            mu, cov_row, cov_col, cov_col_projected,
                            mu_function_est, cov_matrix_est, cov_function_est, eval_arg){
  
  if(functional){
    mu_function_est_eval <- t(eval.fd(evalarg = create_data$eval_arg, fdobj = mu_function_est))
    cov_function_est_eval <- eval.bifd(bifd = cov_function_est, create_data$eval_arg, create_data$eval_arg)
  } else{
    mu_function_est_eval <- mu_function_est
    cov_function_est_eval <- cov_function_est
  }
  
  score_mu <- mean((mu - mu_function_est_eval)^2)
  score_cov_matrix <- mean((cov_row/cov_row[1,1] - cov_matrix_est/cov_matrix_est[1,1])^2)
  score_cov_function <- mean((sum(diag(cov_row)) * cov_col - sum(diag(cov_matrix_est)) * cov_function_est_eval)^2)
  
  cov_function_est_eval_k5 <- cov_col_projected$proj_k5%*%cov_function_est_eval%*%cov_col_projected$proj_k5
  cov_function_est_eval_k10 <- cov_col_projected$proj_k10%*%cov_function_est_eval%*%cov_col_projected$proj_k10
  cov_function_est_eval_k20 <- cov_col_projected$proj_k20%*%cov_function_est_eval%*%cov_col_projected$proj_k20
  score_cov_function_k5 <- mean((sum(diag(cov_row)) * cov_col_projected$cov_k5 - sum(diag(cov_matrix_est)) * cov_function_est_eval_k5)^2)
  score_cov_function_k10 <- mean((sum(diag(cov_row)) * cov_col_projected$cov_k10 - sum(diag(cov_matrix_est)) * cov_function_est_eval_k10)^2)
  score_cov_function_k20 <- mean((sum(diag(cov_row)) * cov_col_projected$cov_k20 - sum(diag(cov_matrix_est)) * cov_function_est_eval_k20)^2)
  score_kronecker <- NA
  if(mutlivariate_score == TRUE){
    score_cov_function_multi <- mean((kronecker(cov_row, cov_col) - kronecker(cov_matrix_est,cov_function_est_eval))^2)
    score_cov_function_multi_k5 <- mean((kronecker(cov_row, cov_col_projected$cov_k5) - kronecker(cov_matrix_est,cov_function_est_eval_k5))^2)
    score_cov_function_multi_k10 <- mean((kronecker(cov_row, cov_col_projected$cov_k10) - kronecker(cov_matrix_est,cov_function_est_eval_k10))^2)
    score_cov_function_multi_k20 <- mean((kronecker(cov_row, cov_col_projected$cov_k20) - kronecker(cov_matrix_est,cov_function_est_eval_k20))^2)
  }
  
  
  n_out <- length(outlier_ind)
  flagged <- outliers_flagged
  n_flagged <- length(flagged)
  
  
  TP <- length(intersect(outlier_ind, flagged)) #(correctly flagged outliers)
  FN <- length(setdiff(outlier_ind, flagged)) #(missed true outliers)
  FP <- length(setdiff(flagged, outlier_ind)) #(normal points incorrectly flagged)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  fscore <- ifelse((2*TP + FP + FN) > 0,
                   (2*TP)/(2*TP + FP + FN),
                   NA)
  
  if(is.null(outlier_score)){
    AUC <- NA
  } else {
    labels <- integer(length(outlier_score)) #Create ground-truth labels (1 = outlier, 0 = inlier)
    labels[outlier_ind] <- 1
    
    #Compute ROC and AUC
    roc_obj <- pROC::roc(labels, outlier_score) 
    AUC <- pROC::auc(roc_obj)
  } 
  
  score_mu_rel <- NA
  score_cov_matrix_rel <- NA
  score_cov_function_rel <- NA
  score_cov_function_rel_k5 <- NA
  score_cov_function_rel_k10 <- NA
  score_cov_function_rel_k20 <- NA
  score_score_cov_function_multi_rel <- NA
  score_score_cov_function_multi_rel_k5 <- NA
  score_score_cov_function_multi_rel_k10 <- NA
  score_score_cov_function_multi_rel_k20 <- NA
  if(!is.null(reference_method)){
    score_mu_rel <- (score_mu)/reference_method$score_mu
    score_cov_matrix_rel <- (score_cov_matrix)/reference_method$score_cov_matrix
    score_cov_function_rel <- (score_cov_function)/reference_method$score_cov_function
    score_cov_function_rel_k5 <- (score_cov_function_k5)/reference_method$score_cov_function_k5
    score_cov_function_rel_k10 <- (score_cov_function_k10)/reference_method$score_cov_function_k10
    score_cov_function_rel_k20 <- (score_cov_function_k20)/reference_method$score_cov_function_k20
    if(mutlivariate_score == TRUE){
      score_score_cov_function_multi_rel <- (score_cov_function_multi)/reference_method$score_cov_function_multi
      score_score_cov_function_multi_rel_k5 <- (score_cov_function_multi_k5)/reference_method$score_cov_function_multi_k5
      score_score_cov_function_multi_rel_k10 <- (score_cov_function_multi_k10)/reference_method$score_cov_function_multi_k10
      score_score_cov_function_multi_rel_k20 <- (score_cov_function_multi_k20)/reference_method$score_cov_function_multi_k20
    }
  }
  # elapsed time, always in seconds
  time <- if (inherits(time, "difftime")) as.numeric(time, units = "secs") else as.numeric(time)

  tmp <- data.frame(
    method,
    "time" = time,
    info,
    "n_flagged" = n_flagged, 
    "TP" = TP,
    "FP" = FP,
    "FN" = FN,
    "precision" = precision, 
    "recall" = recall, 
    "F-score" = fscore,
    "AUC" = AUC,
    "score_mu" = score_mu,
    "score_cov_matrix" = score_cov_matrix,
    "score_cov_function" = score_cov_function,
    "score_cov_function_k5" = score_cov_function_k5,
    "score_cov_function_k10" = score_cov_function_k10,
    "score_cov_function_k20" = score_cov_function_k20,
    "score_cov_function_multi" = score_cov_function_multi, 
    "score_cov_function_multi_k5" = score_cov_function_multi_k5, 
    "score_cov_function_multi_k10" = score_cov_function_multi_k10, 
    "score_cov_function_multi_k20" = score_cov_function_multi_k20, 
    "score_mu_rel" = score_mu_rel, 
    "score_cov_matrix_rel" = score_cov_matrix_rel, 
    "score_cov_function_rel" = score_cov_function_rel, 
    "score_cov_function_rel_k5" = score_cov_function_rel_k5, 
    "score_cov_function_rel_k10" = score_cov_function_rel_k10, 
    "score_cov_function_rel_k20" = score_cov_function_rel_k20, 
    "score_score_cov_function_multi_rel" = score_score_cov_function_multi_rel,
    "score_score_cov_function_multi_rel_k5" = score_score_cov_function_multi_rel_k5,
    "score_score_cov_function_multi_rel_k10" = score_score_cov_function_multi_rel_k10,
    "score_score_cov_function_multi_rel_k20" = score_score_cov_function_multi_rel_k20
  )
  tmp
}

# empty result row, used when a method is not run or fails; all columns atomic
outlier_summary_skeleton <- function(method = "default", info = NULL){
  data.frame(
    "method" = method,
    "time" = NA_real_,
    info,
    t(c("n_flagged" = NA_real_,
        "TP" = NA,
        "FP" = NA,
        "FN" = NA,
        "precision" = NA,
        "recall" = NA,
        "F-score" = NA,
        "AUC" = NA,
        "score_mu" = NA,
        "score_cov_matrix" = NA,
        "score_cov_function" = NA,
        "score_cov_function_k5" = NA,
        "score_cov_function_k10" = NA,
        "score_cov_function_k20" = NA,
        "score_cov_function_multi" = NA, 
        "score_cov_function_multi_k5" = NA, 
        "score_cov_function_multi_k10" = NA, 
        "score_cov_function_multi_k20" = NA, 
        "score_mu_rel" = NA, 
        "score_cov_matrix_rel" = NA, 
        "score_cov_function_rel" = NA, 
        "score_cov_function_rel_k5" = NA, 
        "score_cov_function_rel_k10" = NA, 
        "score_cov_function_rel_k20" = NA, 
        "score_score_cov_function_multi_rel" = NA,
        "score_score_cov_function_multi_rel_k5" = NA,
        "score_score_cov_function_multi_rel_k10" = NA,
        "score_score_cov_function_multi_rel_k20" = NA))
  )
}

