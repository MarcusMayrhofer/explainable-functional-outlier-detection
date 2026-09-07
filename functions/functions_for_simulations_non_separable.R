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
                            mutlivariate_score = FALSE, functional = TRUE, reference_md = NULL,
                            outlier_score = NULL,
                            eval_arg){

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
  
  if(is.null(reference_md)){
    md_error <- NA
    md_error_rel <- NA
    md_cor_pearson <- NA
    md_cor_spearman <- NA
  } else{
    md_error <- sum((outlier_score - reference_md)^2)
    md_error_rel <- md_error/sum(reference_md^2)
    md_cor_pearson <- cor(outlier_score, reference_md, method = "pearson")
    md_cor_spearman <- cor(outlier_score, reference_md, method = "spearman")
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
    "md_error" = md_error,
    "md_error_rel" = md_error_rel,
    "md_cor_pearson" = md_cor_pearson,
    "md_cor_spearman" = md_cor_spearman
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
        "md_error" = NA,
        "md_error_rel" = NA,
        "md_cor_pearson" = NA,
        "md_cor_spearman" = NA
        )
      )
  )
}

