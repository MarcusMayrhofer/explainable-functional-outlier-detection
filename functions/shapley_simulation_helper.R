block_sum <- function(x, k) {
  stopifnot(is.numeric(x), k >= 1)
  
  if (is.null(dim(x))) {
    n <- length(x)
    sizes <- rep(n %/% k, k)
    remainder <- n %% k
    if (remainder > 0) sizes[1:remainder] <- sizes[1:remainder] + 1
    
    ends <- cumsum(sizes)
    starts <- c(1, head(ends, -1) + 1)
    
    sapply(seq_len(k), function(i) sum(x[starts[i]:ends[i]]))
    
  } else {
    apply(x, 1, function(row) block_sum(row, k)) |> t()
  }
}

sliding_sum <- function(x, window, step = 1) {
  stopifnot(is.numeric(x), window >= 1, step >= 1)
  
  if (is.null(dim(x))) {
    n <- length(x)
    if (window > n) stop("window length cannot exceed vector length")
    
    starts <- seq(1, n - window + 1, by = step)
    sapply(starts, function(i) sum(x[i:(i + window - 1)]))
    
  } else {
    t(apply(x, 1, function(row) sliding_sum(row, window, step)))
  }
}


AUC_func <- function(outlier_score, X_diff, n_interval = 1, absolute_value = FALSE, aggregate_time = FALSE, aggregate_coord = FALSE){
  AUC <- NA
  try({
    if(absolute_value){
      outlier_score <- abs(outlier_score)
    }
    labels <- matrix(0, nrow = dim(outlier_score)[1], ncol = dim(outlier_score)[2]) #Create ground-truth labels (1 = outlier, 0 = inlier)
    labels[which(X_diff>0)] <- 1
    if(aggregate_coord){
      labels <- colSums(labels)
      outlier_score <- colSums(outlier_score)
    }
    if(aggregate_time){
      labels <- rowSums(labels)
      outlier_score <- rowSums(outlier_score)
    }
    
    if(n_interval < 1){
      stop("Window must be larger than 1")
    } else if(n_interval == 1 | aggregate_time){
      labels <- +(labels>0)
      
      roc_obj <- pROC::roc(as.numeric(labels), as.numeric(outlier_score))
      AUC <- pROC::auc(roc_obj)
    } else{
      labels_block <- block_sum(labels, k = n_interval)
      labels <- +(labels_block>0)
      
      outlier_score_block <- block_sum(outlier_score,k = n_interval)
      roc_obj <- pROC::roc(as.numeric(labels), as.numeric(outlier_score_block))
      AUC <- pROC::auc(roc_obj)
    }
  })
  AUC
}

AUC_func_advanced <- function(outlier_score, X_diff, n_interval = c(1,10), method = "method",
                              types = c("time_coord", "time", "coord")){
  res <- data.frame()
  if(length(types) == 1 & "coord" %in% types){
    try({
      res <- rbind(res, data.frame(method, "type" = "coord", "absolute_score" = FALSE, "n_interval" = "none", 
                                   AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                  absolute_value = FALSE, aggregate_time = FALSE, aggregate_coord = FALSE)))
    })
    try({
      res <- rbind(res, data.frame(method, "type" = "coord", "absolute_score" = TRUE, "n_interval" = "none", 
                                   AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                  absolute_value = TRUE, aggregate_time = FALSE, aggregate_coord = FALSE)))
    })
  } else{
    for(j in 1:length(n_interval)){
      if(n_interval[j] == 1){
        interval_count = "all"
      } else{
        interval_count = n_interval[j]
      }
      if("time_coord" %in% types){
        try({
          res <- rbind(res, data.frame(method, "type" = "time_coord", "absolute_score" = FALSE, "n_interval" = interval_count, 
                                       AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                      absolute_value = FALSE, aggregate_time = FALSE, aggregate_coord = FALSE)))
        })
        try({
          res <- rbind(res, data.frame(method, "type" = "time_coord", "absolute_score" = TRUE, "n_interval" = interval_count, 
                                       AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                      absolute_value = TRUE, aggregate_time = FALSE, aggregate_coord = FALSE)))
        })
      }
      if("time" %in% types){
        try({
          res <- rbind(res, data.frame(method, "type" = "time", "absolute_score" = FALSE, "n_interval" = interval_count, 
                                       AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                      absolute_value = FALSE, aggregate_time = FALSE, aggregate_coord = TRUE)))
        })
        try({
          res <- rbind(res, data.frame(method, "type" = "time", "absolute_score" = TRUE, "n_interval" = interval_count, 
                                       AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                      absolute_value = TRUE, aggregate_time = FALSE, aggregate_coord = TRUE)))
        })
      }
      if("coord" %in% types){
        try({
          res <- rbind(res, data.frame(method, "type" = "coord", "absolute_score" = FALSE, "n_interval" = "none", 
                                       AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                      absolute_value = FALSE, aggregate_time = TRUE, aggregate_coord = FALSE)))
        })
        try({
          res <- rbind(res, data.frame(method, "type" = "coord", "absolute_score" = TRUE, "n_interval" = "none", 
                                       AUC = AUC_func(outlier_score, X_diff, n_interval = n_interval[j],
                                                      absolute_value = TRUE, aggregate_time = TRUE, aggregate_coord = FALSE)))
        })
      }
    }
  }
  res
}