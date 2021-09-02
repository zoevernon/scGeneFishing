#' 
#' getCutoff
#' 
#' Function to get CFR cutoff from \code{geneFishing()} results. This method 
#' assumes the model is decreasing, flat and increasing. It detects the 
#' boundary ofmthe decreasing part and the flat part, and the boundary of the 
#' flat part and the increasing part. The model is assumed to have a gap between 
#' the left part and the middle part, and a gap between the flat part and the 
#' right part. This algorithm is based on non-parametric MLE.
#' 
#' @param X CFR vector from GeneFishing results from\code{geneFishing()}.
#' @param bait_indices index of bait
#' @param unit searching unit for grid search, defaults to 1 / number of fishing
#' rounds
#' @param left_lb lower bound of left peak of distribution, defaults to minimum
#' of X plus one unit
#' @param left_ub upper bound of the left peak of distribution, defaults to 
#' left_lb + 0.2
#' @param right_ub upper bound of right peak of distribution, defaults to the
#' maximum of X minus one unit
#' @param right_lb lower bound of right peak of distribution, defaults to 
#' right_ub - 0.2
#' @param mu excepted flat part of distribution. defaults to the middle of the
#' left_ub and right_lb
#' @param eps_l the gap between the left part and the middle part.
#' @param eps_r the gap between the right part and the middle part
#' @param topK parameter for finding cutoff
#' @param max_pct pct of genes if clustered with bait 100% we use as cutoff
#' 
#' @return CFR cutoff for labeling fished genes
#' 

getCutoff <- function(X,
                      bait_indices,
                      unit,
                      left_lb = NULL,
                      left_ub = NULL,
                      right_ub = NULL,
                      right_lb = NULL,
                      mu = NULL,
                      eps_l = 0.01,
                      eps_r = 0.01,
                      topK = 1,
                      max_pct = 0.025) {
  
  # first check if there are enough genes that are fished out almost every 
  # time 
  if(sum(X[-bait_indices] >= max(X) - unit) >= 
     (length(X[-bait_indices]) * max_pct)){
    return(max(X) - unit)
  } 
  
  # otherwise we will search for the appropriate cutoff
  if(is.null(left_lb)){
    left_lb <- min(X) + unit
  }
  if(is.null(right_ub)){
    right_ub <- max(X) - unit
  }
  if(is.null(left_ub)){
    left_ub <- left_lb + 0.3
  }
  if(is.null(right_lb)){
    right_lb <- right_ub - 0.1
  }
  if(is.null(mu)){
    mu <- (left_ub + right_lb) / 2
  }

  # check if there are more than 1 value on either side of mu, otherwise it will
  # fail 
  if(length(unique(X[X <= mu])) < 2 | length(unique(X[X > mu])) < 2){
    # in that case we will return a cutoff will all genes that are above mu
    return(max(X) - unit)
  }
  
  n_X <- length(X)
  count_X <- table(sort(X))
  unique_X <- unique(sort(X))
  n_unq_X <- length(unique_X)
  unique_width <- c(min(diff(unique_X)), diff(unique_X))
  
  # g_l, g_r
  left_grenander <-
    fdrtool::grenander(ecdf(X[X <= mu]), "decreasing")
  g_l <- approx(
    x = left_grenander$x.knots,
    y = left_grenander$f.knots,
    xout = unique_X[unique_X <= mu],
    yleft = switch(
      'left',
      'left' = max(left_grenander$f.knots),
      'zero' = 0
    ),
    yright = switch(
      'right',
      'right' = min(left_grenander$f.knots),
      'zero' = 0
    )
  )$y
  g_l <- pmin(1e8, g_l)
  g_l <- g_l / sum(g_l * unique_width[unique_X <= mu])
  
  right_grenander <-
    fdrtool::grenander(ecdf(X[X > mu]), "increasing")
  g_r <- approx(
    x = right_grenander$x.knots,
    y = right_grenander$f.knots,
    xout = unique_X[unique_X > mu],
    yleft = switch(
      'left',
      'left' = min(right_grenander$f.knots),
      'zero' = 0
    ),
    yright = switch(
      'right',
      'right' = max(right_grenander$f.knots),
      'zero' = 0
    )
  )$y
  g_r <- pmin(1e8, g_r)
  g_r <- g_r / sum(g_r * unique_width[unique_X > mu])
  
  
  left_search <- seq(left_lb, min(left_ub, mu), by = unit)
  right_search <- seq(max(right_lb, mu), right_ub, by = unit)
  loglikelihoods <- sapply(left_search, function(c0) {
    logliks <- sapply(right_search, function(c1) {
      N_l <- length(X[X <= c0])
      N_r <- length(X[X > c1])
      N_m <- n_X - N_l - N_r
      
      # estimate the empirical mass
      alpha_l <- N_l / n_X
      alpha_m <- N_m / n_X
      alpha_r <- N_r / n_X
      
      tmp_g_l <-
        g_l[unique_X[unique_X <= mu] <= c0] / 
        sum(g_l[unique_X[unique_X <= mu] <= c0] * 
              unique_width[unique_X <= c0])
      tmp_g_r <-
        g_r[unique_X[unique_X > mu] > c1] / 
        sum(g_r[unique_X[unique_X > mu] > c1] * 
              unique_width[unique_X > c1])
      
      if (min(tmp_g_l) * alpha_l >= eps_l + alpha_m / (c1 - c0) &
          min(tmp_g_r) * alpha_r >= eps_r + alpha_m / (c1 - c0)) {
        loglik <-
          EUILogLikelihood(
            X,
            c0,
            c1,
            alpha_l,
            alpha_r,
            tmp_g_l,
            tmp_g_r,
            N_l = N_l,
            N_r = N_r,
            N_m = N_m,
            n_X = n_X,
            count_X = count_X,
            unique_X = unique_X,
            n_unq_X = n_unq_X
          )
      } else{
        loglik <- -Inf
      }
      return(loglik)
    })
    return(logliks)
  })
  rownames(loglikelihoods) <- right_search
  colnames(loglikelihoods) <- left_search
  
  topK_loglikelihoods <-
    sort(unique(loglikelihoods), decreasing = TRUE)[seq(topK)]
  max_idxs <-
    which(loglikelihoods >= min(topK_loglikelihoods), arr.ind = TRUE)
  est_c0s <-
    sapply(max_idxs[, 2], function(max_idx)
      left_search[max_idx])
  est_c1s <-
    sapply(max_idxs[, 1], function(max_idx)
      right_search[max_idx])
  
  final_N_ls <-
    sapply(est_c0s, function(est_c0)
      length(X[X <= est_c0]))
  final_N_rs <-
    sapply(est_c1s, function(est_c1)
      length(X[X > est_c1]))
  final_N_ms <- n_X - final_N_ls - final_N_rs
  
  
  # estimate the empirical mass
  final_alpha_ls <- final_N_ls / n_X
  final_alpha_ms <- final_N_ms / n_X
  final_alpha_rs <- final_N_rs / n_X
  
  return(max(est_c1s))
}

