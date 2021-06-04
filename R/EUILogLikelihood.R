#' 
#' EUILogLikelihood
#' 
#' Helper function for \code{getCutoff}.
#' 
EUILogLikelihood <- function(X, c0, c1, alpha0, alpha1, g0, g1, ...) {
  X_list <- list(...)
  n_X <- X_list$n_X
  count_X <- X_list$count_X
  unique_X <- X_list$unique_X
  n_unq_X <- length(unique_X)
  N_l <- X_list$N_l
  N_r <- X_list$N_r
  N_m <- X_list$N_m
  
  f <- EUIDensity(
    X,
    c0,
    c1,
    alpha0,
    alpha1,
    g0,
    g1,
    N_l = N_l,
    N_r = N_r,
    N_m = N_m,
    n_X = n_X,
    count_X = count_X,
    unique_X = unique_X,
    n_unq_X = n_unq_X
  )
  
  return(sum(log(f) * count_X))
}
