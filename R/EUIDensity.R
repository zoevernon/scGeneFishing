#' 
#' EUIDensity
#' 
#' Helper function for \code{getCutoff}.
#' 
#' @param X CFR vector
#' @param c0 from getCutoff
#' @param c1 from getCutoff
#' @param alpha0 from getCutoff
#' @param alpha1 from getCutoff
#' @param g0 from getCutoff
#' @param g1 from getCutoff
#' @param ... additional arguments
#' 
EUIDensity <- function(X, c0, c1, alpha0, alpha1, g0, g1, ...) {
  X_list <- list(...)
  n_X <- X_list$n_X
  count_X <- X_list$count_X
  unique_X <- X_list$unique_X
  n_unq_X <- length(unique_X)
  N_l <- X_list$N_l
  N_r <- X_list$N_r
  N_m <- X_list$N_m
  
  
  density_m <- (1 - alpha0 - alpha1) / (c1 - c0)
  f <- c(alpha0 * g0, rep(density_m, length(
    unique_X[unique_X > c0 & unique_X <= c1])), alpha1 * g1)
  return(f)
}