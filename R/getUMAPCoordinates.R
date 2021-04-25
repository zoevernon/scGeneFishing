#' getUMAPCoordinates
#' 
#' Internal function to compute UMAP coordinates from similarity matrix
#'
#' @param A similarity matrix 
#' @param n_eigen_vectors number of eigen vectors to save for clustering
#' 
#' @return \code{coordinate_matrix} - matrix of UMAP coordinates for each gene.
#' 

getUMAPCoordinates <- function(A, n_eigen_vectors = 0, n_neighbors = 15){ 
  diag(A) <- 0
  A <- abs(A)
  d <- apply(A, 1, sum)
  I <- diag(1, nrow = nrow(A))
  L <- I - diag(1/sqrt(d)) %*% A %*% diag(1/sqrt(d))
  
  # L is positive semi-definite and have n non-negative real-valued eigenvalues
  tmp <- eigen(L)
  eigen_values <- tmp$values[length(tmp$values):1]
  eigen_vectors <- tmp$vectors[, length(tmp$values):1]
  
  # automaticaly choose number of eigen-vectors
  q1 <- quantile(eigen_values)[2]
  q3 <- quantile(eigen_values)[4]
  low <- q1 - 3 * (q3 - q1)
  up <- q1 + 3 * (q3 - q1)
  
  # pick out eigen-values that are relatively small
  k <- sum(eigen_values < low) 
  
  if(n_eigen_vectors > 1){
    k <- n_eigen_vectors
  }
  if(k == 1){
    rs <- list(eigen_values = eigen_values, n_eigen_vectors = 1)
  }else{
    # pick out the eigen-vectors associated with the k smallest non-zero 
    # eigen-values
    coordinate_matrix <- eigen_vectors[, 2:(k + 1)]
    rownames(coordinate_matrix) <- rownames(A)
    colnames(coordinate_matrix) <- paste('eigen', 1:(k), sep = "-")
    
    rs <- list(coordinates = coordinate_matrix,
               eigen_values = eigen_values,
               n_eigen_vectors = k)
  }
  
  return(umap(rs[['coordinates']], n_neighbors = n_neighbors)$layout)
}
