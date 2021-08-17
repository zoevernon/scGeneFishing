#' cosineSimMatrix
#' 
#' Internal function to get cosine similarity
#'
#' @param mat matrix to compute similarity
#' 
#' @return similarity matrix between rows of \code{mat}
#' 
cosineSimMatrix <- function(mat){
  sim <- mat / sqrt(rowSums(mat * mat))
  return(sim %*% t(sim))
}
