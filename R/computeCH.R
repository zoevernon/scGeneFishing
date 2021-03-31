#' 
#' geneFishing
#' 
#' Internal function for probing fishability
#' 
computeCH <- function(kmeans_obj, n_genes, n_clusters){
  kmeans_obj$betweenss / kmeans_obj$tot.withinss * 
    ((n_genes - n_clusters) / (n_clusters - 1))
}
