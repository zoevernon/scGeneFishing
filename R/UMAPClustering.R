#' 
#' UMAPClustering
#' 
#' Internal function to do UMAP clustering for \code{geneFishing} and 
#' \code{probeFishability}.
#' 
#' @param cor_mat correlation / similarity matrix
#' @param k number of eigen vectors
#' @param n_neighbors number of neighbors for UMAP
#'

UMAPClustering <- function(cor_mat, k, n_neighbors){
  # get coordinates for spectral clustering
  coordinates <- getUMAPCoordinates(cor_mat, k, n_neighbors)
  
  # do clustering 
  kmeans_obj <-  LICORS::kmeanspp(data = coordinates, k = k)
  
  # merge coordinates with clustering 
  coordinates <- coordinates %>% data.frame 
  coordinates$gene <- rownames(coordinates)
  rownames(coordinates) <- NULL
  coordinates <- coordinates %>%
    merge(data.frame(gene = names(kmeans_obj$cluster), 
                     cluster = kmeans_obj$cluster, row.names = NULL), 
          by = "gene")
  
  
  # return kmeans object and coordinates
  return(list(coordinates = coordinates, kmeans = kmeans_obj))
}
