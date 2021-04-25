#' 
#' spectralClustering
#' 
#' Internal function to do spectral clustering for \code{geneFishing} and 
#' \code{probeFishability}.
#'

spectralClustering <- function(cor_mat, k){
  # get coordinates for spectral clustering
  rs <- getSpectralCoordinates(cor_mat, k)
  coordinates <- rs[['coordinates']]
  
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
