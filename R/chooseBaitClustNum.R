#' 
#' chooseBaitClustNum
#' 
#' Internal function to choose number of clusters in performSplit
#' 

chooseBaitClustNum <- function(cor_mat, potential_bait, max_clusters = 5) {
  # check if max_clusters is greater than the number of genes in the bait
  if(max_clusters > (length(potential_bait) - 1)) {
    max_clusters <- length(potential_bait) - 1
  }
  
  # compute CH index for each potential cluster number
  CH_indices <- foreach(j = 2:max_clusters, .combine = 'c') %do%{
    # do clustering 
    spectral <- spectralClustering(cor_mat, j)
    
    # compute index 
    computeCH(spectral$kmeans, nrow(spectral$coordinates), j)
  }
  
  # make plot to look at (note I should make this into a method)
  plot <- ggplot() + 
    geom_point(aes(x = 2:max_clusters, y = CH_indices)) + 
    geom_line(aes(x = 2:max_clusters, y = CH_indices)) + 
    theme_bw() + labs(x = "Number of clusters", y = "CH")
  
  return(list(cluster_number = which.max(CH_indices) + 1, plot = plot))
}