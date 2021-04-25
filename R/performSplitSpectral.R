#' 
#' performSplitSpectral
#' 
#' Internal function for probing fishability of a set of genes to split up 
#' potential bait set by doing clustering with k = 2.
#' 
performSplitSpectral <- function(potential_bait, exp_mat, n_rounds, round, alpha) {
  cor_mat <- cor(t(exp_mat[potential_bait, ]), method = "spearman")
  # cluster the potential bait
  num_clust <- 2
  
  # get the eigen_space for correct number of clusters
  eigen_space <- spectralClustering(
    cor_mat, 
    k = num_clust
  )$coordinates 
  
  # make a list of genes in each cluster 
  genes_in_clust <- lapply(1:num_clust, function(k){
    eigen_space$gene[eigen_space$cluster == k]
  })
  
  # check if any of the clusters only have one gene and remove
  genes_in_clust <- genes_in_clust[sapply(genes_in_clust, length) > 1]
  
  # compute DB index for each cluster
  db_index <- sapply(1:length(genes_in_clust), function(k){
    computeAvgDBIndexSpectral(genes_in_clust[[k]], exp_mat, n_rounds = n_rounds, 
                              alpha = alpha) %>%
      mean()
  })    
  
  # return data.frame that has the DB index, cluster and genes in cluster 
  db_index_vec <- sapply(eigen_space$cluster, function(i){
    db_index[i]
  })
  
  return(cbind(eigen_space, db_index_vec))
}
