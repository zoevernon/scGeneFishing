#' 
#' performSplitUMAP
#' 
#' Internal function for probing fishability of a set of genes to split up 
#' potential bait set by doing clustering with k = 2.
#' 
#' @param potential_bait potential bait genes
#' @param exp_mat expression matrix
#' @param n_rounds number of rounds of sampling
#' @param round round index
#' @param alpha number of genes to sample
#' @param n_neighbors number of neighbors for UMAP
#' @param min_genes minimum number of genes
#' @param method similarity method
#' 
performSplitUMAP <- function(potential_bait, exp_mat, n_rounds, round, alpha,
                                 n_neighbors, min_genes, method) {
  if(method == "cosine"){
    cor_mat <- cosineSimMatrix(exp_mat[potential_bait, ] %>% as.matrix())
  }else{
    cor_mat <- cor(t(exp_mat[potential_bait, ]) %>% as.matrix(),
                   method = method)
  }
  
  # cluster the potential bait
  num_clust <- 2
  
  # get the eigen_space for correct number of clusters
  if(nrow(cor_mat) > 15){
    eigen_space <- UMAPClustering(
      cor_mat, 
      k = num_clust,
      n_neighbors = n_neighbors
    )$coordinates
  }else{
    eigen_space <- spectralClustering(
      cor_mat, 
      k = num_clust
    )$coordinates
  }
  
  # make a list of genes in each cluster 
  genes_in_clust <- lapply(1:num_clust, function(k){
    eigen_space$gene[eigen_space$cluster == k]
  })
  
  # compute DB index for each cluster
  db_index <- sapply(1:length(genes_in_clust), function(k){
    if(length(genes_in_clust[[k]]) < min_genes){
      999
    }else{
      if(method == "cosine"){
        tmp <- computeAvgDBIndexCosUMAP(genes_in_clust[[k]], exp_mat, 
                                        n_rounds = n_rounds, 
                                        alpha = alpha,
                                        method = method) %>%
          mean(na.rm = TRUE)
      }else{
        tmp <- computeAvgDBIndexUMAP(genes_in_clust[[k]], exp_mat, 
                                     n_rounds = n_rounds, 
                                     alpha = alpha,
                                     method = method) %>%
          mean(na.rm = TRUE)
      }
      tmp
    }
  })    
  
  # return data.frame that has the DB index, cluster and genes in cluster 
  db_index_vec <- sapply(eigen_space$cluster, function(i){
    db_index[i]
  })

  
  return(cbind(eigen_space, db_index_vec))
}
