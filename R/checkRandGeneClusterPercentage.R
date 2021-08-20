#' 
#' checkRandGeneClusterPercentage
#' 
#' Function to check that not too many random genes are clustering with the 
#' discovered bait.  Will return a warning if that is the case.
#' 

checkRandGeneClusterPercentage <- function(results, i, n_rounds, alpha, k = 2){
  # make sure the bait genes are represented 
  bait_genes <- results$bait_sets[[i]]
  
  # get all other genes
  exp_mat <- results$X
  all_genes  <- rownames(exp_mat)       
  pool_genes <- setdiff(all_genes, bait_genes)
  
  # and across different samples of random genes
  avg_perct <- 
    foreach(i = 1:n_rounds, .combine = c) %dopar% {
      # sample random genes 
      rand_genes <- sample(pool_genes, 
                           size = length(bait_genes) * alpha)
      
      # get correlation matrix for those genes
      if(results$method %in% c("euclidean",  "maximum", "manhattan", 
                               "canberra", "binary", "minkowski")){
        cor_mat <- dist(exp_mat[c(bait_genes, rand_genes), ] %>% as.matrix(),
                        method = results$method) %>% as.matrix()
        cor_mat <- 1 / (cor_mat + 1)
      }else if(results$method == "cosine"){
        cor_mat <- cosineSimMatrix(exp_mat[c(bait_genes, rand_genes), ] %>% 
                                     as.matrix())
      }else{
        cor_mat <- cor(t(exp_mat[c(bait_genes, rand_genes), ] %>% 
                           as.matrix()), method = results$method)
      }
      
      # get coordinates for spectral clustering
      if(results$type == "UMAP"){
        coordinates <- getUMAPCoordinates(cor_mat, k)
      }else{
        rs <- getSpectralCoordinates(cor_mat, k)
        coordinates <- rs[['coordinates']]
      }

      # do k-means on the coordinates from the spectral clustering
      kmeans_obj <- 
        LICORS::kmeanspp(data = coordinates, k = k)
      
      # extract clusters
      cluster_vec <- kmeans_obj$cluster
      names(cluster_vec) <- rownames(coordinates)
      
      # figure out which genes, if any, cluster with the bait
      cluster_freq_df <- table(cluster_vec[bait_genes]) %>% 
        as.data.frame
      cluster_freq_df <- cluster_freq_df[order(cluster_freq_df$Freq, 
                                               decreasing = TRUE),]
      value <- cluster_freq_df$Var1[1] %>% 
        as.character %>% 
        as.integer
      
      # all genes in cluster with majority of bait
      f1 <- cluster_vec == value
      
      sum(f1) / length(f1)
    } 
  
  return(avg_perct %>% mean())
}
