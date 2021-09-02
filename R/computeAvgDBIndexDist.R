#' 
#' computeAvgDBIndexDist
#' 
#' Internal function to compute tightness of gene set when probing for fishability
#' 
#' @param bait_genes bait or potential bait
#' @param exp_mat expression matrix
#' @param method similarity metric
#' @param k clusters
#' @param alpha number of genes to sample
#' @param n_rounds number of rounds of sampling
#' 
computeAvgDBIndexDist <- function(bait_genes, exp_mat, method, 
                                    k = 2, alpha = 5, n_rounds = 50){
  
  # make sure the bait genes are represented 
  bait_genes <- bait_genes[bait_genes %in% rownames(exp_mat)] %>%
    as.character()
  
  # get all other genes
  all_genes  <- rownames(exp_mat)       
  pool_genes <- setdiff(all_genes, bait_genes)
  
  # and across different samples of random genes
  bait_tightness <- 
    foreach(i = 1:n_rounds, .combine = c) %dopar% {
      # sample random genes 
      rand_genes <- sample(pool_genes, 
                           size = length(bait_genes) * alpha)
      
      # get correlation matrix for those genes
      dist_subset <- dist(exp_mat[c(bait_genes, rand_genes), ] %>% as.matrix(),
                          method = method) %>% as.matrix()
      
      # get coordinates for spectral clustering
      rs <- getSpectralCoordinates(1 / (dist_subset + 1), k)
      
      coordinates <- rs[['coordinates']]
      
      # compute modified DB index, first compute centroids of clusters
      centroids <- rbind(Rfast::colMedians(coordinates[bait_genes, ]),
                         Rfast::colMedians(coordinates[rand_genes, ]))
      
      
      avg_bait_dist <- Rfast::Dist(coordinates[bait_genes, ]) %>%
        as.matrix() %>%
        gdata::upperTriangle() %>%
        median()
      
      # compute index
      avg_bait_dist / (Rfast::Dist(centroids, vector = TRUE) %>% as.numeric())
    }
  
  
  return(bait_tightness)
}

