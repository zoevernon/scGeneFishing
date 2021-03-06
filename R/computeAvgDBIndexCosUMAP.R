#' 
#' computeAvgDBIndexCosUMAP
#' 
#' Internal function to compute tightness of gene set when probing for fishability
#' 
#' @param bait_genes bait or potential bait
#' @param exp_mat expression matrix
#' @param method similarity metric
#' @param k clusters
#' @param alpha number of genes to sample
#' @param n_rounds number of rounds of sampling
#' @param n_neighbors number of neighbors for UMAP
#' 
computeAvgDBIndexCosUMAP <- function(bait_genes, exp_mat, method, 
                                  k = 2, alpha = 5, n_rounds = 50, 
                                  n_neighbors = 15){
  # make sure the bait genes are represented 
  bait_genes <- bait_genes[bait_genes %in% rownames(exp_mat)] %>%
    as.character()
  
  # get all other genes
  all_genes  <- rownames(exp_mat)       
  pool_genes <- setdiff(all_genes, bait_genes)
  
  # and across different samples of random genes
  bait_tightness <- 
    foreach(i = 1:n_rounds, .combine = c) %dopar% {
      result <- tryCatch({
        # sample random genes 
        rand_genes <- sample(pool_genes, 
                             size = length(bait_genes) * alpha)
        
        # get correlation matrix for those genes
        cor_subset <- cosineSimMatrix(exp_mat[c(bait_genes, rand_genes), ] %>% 
                                        as.matrix())
        
        # get coordinates for spectral clustering
        coordinates <- getUMAPCoordinates(cor_subset, k, n_neighbors)
        
        # compute modified DB index, first compute centroids of clusters
        centroids <- rbind(Rfast::colMedians(coordinates[bait_genes, ]),
                           Rfast::colMedians(coordinates[rand_genes, ]))
        
        
        avg_bait_dist <- Rfast::Dist(coordinates[bait_genes, ]) %>%
          as.matrix() %>%
          gdata::upperTriangle() %>%
          median()
        
        # compute index
        avg_bait_dist / 
          sqrt((Rfast::Dist(centroids, vector = TRUE) %>% as.numeric()))
      },
      error = function(err) {
        NA
      })
      result
    }
  
  
  return(bait_tightness[!is.na(bait_tightness)])
}

