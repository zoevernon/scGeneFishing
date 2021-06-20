#' 
#' geneFishingUMAP
#' 
#' Function to do UMAP version of geneFishing
#' 
#' @param exp_mat matrix of gene expression where the rows are genes and the 
#' columns are samples (cells or individuals). 
#' @param bait_genes genes to use as bait to fish out similar genes from the 
#' \code{exp_mat}.  These genes should be a subset of the rownames of \code{exp_mat}. 
#' @param alpha controls number of random genes that are sampled in each round 
#' of fishing.  It will sample \code{alpha} times the number of genes in 
#' \code{bait_genes}.  The default is 5.  The stronger the bait separates from 
#' samples of random genes the larger \code{alpha} the algorithm can handle.
#' @param fishing_rounds number of rounds of fishing to do.  The default is 1,000.
#' @param k number of clusters to use in geneFishing.  The default is 2, however
#' if you fished out too many genes with \code{k = 2} it can help to increase to 
#' \code{k = 3}. 
#' 
geneFishingUMAP <- function(exp_mat, bait_genes, pool_genes, 
                                alpha, fishing_rounds, k, method){
  
  # create lists of genes to be fished together for each of the rounds of gene
  # fishing
  tmp <- foreach(i = 1:fishing_rounds) %dopar%{
    # shuffle the genes to do sampling
    shuffle_pool_genes <- sample(pool_genes, length(pool_genes))
    sub_pool_size <- length(bait_genes) * alpha
    
    # split the shuffled gene vector into lists so that each gene is in one
    # of the lists 
    l <- foreach(j = seq(from = 1, to = length(shuffle_pool_genes), 
                         by = sub_pool_size)) %do% {
                           up <- min(j + sub_pool_size - 1, length(shuffle_pool_genes))
                           shuffle_pool_genes[j:up]
                         }  
    
    s <- l[[length(l)]]
    d <- l[[length(l)-1]]
    d <- c(s,d)
    l[[length(l)-1]] <- d
    l <- l[1:(length(l) - 1)]
    l
  }
  
  # flatten the list, so that we just have one long list
  sub_pool_list <- unlist(tmp, recursive = FALSE)
  
  fish_vec <- foreach(sub_pool = sub_pool_list, .combine = 'c') %dopar% {
    result <- tryCatch({
      bait_indices <- which(rownames(exp_mat) %in% bait_genes)
      sub_indices <- which(rownames(exp_mat) %in% sub_pool)
      cor_mat_bait_and_pool <- cor(t(exp_mat[c(bait_indices, 
                                               sub_indices), ]) %>%
                                     as.matrix(), 
                                   method = method)
      
      # get coordinates to do the spectral clustering on the coexpression matrix
      # above 
      coordinates <- getUMAPCoordinates(cor_mat_bait_and_pool, k)
      
      # do k-means on the coordinates from the spectral clustering
      kmeans_obj <- LICORS::kmeanspp(data = coordinates, k = k)
      
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
      
      # all bait genes 
      f2 <- names(cluster_vec) %in% bait_genes
      
      # print out all the genes that are in the cluster with majority of the 
      # bait including the bait genes
      names(cluster_vec)[f1]
    },
    error = function(err) {
      paste0("ERROR_", c(sub_pool, bait_genes))
    })
    result
  }
  
  # update the frequency that genes are fished
  fish_freq_df <- table(fish_vec) %>% as.data.frame 
  colnames(fish_freq_df) <- c('fish_id', 'fish_freq')
  fish_freq_df <-  fish_freq_df %>% filter(fish_id != "NOTHING")
  
  # add column with the error and remove error tag from fish id
  fish_freq_df$ERROR <- ifelse(stringr::str_detect(
    fish_freq_df$fish_id, "ERROR"), 1, 0)
  fish_freq_df$fish_id <- stringr::str_remove(fish_freq_df$fish_id, "ERROR_")
  
  # combine genes were there was an error and change denominator to match
  # the error 
  fish_freq_df <- fish_freq_df %>% 
    group_by(fish_id) %>%
    mutate(fish_freq_error = sum(ERROR)) %>%
    filter(ERROR == 0) %>%
    dplyr::select(-ERROR) %>%
    ungroup() %>%
    mutate(denom = ifelse(fish_id %in% bait_genes, 
                          length(sub_pool_list) - fish_freq_error,
                          fishing_rounds - fish_freq_error)) %>% 
    mutate(fish_freq = fish_freq / denom) %>%
    arrange(-fish_freq) %>%
    mutate(fish_id = as.character(fish_id)) %>%
    dplyr::select(-denom, -fish_freq_error)
  
  # for genes that were not fished at all save them has having been fished 
  # out none of the times.  
  if(length(setdiff(pool_genes, fish_freq_df$fish_id)) > 0){
    df <- data.frame(fish_id = setdiff(pool_genes, fish_freq_df$fish_id),
                     fish_freq = 0)
    fish_freq_df <- rbind(fish_freq_df, df)
  }
  
  # compute the p-value for the gene fishing and adjust
  p_hat <- mean(fish_freq_df$fish_freq)
  fish_freq_df$p_value <- pbinom(fish_freq_df$fish_freq * fishing_rounds, 
                                 fishing_rounds, p_hat, lower.tail = FALSE)
  fish_freq_df$adj_p_value <- p.adjust(fish_freq_df$p_value, 
                                       method = 'bonferroni')
  fish_freq_df$bait <- ifelse(fish_freq_df$fish_id %in% bait_genes, "yes", "no")
  colnames(fish_freq_df)  <- c('gene_id',
                               'CFR', 
                               'p_value', 
                               'adj_p_value',
                               'bait')
  return(fish_freq_df)
}
