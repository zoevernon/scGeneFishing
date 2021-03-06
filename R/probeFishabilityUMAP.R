#' 
#' probeFishabilityUMAP
#' 
#' Function to find a set of genes that can be used as bait in geneFishing from 
#' a larger set of potential bait genes. 
#' 
#' @param exp_mat matrix of gene expression where the rows are genes and the 
#' columns are samples (cells or individuals). 
#' @param potential_bait set of genes for which to check if there is a 
#' subset that can be used as bait.  These genes should be a subset of the 
#' rownames of \code{exp_mat}. 
#' @param n_rounds number of random samples to be used in assessing fishability
#' of a set of genes.
#' @param alpha controls number of random genes that are sampled in each round 
#' of fishing.  It will sample \code{alpha} times the number of genes in 
#' \code{bait_genes}. 
#' @param min_tightness parameter that controls how tight the set of genes needs
#' to be to be considered as bait.  Tightness is measured computing the ratio of
#' the median distance between all potential bait genes to the distance between the
#' medoid of the potential bait genes and a set of randomly sampled genes. The
#' average of these tightness measures are taken over \code{n_rounds} random samples
#' of genes.  
#' @param min_genes minimum size that the bait set can be. 
#' @param n_neighbors number of neighbors for UMAP
#' @param method similarity metric
#'  
probeFishabilityUMAP <- function(exp_mat, potential_bait, n_rounds, alpha, 
                                 min_tightness, min_genes, n_neighbors,
                                 method){
  
  # split the potential bait until we find groups of bait 
  continue_processing <- TRUE
  bait <- list()
  tightness_index <- list()
  cnt <- 1
  round <- 1
  
  # check tightness of whole set
  if(method == "cosine"){
    db_index <- computeAvgDBIndexCosUMAP(potential_bait, exp_mat, 
                                         n_rounds = n_rounds, 
                                         alpha = alpha, 
                                         n_neighbors = n_neighbors,
                                         method = method) %>%
      mean(na.rm = TRUE)
  }else{
    db_index <- computeAvgDBIndexUMAP(potential_bait, exp_mat, 
                                      n_rounds = n_rounds, 
                                      alpha = alpha, 
                                      n_neighbors = n_neighbors,
                                      method = method) %>%
      mean(na.rm = TRUE)
  }
  
  
  if(db_index < min_tightness) {
    continue_processing <- FALSE
    bait <- list(potential_bait)
    tightness_index <- list(db_index)
  }
  
  # list of potential bait to update and loop through
  potential_bait_list <- list(potential_bait)
  while(continue_processing) {
    potential_bait_temp <- list()
    cnt2 <- 1
    for(i in seq(length(potential_bait_list))) {
      # split each element 
      # n_neighbors_tmp <- ifelse(
      #   n_neighbors > round(length(potential_bait_temp) / 1.5), 
      #   round(length(potential_bait_temp) / 1.5),
      #   n_neighbors)
      # n_neighbors_tmp <- ifelse(n_neighbors_tmp < 2, 2, n_neighbors_tmp)
      
      split_data <- performSplitUMAP(potential_bait_list[[i]], 
                                     exp_mat, 
                                     n_rounds = n_rounds, 
                                     round = round, 
                                     alpha = alpha, 
                                     n_neighbors = n_neighbors,
                                     min_genes = min_genes,
                                     method = method)
      
      # check if any of these splits result in a DB index of less than the 
      # cutoff
      subsetted <- split_data %>% filter(db_index_vec <= min_tightness)
      if(nrow(subsetted) > 0) {
        # save the genes in the bait list
        clusters <- unique(subsetted$cluster)
        for(k in clusters) {
          bait_genes <- subsetted$gene[subsetted$cluster == k] %>% 
            as.character()
          
          # only keep the bait if there are more than 5 genes (default), will
          # change for cross validation procedure.  
          if(length(bait_genes) >= min_genes){
            bait[[cnt]] <- bait_genes
            tightness_index[[cnt]] <- 
              subsetted$db_index_vec[subsetted$cluster == k] %>%
              unique()
            cnt <- cnt + 1
          }
        }
      }
      
      # update the potential bait list 
      subsetted <- split_data %>%
        filter(!(gene %in% unlist(bait)))
      for(k in unique(subsetted$cluster)) {
        potential_bait_temp[[cnt2]] <- 
          subsetted$gene[subsetted$cluster == k] %>% as.character()
        cnt2 <- cnt2 + 1
      }
    }
    
    # remove any of the elements that have less than min_genes + 2
    vec_lengths <- sapply(potential_bait_temp, length)
    potential_bait_temp <- potential_bait_temp[vec_lengths >= (min_genes + 2)]
    
    # reset the potential bait list
    potential_bait_list <- potential_bait_temp
    
    # update overall round (used for the assessing number of clusters)
    round <- round + 1
    
    # potential_bait_list
    # bait
    # stop processing if there are no genes left 
    if(length(potential_bait_list) == 0) {
      continue_processing = FALSE
    }
  }
  
  # make bait into a more usable list 
  if(length(bait) > 1){
    bait_final <- lapply(1:length(bait), function(i){
      list(bait = bait[[i]], tightness = tightness_index[[i]])
    })
  }else if(length(bait) == 1){
    bait_final <- list(list(bait = bait[[1]], 
                            tightness = tightness_index[[1]]))
  }else{
    bait_final <- bait
  }
  
  # reorder original bait list based on tightness
  if(length(bait_final) > 0){
    # save a data frame with the tightness and number of genes of each bait 
    # set
    number <- sapply(1:length(bait_final), function(i) length(bait_final[[i]]$bait))
    tightness <- sapply(1:length(bait_final), function(i) bait_final[[i]]$tightness)
    bait_df <- data.frame(bait_index = 1:length(number),
                          bait_length = number,
                          tightness) %>% arrange(tightness)
    
    bait_final <- lapply(1:nrow(bait_df), function(i){
      bait_final[[bait_df$bait_index[i]]]$bait
    })
    
    # fix bait ID (bait now in order from tightest to least tight)
    bait_df$bait_index <- 1:nrow(bait_df)
  }else{
    bait_final <- list(NULL)
    bait_df <- NULL
    warning("No bait discovered.")
  }
  
 
  # save output as a class to be used for plotting etc. 
  final_output <-list(best_bait = bait_final[[1]], 
                      bait_sets = bait_final, 
                      bait_info = bait_df,
                      X = exp_mat,
                      method = method,
                      type = "UMAP")
  
  class(final_output) <- "gene_fishing_probe"
  
  return(final_output)
}
