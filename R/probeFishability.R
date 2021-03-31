#' 
#' probeFishability
#' 
#' Function to find a set of genes that can be used as bait in geneFishing from 
#' a larger set of potential bait genes. 
#' 
#' @param exp_mat matrix of gene expression where the rows are genes and the 
#' columns are samples (cells or individuals). 
#' @param potential_bait_genes set of genes for which to check if there is a 
#' subset that can be used as bait.  These genes should be a subset of the 
#' rownames of \code{exp_mat}. 
#' @param n_rounds number of random samples to be used in assessing fishability
#' of a set of genes.  The default is 100.  
#' @param alpha controls number of random genes that are sampled in each round 
#' of fishing.  It will sample \code{alpha} times the number of genes in 
#' \code{bait_genes}.  The default is 5. 
#' @param min_tightness parameter that controls how tight the set of genes needs
#' to be to be considered as bait.  Tightness is measured computing the ratio of
#' the median distance between all potential bait genes to the distance between the
#' medoid of the potential bait genes and a set of randomly sampled genes. The
#' average of these tightness measures are taken over \code{n_rounds} random samples
#' of genes.  
#' @param min_genes minimum size that the bait set can be.  The default is set to
#' 5.  
#' 
#' @return List with the following components:
#' \itemize{
#'   \item \code{best_bait} - tightest set of genes found.  
#'   \item \code{bait_sets} - list with all tight sets of genes and the corresponding
#'   tightness metric.  
#'   \item \code{bait_info_df} - data.frame with the number and tightness metric
#'   of each of the discovered bait sets. 
#' }
#' 
#' @examples
#' 
#' @import foreach
#' @import doParallel
#' 

probeFishability <- function(exp_mat, potential_bait, n_rounds = 100, alpha = 5,
                             min_tightness = 0.5, min_genes = 5){
  
  # split the potential bait until we find groups of bait 
  continue_processing <- TRUE
  bait <- list()
  tightness_index <- list()
  cnt <- 1
  round <- 1
  
  # list of potential bait to update and loop through
  potential_bait_list <- list(potential_bait)
  while(continue_processing) {
    potential_bait_temp <- list()
    cnt2 <- 1
    for(i in seq(length(potential_bait_list))) {
      # split each element 
      split_data <- performSplit(potential_bait_list[[i]], 
                                 exp_mat, 
                                 n_rounds = n_rounds, 
                                 round = round, 
                                 alpha = alpha)
      
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
    
    # remove any of the elements that have less than 5 genes 
    vec_lengths <- sapply(potential_bait_temp, length)
    potential_bait_temp <- potential_bait_temp[vec_lengths >= 10]
    
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
  if(length(bait) > 0){
    bait_final <- lapply(1:length(bait), function(i){
      list(bait = bait[[i]], tightness = tightness_index[[i]])
    })
  }else {
    bait_final <- bait
  }
  
  return(bait_final)
}
