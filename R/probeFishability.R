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
#' @param parallel indicator of whether the computation should be done in parallel.
#' The default is FALSE.  
#' @param ncore number of cores to do parallel computation on if parallel = TRUE.
#' 
#' @return List with the following components:
#' \itemize{
#'   \item \code{best_bait} - tightest set of genes found.  
#'   \item \code{bait_sets} - list with all tight sets of genes and the corresponding
#'   tightness metric.  
#'   \item \code{bait_info} - data.frame with the number and tightness metric
#'   of each of the discovered bait sets. The \code{bait_index} column corresponds
#'   to the slot in the \code{bait_sets} list (e.g. \code{bait_index} = 1 has
#'   the metadata for \code{bait_sets[[1]]})
#' }
#' 
#' @examples
#' 
#' @import foreach
#' @import ggplot2
#' @export

probeFishability <- function(exp_mat, potential_bait, n_rounds = 100, alpha = 5,
                             min_tightness = 0.5, min_genes = 5,
                             ncores = 2){
  doParallel::registerDoParallel(ncores)
  
  # split the potential bait until we find groups of bait 
  continue_processing <- TRUE
  bait <- list()
  tightness_index <- list()
  cnt <- 1
  round <- 1
  
  # check tightness of whole set
  db_index <- computeAvgDBIndex(potential_bait, exp_mat, 
                                n_rounds = n_rounds, alpha = alpha) %>%
    mean()
  
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
    potential_bait_temp <- potential_bait_temp[vec_lengths >= min_genes]
    
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
  }else {
    bait_final <- list(list(bait = bait[[1]], 
                            tightness = tightness_index[[1]]))
  }
  
  # save a data frame with the tightness and number of genes of each bait 
  # set
  number <- sapply(1:length(bait_final), function(i) length(bait_final[[i]]$bait))
  tightness <- sapply(1:length(bait_final), function(i) bait_final[[i]]$tightness)
  bait_df <- data.frame(bait_index = 1:length(number),
                        bait_length = number,
                        tightness) %>% arrange(tightness)
  
  # reorder original bait list based on tightness
  if(length(bait_final) > 0){
    bait_final <- lapply(1:nrow(bait_df), function(i){
      bait_final[[bait_df$bait_index[i]]]$bait
    })
    
    # fix bait ID (bait now in order from tightest to least tight)
    bait_df$bait_index <- 1:nrow(bait_df)
  }else{
    best_final <- NULL
    bait_df <- NULL
  }
  
  # save output as a class to be used for plotting etc. 
  final_output <-list(best_bait = bait_final[[1]], 
                      bait_sets = bait_final, 
                      bait_info = bait_df,
                      exp_mat = exp_mat)
  
  class(final_output) <- "gene_fishing_probe"
  
  return(final_output)
}

#' @export
print.gene_fishing_probe <- function(x, ...){
  if(length(x$best_bait) > 5){
    cat(paste(length(x$best_bait), "in tightest bait set:\n"))
    cat(sort(x$best_bait)[1:5], ", ...\n\n")
    cat(paste("Found", nrow(x$bait_info) - 1, "additional bait sets"))
  }else if(length(x$best_bait) == 0){
    cat("No bait found, try again with higher min_tightness.\n")
    cat("Note, it is not recommended to use min_tightness > 0.5.")
  }else{
    cat("Tightest bait set:\n")
    cat(sort(x$best_bait), ", ...\n\n")
    cat(paste("Found", nrow(x$bait_info) - 1, "additional bait sets"))
  }
  
}

#' @export
plot.gene_fishing_probe <- function(x, n_random = 50, ...){
  if(length(x$best_bait) == 0){
    cat("No bait found, try again with higher min_tightness.\n")
    cat("Note, it is not recommended to use min_tightness > 0.5.")
  }else{
    all_bait <- sapply(1:nrow(x$bait_info), 
                       function(i) x$bait_sets[[i]]) %>% unlist()
    bait_df <- x$bait_info
    
    # look at plots 
    rand_genes <- sample(setdiff(rownames(x$exp_mat), all_bait), n_random)
    genes <- c(rand_genes, all_bait)
    cor_mat <- cor(t(x$exp_mat[genes, ]), method = "spearman")
    
    eigen_df <- foreach(i = 1:nrow(bait_df), .combine = "rbind") %do% {
      bait <- x$bait_sets[[i]]
      label <- paste0("Bait ", i, ": ", round(bait_df$tightness[i], 2))
      
      # remove any rows that are missing 
      cor_mat_tmp <- cor_mat[c(rand_genes, bait), 
                             c(rand_genes, bait)]
      
      # see if any columns have zero variance 
      cols <- which(is.na(colSds(cor_mat_tmp, na.rm = TRUE)))
      if(length(cols) > 0){
        cor_mat_tmp <- cor_mat_tmp[-cols, -cols]
      }
      
      eigen_space <- getSpectralCoordinates(cor_mat_tmp, 2)$coordinates %>% 
        data.frame() 
      
      
      eigen_space %>% dplyr::mutate(bait = ifelse(rownames(eigen_space) %in% 
                                                    as.character(bait), "bait", "random")) %>%
        dplyr::mutate(gene = rownames(eigen_space[['coordinates']]),
                      label = label)
    }
    
    p <- ggplot(eigen_df) + 
      geom_point(aes(x = eigen.1, y = eigen.2, color = bait), 
                 alpha = 0.5) + 
      theme_bw() + scale_color_discrete(name = "Gene type") +
      facet_wrap(~label) + theme(legend.position = "top")
    
    return(p)
  }
 
}

