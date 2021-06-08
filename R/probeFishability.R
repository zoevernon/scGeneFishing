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
#' @param umap indicator of whether to do UMAP on spectral coordinates for 
#' clustering.  Default is TRUE, which is recommended for single cell data, despite
#' a reduction in speed.
#' @param ncores number of cores to do parallel computation on if parallel = TRUE.
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
#' @import Matrix
#' @export

probeFishability <- function(X, potential_bait, n_rounds = 100, alpha = 5,
                             min_tightness = 0.5, min_genes = 5, n_neighbors = 15,
                             umap = TRUE, ncores = 2){
  doParallel::registerDoParallel(ncores)
  
  # assertion to make sure the entry is a matrix, sparse matrix or sce
  assertthat::assert_that(
    is.matrix(X) | any(class(X) %in% c("SingleCellExperiment", "dgCMatrix", 
                                       "dgTMatrix", "dgRMatrix")), 
    msg = paste0("X must be a matrix, or in one of the following classes:\n", 
                 "SingleCellExperiment, dgCMatrix, dgTMatrix, dgRMatrix"))
  
  # check what genes overlap with rownames of expression matrix
  potential_bait_orig <- unique(potential_bait)
  potential_bait <- intersect(rownames(X), potential_bait)
  
  # make sure there are some genes in the potential bait
  assertthat::assert_that(
    length(potential_bait) > 0, 
    msg = "No intersection of potential_bait and row names of expression matrix")
  
  # make sure there are at least min_genes in potential bait
  assertthat::assert_that(
    length(potential_bait) >= min_genes, 
    msg = paste0("The interesection of potential_bait and row names of ",
                 "expression matrix is less than min_genes.", 
                 "\nConsider lower min_genes or using new potential bait."))
  
  # return a warning that some genes were removed as they are not in the row-
  # names of the expression matrix
  if(length(potential_bait) < length(potential_bait_orig)) {
    warning(paste0("Inputted potential bait had ", length(potential_bait_orig),
                   " genes, only ", length(potential_bait), 
                   " are in the row names of the inputted expression matrix."))
  }
    
  if(umap){
    if(is.matrix(X)){
      results <- mat_probeFishabilityUMAP(X, 
                                          potential_bait,
                                          n_rounds = n_rounds, 
                                          alpha = alpha, 
                                          min_tightness = min_tightness, 
                                          min_genes = min_genes,
                                          n_neighbors = n_neighbors)
    }else if(class(X) == "SingleCellExperiment"){
      assertthat::assert_that(
        any(assayNames(X) == "logcounts"), 
        msg = paste0("There must be a logcounts assay in SCE."))
      
      results <- sce_probeFishabilityUMAP(X, 
                                          potential_bait,
                                          n_rounds = n_rounds, 
                                          alpha = alpha, 
                                          min_tightness = min_tightness, 
                                          min_genes = min_genes,
                                          n_neighbors = n_neighbors)
    }else if(class(X) %in% c("dgCMatrix", "dgTMatrix", "dgRMatrix")) {
      results <- mat_probeFishabilityUMAP(X, 
                                          potential_bait,
                                          n_rounds = n_rounds, 
                                          alpha = alpha, 
                                          min_tightness = min_tightness, 
                                          min_genes = min_genes,
                                          n_neighbors = n_neighbors)
    }
    
  }else{
    if(is.matrix(X)){
      results <- mat_probeFishabilitySpectral(X, 
                                              potential_bait,
                                              n_rounds = n_rounds, 
                                              alpha = alpha, 
                                              min_tightness = min_tightness, 
                                              min_genes = min_genes)
    }else if(class(X) == "SingleCellExperiment"){
      assertthat::assert_that(
        any(assayNames(X) == "logcounts"), 
        msg = paste0("There must be a logcounts assay in SCE."))
      results <- sce_probeFishabilitySpectral(X, 
                                              potential_bait,
                                              n_rounds = n_rounds, 
                                              alpha = alpha, 
                                              min_tightness = min_tightness, 
                                              min_genes = min_genes)
    } else if(class(X) %in% c("dgCMatrix", "dgTMatrix", "dgRMatrix")) {
      results <- mat_probeFishabilitySpectral(X, 
                                              potential_bait,
                                              n_rounds = n_rounds, 
                                              alpha = alpha, 
                                              min_tightness = min_tightness, 
                                              min_genes = min_genes)
    }
  }
  
  return(results)
}

#' @export
print.gene_fishing_probe_spectral <- function(x, ...){
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
plot.gene_fishing_probe_spectral <- function(x, n_random = 50, ...){
  if(length(x$best_bait) == 0){
    cat("No bait found, try again with higher min_tightness.\n")
    cat("Note, it is not recommended to use min_tightness > 0.5.")
  }else{
    all_bait <- sapply(1:nrow(x$bait_info), 
                       function(i) x$bait_sets[[i]]) %>% unlist()
    bait_df <- x$bait_info
    
    # look at plots 
    rand_genes <- sample(setdiff(rownames(x$X), all_bait), n_random)
    genes <- c(rand_genes, all_bait)
    if(is.matrix(x$X)){
      cor_mat <- cor(t(x$X[genes, ]), method = "spearman")
    } else if(class(x$X) == "SingleCellExperiment"){
      cor_mat <- cor(t(logcounts(x$X)[genes, ])  %>% as.matrix(),
                     method = "spearman")
    } else if(class(x$X) %in% c("dgCMatrix", "dgTMatrix", "dgRMatrix")) {
      cor_mat <- cor(t(x$X[genes, ]) %>% as.matrix(), method = "spearman")
    }
    
    eigen_df <- foreach(i = 1:nrow(bait_df), .combine = "rbind") %do% {
      bait <- x$bait_sets[[i]]
      label <- paste0("Bait ", i, ": ", round(bait_df$tightness[i], 2))
      
      # remove any rows that are missing 
      cor_mat_tmp <- cor_mat[c(rand_genes, bait), 
                             c(rand_genes, bait)]
      
      # see if any columns have zero variance 
      cols <- which(is.na(matrixStats::colSds(cor_mat_tmp, na.rm = TRUE)))
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


#' @export
print.gene_fishing_probe_umap <- function(x, ...){
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
plot.gene_fishing_probe_umap <- function(x, n_random = 50, ...){
  if(length(x$best_bait) == 0){
    cat("No bait found, try again with higher min_tightness.\n")
    cat("Note, it is not recommended to use min_tightness > 0.5.")
  }else{
    all_bait <- sapply(1:nrow(x$bait_info), 
                       function(i) x$bait_sets[[i]]) %>% unlist()
    bait_df <- x$bait_info
    
    # look at plots 
    rand_genes <- sample(setdiff(rownames(x$X), all_bait), n_random)
    genes <- c(rand_genes, all_bait)
    if(is.matrix(x$X)){
      cor_mat <- cor(t(x$X[genes, ]), method = "spearman")
    } else if(class(x$X) == "SingleCellExperiment"){
      cor_mat <- cor(t(logcounts(x$X)[genes, ]) %>% as.matrix(),
                     method = "spearman")
    } else if(class(x$X) %in% c("dgCMatrix", "dgTMatrix", "dgRMatrix")) {
      cor_mat <- cor(t(x$X[genes, ]) %>% as.matrix(), method = "spearman")
    }
    
    eigen_df <- foreach(i = 1:nrow(bait_df), .combine = "rbind") %do% {
      bait <- x$bait_sets[[i]]
      label <- paste0("Bait ", i, ": ", round(bait_df$tightness[i], 2))
      
      # remove any rows that are missing 
      cor_mat_tmp <- cor_mat[c(rand_genes, bait), 
                             c(rand_genes, bait)]
      
      # see if any columns have zero variance 
      cols <- which(is.na(matrixStats::colSds(cor_mat_tmp, na.rm = TRUE)))
      if(length(cols) > 0){
        cor_mat_tmp <- cor_mat_tmp[-cols, -cols]
      }
      
      eigen_space <- getUMAPCoordinates(cor_mat_tmp, 2) %>% 
        data.frame() 
      
      
      eigen_space %>% dplyr::mutate(bait = ifelse(rownames(eigen_space) %in% 
                                                    as.character(bait), "bait", "random")) %>%
        dplyr::mutate(gene = rownames(eigen_space[['coordinates']]),
                      label = label)
    }
    
    p <- ggplot(eigen_df) + 
      geom_point(aes(x = X1, y = X2, color = bait), 
                 alpha = 0.5) + 
      theme_bw() + scale_color_discrete(name = "Gene type") +
      facet_wrap(~label) + theme(legend.position = "top") + 
      labs(x = "UMAP.1", y = "UMAP.2")
    
    return(p)
  }
  
}

