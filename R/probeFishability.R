#' 
#' probeFishability
#' 
#' Function to find a set of genes that can be used as bait in geneFishing from 
#' a larger set of potential bait genes. 
#' 
#' @param X matrix of gene expression where the rows are genes and the 
#' columns are samples (cells or individuals). 
#' @param potential_bait_genes set of genes for which to check if there is a 
#' subset that can be used as bait.  These genes should be a subset of the 
#' rownames of \code{X}. 
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
#' @param method input to \code{cor()} function telling what type of correlation
#' to use. Defaults to spearman. 
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
#' @import dplyr
#' @import ggplot2
#' @importClassesFrom Matrix dgCMatrix dgeMatrix dgTMatrix dgRMatrix
#' @importFrom SummarizedExperiment assay assayNames
#' @export

probeFishability <- function(X, potential_bait, n_rounds = 100, alpha = 5,
                             min_tightness = 0.5, min_genes = 5, n_neighbors = 15,
                             umap = TRUE, ncores = 2, 
                             method = c("spearman", "pearson", "cosine",
                                        "euclidean",  "maximum", "manhattan", 
                                        "canberra", "binary", "minkowski")){
  doParallel::registerDoParallel(ncores)
  
  # check correlation method provided is correct 
  method <- match.arg(method)
  
  # assertion to make sure the entry is a matrix, sparse matrix or sce
  assertthat::assert_that(
    is.matrix(X) | any(class(X) %in% c("SingleCellExperiment", "dgCMatrix", 
                                       "dgTMatrix", "dgRMatrix", "dgeMatrix",
                                       "SummarizedExperiment")), 
    msg = paste0("X must be a matrix, or in one of the following classes:\n", 
                 "SingleCellExperiment, dgCMatrix, dgTMatrix, dgRMatrix, dgeMatrix"))
  
  if(any(class(X) %in% c("SingleCellExperiment", "SummarizedExperiment"))){
    assertthat::assert_that(
      any(assayNames(X) == "logcounts"), 
      msg = paste0("There must be a logcounts assay."))
    X_sce <- X
    X <- assay(X, expr_values = "logcounts")
  }
  
  # check that there are alpha * length(potential_bait_genes) in the data
  assertthat::assert_that(
    (length(potential_bait) * alpha) <= (nrow(X) - length(potential_bait)), 
    msg = paste0("There are not enough genes (rows) provided.  Need to have ", 
                 "at least alpha * length(potential_bait) rows in X."))
  
  # needs to be more than 1 column
  assertthat::assert_that(
    ncol(X) > 1, msg = "X needs to have more than 1 column")
  
  # remove any rows and columns that are all 0 
  all_zero_rows <- which(rowSums(X) == 0)
  all_zero_cols <- which(colSums(X) == 0)
  if(length(all_zero_rows) > 0){
    message(paste0("There were ", length(all_zero_rows), 
                   " rows in X with all zeroes.  They will be removed."))
    X <- X[-all_zero_rows, ]
  }
  if(length(all_zero_cols) > 0){
    message(paste0("There were ", length(all_zero_rows), 
                   " columns in X with all zeroes.  They will be removed."))
    X <- X[, -all_zero_cols]
  }

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
    results <- probeFishabilityUMAP(X, 
                                    potential_bait,
                                    n_rounds = n_rounds, 
                                    alpha = alpha, 
                                    min_tightness = min_tightness, 
                                    min_genes = min_genes,
                                    n_neighbors = n_neighbors,
                                    method = method)
  }else{
    results <- probeFishabilitySpectral(X, 
                                        potential_bait,
                                        n_rounds = n_rounds, 
                                        alpha = alpha, 
                                        min_tightness = min_tightness, 
                                        min_genes = min_genes,
                                        method = method)
  }
  
  
  return(results)
}

#' @export
print.gene_fishing_probe <- function(x, ...){
  if(length(x$best_bait) > 5){
    cat(paste(length(x$best_bait), "in tightest bait set:\n"))
    cat(sort(x$best_bait)[1:5], ", ...\n\n")
    cat(paste("Found", nrow(x$bait_info) - 1, "additional bait sets.\n"))
  }else if(length(x$best_bait) == 0){
    cat("No bait found, try again with higher min_tightness.\n")
    cat("Note, it is not recommended to use min_tightness > 0.5.\n")
  }else{
    cat("Tightest bait set:\n")
    cat(sort(x$best_bait), "\n\n")
    cat(paste("Found", nrow(x$bait_info) - 1, "additional bait sets.\n"))
  }
  
}

#' @export
plot.gene_fishing_probe <- function(x, alpha = 5, 
                                    bait_indices = "all", ...){
  if(length(x$best_bait) == 0){
    cat("No bait found, try again with higher min_tightness.\n")
    cat("Note, it is not recommended to use min_tightness > 0.5.\n")
  }else{
    if(any(bait_indices == "all") | any(!bait_indices %in% 1:length(x$bait_sets))){
      bait_indices <- 1:length(x$bait_sets)
    }
    
    all_bait <- sapply(1:nrow(x$bait_info), 
                       function(i) x$bait_sets[[i]]) %>% unlist()
    bait_df <- x$bait_info
    
    # look at plots (will plot the number of random genes corresponding to 
    # alpha * number of genes in tightest bait)
    n_random <- alpha * length(x$bait_sets[[min(bait_indices)]])
    rand_genes <- sample(setdiff(rownames(x$X), all_bait), n_random)
    genes <- c(rand_genes, all_bait)
    
    if(is.matrix(x$X)){
      X <- x$X[genes, ]
    } else if(class(x$X) == "SingleCellExperiment"){
      X <- logcounts(x$X)[genes, ] %>% as.matrix()
    } else if(class(x$X) %in% c("dgCMatrix", "dgTMatrix", "dgRMatrix", "dgeMatrix")) {
      X <- x$X[genes, ] %>% as.matrix()
    }
    
    if(x$method %in% c("euclidean",  "maximum", "manhattan", "canberra", 
                       "binary", "minkowski")){
      cor_mat <- dist(X, method = x$method) %>% as.matrix()
      cor_mat <- 1 / (cor_mat + 1)
    } else if(x$method == "cosine"){
      cor_mat <- cosineSimMatrix(X)
    } else {
      cor_mat <- cor(t(X), method = x$method)
    }
    
    eigen_df <- foreach(i = bait_indices, .combine = "rbind") %do% {
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
      
      if(x$type == "spectral"){
        eigen_space <- getSpectralCoordinates(cor_mat_tmp, 2)$coordinates %>% 
          data.frame()
        
      }else{
        eigen_space <- getUMAPCoordinates(cor_mat_tmp, 2) %>% 
          data.frame() %>% rename(eigen.1 = X1, eigen.2 = X2)
      }

      
      
      eigen_space %>% dplyr::mutate(bait = ifelse(rownames(eigen_space) %in% 
                                                    as.character(bait), "bait", "random")) %>%
        dplyr::mutate(gene = rownames(eigen_space[['coordinates']]),
                      label = label)
    }
    
    if(x$type == "spectral"){
      x_lab <- "eigen.1"
      y_lab <- "eigen.2"
    }else{
      x_lab <- "UMAP.1"
      y_lab <- "UMAP.2"
    }
    
    p <- ggplot(eigen_df) + 
      geom_point(aes(x = eigen.1, y = eigen.2, color = bait), 
                 alpha = 0.5) + 
      theme_bw() + scale_color_discrete(name = "Gene type") +
      facet_wrap(~label) + theme(legend.position = "top") + 
      labs(x = x_lab, y = y_lab)
    
    return(p)
  }
 
}
