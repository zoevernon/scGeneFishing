#' 
#' plotSpectralCoordinates
#' 
#' Function to plot the eigen decomposition of a set of potential bait genes
#' and a randomly sample group of additional genes. Will plot the UMAP of those
#' coordinates if 
#' 
#' @param X matrix of gene expression where the rows are genes and the 
#' columns are samples (cells or individuals). 
#' @param potential_bait_genes set of genes for which to check if there is a 
#' subset that can be used as bait.  These genes should be a subset of the 
#' rownames of \code{X}. 
#' @param alpha controls number of random genes that are sampled in each round 
#' of fishing.  It will sample \code{alpha} times the number of genes in 
#' \code{bait_genes}.  The default is 5.  The stronger the bait separates from 
#' samples of random genes the larger \code{alpha} the algorithm can handle.
#' @param umap indicator of whether the computation should be done using UMAP.
#' The default is TRUE.  When umap = FALSE it will use spectral coordinates. 
#' @param k number of eigen vectors to use in UMAP.
#' @param method input to \code{cor()} function telling what type of correlation
#' to use. Defaults to spearman. 
#' 
#' @return Plot of gene projection.
#' 
#' @examples
#' 
#' @import ggplot2
#' @export
#' 
plotSpectralCoordinates <- function(X, potential_bait, alpha = 5, umap = TRUE, 
                                    method = c("spearman", "pearson", "cosine",
                                               "euclidean",  "maximum", "manhattan", 
                                               "canberra", "binary", "minkowski")){
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
  
  # look at plots (will plot the number of random genes corresponding to 
  # alpha * number of genes in tightest bait)
  n_random <- alpha * length(potential_bait)
  rand_genes <- sample(setdiff(rownames(X), potential_bait), n_random)
  genes <- c(rand_genes,potential_bait)
  
  if(is.matrix(X)){
    X <- X[genes, ]
  } else if(class(X) == "SingleCellExperiment"){
    X <- logcounts(X)[genes, ] %>% as.matrix()
  } else if(class(X) %in% c("dgCMatrix", "dgTMatrix", "dgRMatrix", "dgeMatrix")) {
    X <- X[genes, ] %>% as.matrix()
  }
  
  if(method %in% c("euclidean",  "maximum", "manhattan", "canberra", 
                   "binary", "minkowski")){
    cor_mat <- dist(X, method = method) %>% as.matrix()
    cor_mat <- 1 / (cor_mat + 1)
  } else if(method == "cosine"){
    cor_mat <- cosineSimMatrix(X)
  } else {
    cor_mat <- cor(t(X), method = method)
  }
  
  # remove any rows that are missing 
  cor_mat_tmp <- cor_mat[c(rand_genes, potential_bait), 
                         c(rand_genes, potential_bait)]
  
  # see if any columns have zero variance 
  cols <- which(is.na(matrixStats::colSds(cor_mat_tmp, na.rm = TRUE)))
  if(length(cols) > 0){
    cor_mat_tmp <- cor_mat_tmp[-cols, -cols]
  }
  
  if(!umap){
    eigen_df <- getSpectralCoordinates(cor_mat_tmp, 2)$coordinates %>% 
      data.frame()
    
  }else{
    eigen_df <- getUMAPCoordinates(cor_mat_tmp, 2) %>% 
      data.frame() %>% rename(eigen.1 = X1, eigen.2 = X2)
  }
  
  
  
  eigen_df <- eigen_df %>%
    dplyr::mutate(bait = ifelse(rownames(eigen_df) %in% 
                                  as.character(potential_bait), 
                                "Potential bait", "Non-bait")) %>%
    dplyr::mutate(gene = rownames(eigen_df[['coordinates']]))
  
  if(!umap){
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
    theme(legend.position = "top") + 
    labs(x = x_lab, y = y_lab)
  
  return(p)
}
