#' 
#' geneFishing
#' 
#' Function to perform gene fishing on a gene expression matrix. 
#' 
#' @param X matrix of gene expression where the rows are genes and the 
#' columns are samples (cells or individuals). 
#' @param bait_genes genes to use as bait to fish out similar genes from the 
#' \code{X}.  These genes should be a subset of the rownames of \code{X}. It 
#' can also be the output of the \code{probeFishability} function.  In that case
#' it will by default use the tightest bait set.  To change bait set you can use
#' the \code{bait_index} parameter.
#' @param alpha controls number of random genes that are sampled in each round 
#' of fishing.  It will sample \code{alpha} times the number of genes in 
#' \code{bait_genes}.  The default is 5.  The stronger the bait separates from 
#' samples of random genes the larger \code{alpha} the algorithm can handle.
#' @param fishing_rounds number of rounds of fishing to do.  The default is 1,000.
#' @param k number of clusters to use in geneFishing.  The default is 2, however
#' if you fished out too many genes with \code{k = 2} it can help to increase to 
#' \code{k = 3}. 
#' @param n_probing_rounds number of random samples to be used in assessing fishability
#' of a set of genes.  The default is 100.  . 
#' @param min_tightness parameter that controls how tight the set of genes needs
#' to be to be considered as bait.  Tightness is measured computing the ratio of
#' the median distance between all potential bait genes to the distance between the
#' medoid of the potential bait genes and a set of randomly sampled genes. The
#' average of these tightness measures are taken over \code{n_probing_rounds} random samples
#' of genes.  
#' @param umap indicator of whether the computation should be done using UMAP.
#' The default is FALSE.  When umap = TRUE the computation will take longer.
#' @param ncores number of cores to do parallel computation on if parallel = TRUE.
#' @param min_bait_genes minimum number of genes for bait set to be viable. 
#' @param bait_index index of \code{bait_sets} list from \code{probeFishability()}, to 
#' use.  Defaults to 1, which corresponds to the tightest bait.
#' @param method input to \code{cor()} or \code{dist()} functions telling what 
#' type of correlation to use. Cosine similarity also a possibility. 
#' Defaults to spearman.  This must be one of "spearman", "pearson", "cosine",
#' "euclidean",  "maximum", "manhattan", "canberra", "binary", or "minkowski". 
#' 
#' @return A list with the following items 
#' \itemize{
#'   \item \code{results}: data.frame with the capture frequency rate for all
#'   genes. 
#'   \item \code{fished_genes}: vector of fished out genes (includes bait genes).
#'   \item \code{CFR_cutoff}: cutoff on CFR that was used to define the fished
#'   out gene set.  One can always choose their own cutoff and determine fished
#'   out genes from the \code{results} data.frame. 
#'   \item \code{bait}: vector of bait genes that were used for doing the GeneFishing.
#'   \item \code{parameters}: list of parameters used in the fishing.  
#' }
#' 
#' @examples
#' 
#' @import foreach
#' @import dplyr
#' @importClassesFrom Matrix dgCMatrix dgeMatrix dgTMatrix dgRMatrix
#' @importFrom SummarizedExperiment assay
#' @export

geneFishing <- function(X, bait_genes, alpha = 5, fishing_rounds = 1000, 
                        k = 2, min_tightness = 0.5, n_probing_rounds = 100, 
                        umap = TRUE, ncores = 2, min_bait_genes = 5, 
                        bait_index = 1, 
                        method = c("spearman", "pearson", "cosine",
                                   "euclidean",  "maximum", "manhattan", 
                                   "canberra", "binary", "minkowski")){
  
  doParallel::registerDoParallel(ncores)
  
  # check correlation method provided is correct 
  method <- match.arg(method)
  using_dist <- ifelse(method %in% c("euclidean", "maximum", "manhattan", 
                                     "canberra", "binary", "minkowski"), 
                       TRUE, FALSE)
  
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
  
  # check that there are alpha * length(bait_genes) in the data
  assertthat::assert_that(
    (length(bait_genes) * alpha) <= (nrow(X) - length(bait_genes)), 
    msg = paste0("There are not enough genes (rows) provided.  Need to have ", 
                 "at least alpha * length(bait_genes) rows in X."))
  
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
  
  # check class of bait and if it is from bait probing just use the best 
  # bait 
  if(class(bait_genes) == "gene_fishing_probe") {
    used_probe_output <- TRUE
    if(length(bait_genes$bait_sets) > 1){
      message(paste0("Bait genes provided from probeFishability. Using bait set ", 
                     bait_index, ". To use a different bait set, change bait_index."))
    }
    
    # update method 
    method <- bait_genes$method
    bait_genes <- bait_genes$bait_sets[[bait_index]]
  }else{
    used_probe_output <- FALSE
  }
  
  # make sure bait genes are in expression matrix 
  bait_genes_orig <- bait_genes
  bait_genes <- bait_genes[bait_genes %in% rownames(X)] %>% 
    as.character()
  
  # make sure there are at least min_bait_genes in bait
  assertthat::assert_that(
    length(bait_genes) >= min_bait_genes, 
    msg = paste0("The interesection of bait_genes and row names of ",
                 "expression matrix is less than ", min_bait_genes, 
                 "\nConsider lowering min_bait_genes or using different bait."))
  
  # return a warning that some genes were removed as they are not in the row-
  # names of the expression matrix
  if(length(bait_genes) < length(bait_genes_orig)) {
    warning(paste0("Inputted bait had ", length(bait_genes_orig),
                   " genes, only ", length(bait_genes), 
                   " are in the row names of the inputted expression matrix."))
  }
  
  # extract names of all genes other than bait
  all_genes <- rownames(X)
  pool_genes <- setdiff(all_genes, bait_genes)
  
  # check what type of data was inputted
  # if user inputted bait not from probeFishability, make sure that the 
  # provided bait set is tight enough
  # note there is this long ifelse statement that calls different DB index 
  # functions because it is slower to have the ifelse statements inside the loop
  # that computes the correlation matrix which is in that function
  if(used_probe_output == FALSE){
    if(using_dist){
      db_index <- computeAvgDBIndexDist(bait_genes, X, 
                                       k = k,
                                       n_rounds = n_probing_rounds, 
                                       alpha = alpha,
                                       method = method) %>% mean(na.rm = TRUE)
    }else if(method == "cosine"){
      if(umap){
        db_index <- computeAvgDBIndexCosUMAP(bait_genes, X, 
                                             k = k,
                                             n_rounds = n_probing_rounds, 
                                             alpha = alpha,
                                             method = method) %>% 
          mean(na.rm = TRUE)
      }else{
        db_index <- computeAvgDBIndexCosSpectral(bait_genes, X, 
                                                 k = k,
                                                 n_rounds = n_probing_rounds, 
                                                 alpha = alpha,
                                                 method = method) %>% 
          mean(na.rm = TRUE)
      }
    }else if(umap){
      db_index <- computeAvgDBIndexUMAP(bait_genes, X, 
                                        k = k,
                                        n_rounds = n_probing_rounds, 
                                        alpha = alpha,
                                        method = method) %>% 
        mean(na.rm = TRUE)
    }else{
      db_index <- computeAvgDBIndexSpectral(bait_genes, X, 
                                            k = k,
                                            n_rounds = n_probing_rounds, 
                                            alpha = alpha,
                                            method = method) %>% 
        mean(na.rm = TRUE)
    }

    # check that the bait is tight enough
    assertthat::assert_that(db_index < min_tightness, msg = paste0(
      "Inputted bait has tightness less than min_tightness, ", 
      "consider using probeFishability function. ", 
      "\nYou can also increase min_tightness, and try using the same bait.\n"))
  }
  
  # do gene fishing in UMAP 
  if(using_dist){
    fish_freq_df <- geneFishingDist(X, bait_genes, pool_genes,
                                    alpha, fishing_rounds, k, method)
  }else if(method == "cosine"){
    if(umap){
      fish_freq_df <- geneFishingCosUMAP(X, bait_genes, pool_genes,
                                         alpha, fishing_rounds, k, method)
    }else{
      fish_freq_df <- geneFishingCosSpectral(X, bait_genes, pool_genes,
                                             alpha, fishing_rounds, k, method)
    }
  }else if(umap){
    fish_freq_df <- geneFishingUMAP(X, bait_genes, pool_genes,
                                    alpha, fishing_rounds, k, method)
  }else{
    fish_freq_df <- geneFishingSpectral(X, bait_genes, pool_genes,
                                        alpha, fishing_rounds, k, method)
  }
  
  
  # automatically select cutoff
  cutoff <- getCutoff(fish_freq_df$CFR, 
                      which(fish_freq_df$gene_id %in% bait_genes),
                      unit = 1 / fishing_rounds)
  
  # if cutoff is less than 0.9 warn them it was probably not a good 
  # bait set 
  if(cutoff < 0.9){
    warning(paste0("CFR cutoff is less than ", 0.9, 
                   ". This typically indicates a bait set that is not tight enough.", 
                   "\nTry increasing k or using a new bait set."))
  }
  
  # check if there were any non-bait genes that were fished out above the cutoff
  if(all(fish_freq_df$bait[fish_freq_df$CFR >= cutoff] == "yes")){
    warning(paste0("All fished_genes are bait. You can still use results", 
                   " data.frame to rank non-bait genes \nwith CFR, but none have", 
                   " CFR above CFR_cutoff"))
  }
  
  # get list to return as gene_fishing class 
  final_output <- list(
    results = fish_freq_df, 
    fished_genes = fish_freq_df$gene_id[fish_freq_df$CFR >= cutoff],
    CFR_cutoff = cutoff, 
    bait = bait_genes,
    parameters = list(k = k, fishing_rounds = fishing_rounds, 
                      n_probing_rounds = n_probing_rounds,
                      min_tightness = min_tightness,
                      min_bait_genes = min_bait_genes,
                      method = method)
  )
  class(final_output) <- "gene_fishing"
  
  return(final_output)
}

# print and plotting methods
#' @export
print.gene_fishing <- function(x, ...){
  cat(paste("Fished out", length(x$fished_genes),
            "with a CFR cutoff of", x$CFR_cutoff, "\n\n"))
  cat("$bait\n")
  cat(sort(x$bait), "\n\n")
  cat(paste("$results\n"))
  print.data.frame(head(x$results, 10))
}

#' @export
plot.gene_fishing<- function(x, ...){
  
  p <- ggplot(x$results) + 
    geom_histogram(aes(CFR), color = "grey", size = 0.1, bins = 50) + 
    theme_bw() 
  
  return(p)
}

