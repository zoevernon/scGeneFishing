#' 
#' geneFishing
#' 
#' Function to perform gene fishing on a gene expression matrix. 
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
#' @param n_rounds number of random samples to be used in assessing fishability
#' of a set of genes.  The default is 100.  . 
#' @param min_tightness parameter that controls how tight the set of genes needs
#' to be to be considered as bait.  Tightness is measured computing the ratio of
#' the median distance between all potential bait genes to the distance between the
#' medoid of the potential bait genes and a set of randomly sampled genes. The
#' average of these tightness measures are taken over \code{n_rounds} random samples
#' of genes.  
#' @param umap indicator of whether the computation should be done using UMAP.
#' The default is TRUE.  When umap = FALSE it will use spectral coordinates. 
#' @param ncores number of cores to do parallel computation on if parallel = TRUE.
#' @param min_bait_genes minimum number of genes for bait set to be viable. 
#' 
#' @return A data.frame with the capture frequency rate for all genes
#' 
#' @examples
#' 
#' @import foreach
#' @import dplyr
#' @export

geneFishing <- function(exp_mat, bait_genes, alpha = 5, fishing_rounds = 1000, 
                        k = 2, min_tightness = 0.5, n_rounds = 100, umap = TRUE,
                        ncores = 2, min_bait_genes = 5){
  
  doParallel::registerDoParallel(ncores)
  
  # make sure bait genes are in expression matrix 
  bait_genes_orig <- bait_genes
  bait_genes <- bait_genes[bait_genes %in% rownames(exp_mat)] %>% 
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
  all_genes <- rownames(exp_mat)
  pool_genes <- setdiff(all_genes, bait_genes)
  
  # make sure that the provided bait set is tight enough 
  db_index <- ifelse(umap, 
                     computeAvgDBIndexUMAP(bait_genes, exp_mat, 
                                           n_rounds = n_rounds, 
                                           alpha = alpha) %>%
                       mean(),
                     computeAvgDBIndexSpectral(bait_genes, exp_mat, 
                                               n_rounds = n_rounds, 
                                               alpha = alpha) %>%
                       mean())
  
  # check that the bait is tight enough
  assertthat::assert_that(db_index < min_tightness, msg = paste(
    "Inputted bait is", round(db_index, 2), "which is less than min_tightness,", 
    "consider using probeFishability().", 
    "\nYou can also increase min_tightness, and use the same bait."))
    
  # do gene fishing in UMAP 
  if(umap){
    fish_freq_df <- geneFishingUMAP(exp_mat, bait_genes, pool_genes,
                                    alpha, fishing_rounds, k)
  }else{
    fish_freq_df <- geneFishingSpectral(exp_mat, bait_genes, pool_genes,
                                        alpha, fishing_rounds, k)
  }
  
  # CHANGE TO YUTINGS thing
  cutoff <- 0.99
  
  # get list to return as gene_fishing class 
  final_output <- list(
    results = fish_freq_df, 
    fished_genes = fish_freq_df$gene_id[fish_freq_df$CFR >= cutoff],
    CFR_cutoff = cutoff, 
    bait = bait_genes, 
    exp_mat = exp_mat
  )
  class(final_output) <- "gene_fishing"
  
  return(final_output)
}

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

