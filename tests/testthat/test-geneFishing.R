# error handling
test_that("geneFishing returns errors when appropriate", {
  expect_error(geneFishing("A", "B"))
  
  # no overlap in gene names with row names of matrix 
  X <- matrix(rnorm(200), nrow = 20)
  expect_error(geneFishing(X, c("A", "B", "C", "D"), alpha = 3, 
                           min_bait_genes = 5), 
               paste0("The interesection of bait_genes and row names of ",
                      "expression matrix is less than ", 5, 
                      "\nConsider lowering min_bait_genes or using different bait."),
               fixed = TRUE)
  
  # non-matrix input
  X <- data.frame(matrix(rnorm(200), nrow = 20))
  rownames(X) <- paste0("R", 1:20)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 3),
               paste0("X must be a matrix, or in one of the following classes:\n", 
                      "SingleCellExperiment, dgCMatrix, dgTMatrix, dgRMatrix, dgeMatrix"),
               fixed = TRUE)
  
  # get error when there are not enough genes 
  X <- matrix(rnorm(200), nrow = 20)
  rownames(X) <- paste0("R", 1:20)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 5),
               paste0("There are not enough genes (rows) provided.  Need to have ", 
                      "at least alpha * length(bait_genes) rows in X."),
               fixed = TRUE)
  
  # produce error when there is only one column
  X <- matrix(rnorm(200), nrow = 200)
  rownames(X) <- paste0("R", 1:200)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 3, umap = FALSE), 
               "X needs to have more than 1 column",
               fixed = TRUE)
  
  # get error when the bait is not tight enough
  X <- matrix(rnorm(200 * 10), nrow = 200)
  rownames(X) <- paste0("R", 1:200)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 3, umap = FALSE),
               paste0(
                 "Inputted bait has tightness less than min_tightness, ", 
                 "consider using probeFishability function. ", 
                 "\nYou can also increase min_tightness, and try using the same bait.\n"),
               fixed = TRUE)
  
})


# appropriate output types for good input
test_that("geneFishing works with valid input for spectral coordinates", {
  # generate two random clouds of data
  X <- rbind(mvrnorm(100, mu = rep(4, 10), 
                     Sigma = diag(rep(0.25, 10))), 
             mvrnorm(1000, mu = rep(1, 10), 
                     Sigma = diag(rep(1, 10))))
  
  rownames(X) <- paste0("R", 1:nrow(X))
  
  # use first 20 genes as bait
  bait <- paste0("R", 1:20)
  
  # euclidean distance
  output <- geneFishing(X, bait, alpha = 5, method = "euclidean", umap = FALSE, 
                        fishing_rounds = 50, n_probing_rounds = 50)
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # cosine similarity
  # generate data
  m <- 110
  n <- 2000
  
  # bait gene correlation 
  bait_gene_cor <- matrix(runif(10 * 10, 0.5, 1), nrow = 10) 
  diag(bait_gene_cor) <- 1
  bait_gene_cor[lower.tri(bait_gene_cor)] = t(bait_gene_cor)[
    lower.tri(bait_gene_cor)]
  non_bait_gene_cor <- matrix(runif(100 * 100, -1, 1), nrow = 100)
  diag(non_bait_gene_cor) <- 1
  non_bait_gene_cor[lower.tri(non_bait_gene_cor)] = t(non_bait_gene_cor)[
    lower.tri(non_bait_gene_cor)]
  
  sigma <- rbind(cbind(bait_gene_cor, matrix(runif(10 * 100, -0.1, 0.1), nrow = 10)),
                 cbind(matrix(runif(10 * 100, -1, 1), nrow = 100), non_bait_gene_cor)) 
  sigma <- as.matrix(Matrix::nearPD(sigma, corr = TRUE, keepDiag = TRUE)$mat)
  
  X <- mvrnorm(200, mu = rep(0, m), Sigma = sigma, empirical = T)
  X <- t(X)[, 1:50]
  rownames(X) <- paste0("R", 1:nrow(X))
  
  # use first 20 genes as bait
  bait <- paste0("R", 1:7)
  
  output <- geneFishing(X, bait, alpha = 5, method = "cosine", umap = FALSE, 
                        fishing_rounds = 50, n_probing_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # spearman 
  output <- geneFishing(X, bait, alpha = 5, method = "spearman", umap = FALSE, 
                        fishing_rounds = 50, n_probing_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # pearson
  output <- geneFishing(X, bait, alpha = 5, method = "pearson", umap = FALSE, 
                        fishing_rounds = 50, n_probing_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
})

test_that("geneFishing works with valid input for UMAP coordinates", {
  # cosine similarity
  # generate data
  m <- 110
  n <- 2000
  
  # bait gene correlation 
  bait_gene_cor <- matrix(runif(10 * 10, 0.5, 1), nrow = 10) 
  diag(bait_gene_cor) <- 1
  bait_gene_cor[lower.tri(bait_gene_cor)] = t(bait_gene_cor)[
    lower.tri(bait_gene_cor)]
  non_bait_gene_cor <- matrix(runif(100 * 100, -1, 1), nrow = 100)
  diag(non_bait_gene_cor) <- 1
  non_bait_gene_cor[lower.tri(non_bait_gene_cor)] = t(non_bait_gene_cor)[
    lower.tri(non_bait_gene_cor)]
  
  sigma <- rbind(cbind(bait_gene_cor, matrix(runif(10 * 100, -0.1, 0.1), nrow = 10)),
                 cbind(matrix(runif(10 * 100, -1, 1), nrow = 100), non_bait_gene_cor)) 
  sigma <- as.matrix(Matrix::nearPD(sigma, corr = TRUE, keepDiag = TRUE)$mat)
  
  X <- mvrnorm(200, mu = rep(0, m), Sigma = sigma, empirical = T)
  X <- t(X)[, 1:50]
  rownames(X) <- paste0("R", 1:nrow(X))
  
  # use first 20 genes as bait
  bait <- paste0("R", 1:7)
  
  output <- geneFishing(X, bait, alpha = 5, method = "cosine", umap = TRUE, 
                        fishing_rounds = 50, n_probing_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # spearman 
  output <- geneFishing(X, bait, alpha = 5, method = "spearman", umap = TRUE, 
                        fishing_rounds = 50, n_probing_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # pearson
  output <- geneFishing(X, bait, alpha = 5, method = "pearson", umap = TRUE, 
                        fishing_rounds = 50, n_probing_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing")
  expect_equal(names(output), c("results", "fished_genes", 
                                "CFR_cutoff", "bait", "parameters"))
  # check all genes are in the output
  expect_equal(length(intersect(rownames(X), output$results$gene_id)), nrow(X))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
})

