# error handling
test_that("geneFishing returns errors when appropriate", {
  expect_error(probeFishability("A", "B"))
  
  # no overlap in gene names with row names of matrix 
  X <- matrix(rnorm(200), nrow = 20)
  expect_error(probeFishability(X, c("A", "B", "C", "D"), alpha = 3, 
                                min_genes = 5), 
               paste0("No intersection of potential_bait and row names of expression matrix"),
               fixed = TRUE)
  
  # non-matrix input
  X <- data.frame(matrix(rnorm(200), nrow = 20))
  rownames(X) <- paste0("R", 1:20)
  expect_error(probeFishability(X, paste0("R", 1:5), alpha = 3),
               paste0("X must be a matrix, or in one of the following classes:\n", 
                      "SingleCellExperiment, dgCMatrix, dgTMatrix, dgRMatrix, dgeMatrix"),
               fixed = TRUE)
  
  # get error when there are not enough genes 
  X <- matrix(rnorm(200), nrow = 20)
  rownames(X) <- paste0("R", 1:20)
  expect_error(probeFishability(X, paste0("R", 1:5), alpha = 5),
               paste0("There are not enough genes (rows) provided.  Need to have ", 
                      "at least alpha * length(potential_bait) rows in X."),
               fixed = TRUE)
  
  # produce error when there is only one column
  X <- matrix(rnorm(200), nrow = 200)
  rownames(X) <- paste0("R", 1:200)
  expect_error(probeFishability(X, paste0("R", 1:5), alpha = 3, umap = FALSE), 
               "X needs to have more than 1 column",
               fixed = TRUE)
  
  # get warning if there is no bait
  X <- matrix(rnorm(200 * 10), nrow = 200)
  rownames(X) <- paste0("R", 1:200)
  expect_warning(probeFishability(X, paste0("R", 1:5), alpha = 3, umap = FALSE),
               "No bait discovered.",
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
  bait <- paste0("R", 1:10)
  
  # euclidean distance
  output <- probeFishability(X, bait, alpha = 5, method = "euclidean", umap = FALSE, 
                             n_rounds = 50)
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
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
  
  # use first 20 genes as potential bait
  bait <- paste0("R", 1:20)
  
  output <- probeFishability(X, bait, alpha = 3, method = "cosine", umap = FALSE, 
                        n_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # spearman 
  output <- probeFishability(X, bait, alpha = 3, method = "spearman", umap = FALSE, 
                             n_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # pearson
  output <- probeFishability(X, bait, alpha = 3, method = "pearson", umap = FALSE, 
                             n_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
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
  
  # use first 20 genes as potential bait
  bait <- paste0("R", 1:20)
  
  output <- probeFishability(X, bait, alpha = 3, method = "cosine", umap = TRUE, 
                             n_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # spearman 
  output <- probeFishability(X, bait, alpha = 3, method = "spearman", umap = TRUE, 
                             n_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
  # pearson
  output <- probeFishability(X, bait, alpha = 3, method = "pearson", umap = TRUE, 
                             n_rounds = 50) 
  expect_type(output, "list")
  expect_s3_class(output, "gene_fishing_probe")
  expect_equal(names(output), c("best_bait", "bait_sets", "bait_info", "X",  
                                "method", "type", "rand_perct"))
  # check that we are returning a plot
  expect_s3_class(plot(output), c("gg", "ggplot"))
  
})

