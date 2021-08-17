# appropriate output types for good input
test_that("geneFishing returns list with appropriate arguments", {
  X <- matrix(rnorm(200), nrow = 20)
  rownames(X) <- paste0("R", 1:20)
  bait <- paste0("R", 1:5)
  geneFishing(X, bait)
})


# appropriate output for good input 
test_that("geneFishing returns list with appropriate values", {
  
})

# error handling
test_that("geneFishing returns errors when appropriate", {
  expect_error(geneFishing("A", "B"))
  
  # no overlap in gene names with row names of matrix 
  X <- matrix(rnorm(200), nrow = 20)
  expect_error(geneFishing(X, c("A", "B", "C", "D"), alpha = 3, 
                           min_bait_genes = 5), 
               paste0("The interesection of bait_genes and row names of ",
                      "expression matrix is less than ", 5, 
                      "\nConsider lowering min_bait_genes or using different bait."))
  
  # non-matrix input
  X <- data.frame(matrix(rnorm(200), nrow = 20))
  rownames(X) <- paste0("R", 1:20)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 3),
               paste0("X must be a matrix, or in one of the following classes:\n", 
                      "SingleCellExperiment, dgCMatrix, dgTMatrix, dgRMatrix, dgeMatrix"))
  
  # get error when there are not enough genes 
  X <- matrix(rnorm(200), nrow = 20)
  rownames(X) <- paste0("R", 1:20)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 5),
               paste0("There are not enough genes (rows) provided.  Need to have ", 
                      "at least alpha * length(bait_genes) rows in X."))
  
  # produce error when there is only one column
  X <- matrix(rnorm(200), nrow = 200)
  rownames(X) <- paste0("R", 1:200)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 3, umap = FALSE), 
               "X needs to have more than 1 column")
  
  # get error when the bait is not tight enough
  X <- matrix(rnorm(200 * 10), nrow = 200)
  rownames(X) <- paste0("R", 1:200)
  expect_error(geneFishing(X, paste0("R", 1:5), alpha = 3, umap = FALSE),
               paste0(
                 "Inputted bait has tightness less than min_tightness, ", 
                 "consider using probeFishability function. ", 
                 "\nYou can also increase min_tightness, and try using the same bait.\n"))
  
})


