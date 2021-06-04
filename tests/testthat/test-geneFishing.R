test_that("geneFishing returns errors when appropriate", {
  expect_error(geneFishing("A", "B"))
})