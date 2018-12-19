library(Gviz)

context("Gviz classes")

test_that("checking the class and structure works", {
  mat <- matrix(1, ncol=4, dimnames=list("a",NULL))
  tags <- list(a=c("a"="tag"))
  expect_error(Gviz:::ImageMap(matrix(1, ncol=3), tags=tags), "must be a numeric matrix with 4 column")
  expect_error(Gviz:::ImageMap(matrix(1, ncol=4), tags=tags), "Rownames must be set for the matrix in")
  expect_identical(Gviz:::ImageMap(mat, tags=tags), 
               new("ImageMap", coords = mat, tags=tags))
})
