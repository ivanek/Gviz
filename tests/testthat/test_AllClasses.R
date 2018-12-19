library(Gviz)

context("Gviz classes")

test_that("building ImageMap works", {
  mat <- matrix(1, ncol=4, dimnames=list("a",NULL))
  tags <- list(a=c(a="tag"))
  expect_error(Gviz:::ImageMap(matrix(1, ncol=3), tags=tags), "must be a numeric matrix with 4 column")
  expect_error(Gviz:::ImageMap(matrix(1, ncol=4), tags=tags), "Rownames must be set for the matrix in")
  expect_error(Gviz:::ImageMap(mat, tags=list(c(a="tag1"))), "must be a named list with character vector items.")
  expect_error(Gviz:::ImageMap(mat, tags=list(a=c(a="tag1", b="tag2"))), "following values in the")
  expect_error(Gviz:::ImageMap(mat, tags=list(a=c("tag1"))), "items in the 'tags' list must be named character")
  # check the class
  expect_s4_class(Gviz:::ImageMap(mat, tags=tags), "ImageMap")
  # check the slots
  expect_identical(Gviz:::ImageMap(mat, tags=tags)@coords,mat) 
  expect_identical(Gviz:::ImageMap(mat, tags=tags)@tags,tags)
})
