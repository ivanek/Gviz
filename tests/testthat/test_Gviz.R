library(Gviz)

test_that("checking the class and structure works", {
  expect_null(Gviz:::.checkClass(data.frame(), "data.frame"))
  expect_error(Gviz:::.checkClass(data.frame(), "matrix"), "of class")
  expect_error(Gviz:::.checkClass(class="matrix"), "is missing with no default")
  expect_error(Gviz:::.checkClass(data.frame(), "data.frame", length=2), "of length")
})

test_that("conversion of chromosome names works", {
  expect_identical(Gviz:::.chrName(character()), character())
  options("ucscChromosomeNames" = TRUE)
  expect_identical(Gviz:::.chrName(c("1", "chr1", "M", "MT", "X")), c("chr1","chr1","chrM","chrM","chrX"))
  expect_identical(Gviz:::.chrName("foo", force=TRUE), "chrfoo")
  expect_error(Gviz:::.chrName("foo"), "Please consider setting options\\(ucscChromosomeNames=FALSE\\)")
  options("ucscChromosomeNames" = FALSE)
  expect_identical(Gviz:::.chrName(c("1", "chr1", "M", "X")), c("1","chr1","M","X"))
})