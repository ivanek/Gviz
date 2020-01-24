library(Gviz)

test_that("ImageMap works", {
  mat <- matrix(1, ncol=4, dimnames=list("a",NULL))
  tags <- list(a=c(a="tag"))
  expect_s4_class(Gviz:::ImageMap(mat, tags=tags), "ImageMap")
  expect_identical(Gviz:::ImageMap(mat, tags=tags)@coords,mat) 
  expect_identical(Gviz:::ImageMap(mat, tags=tags)@tags,tags)
  expect_error(Gviz:::ImageMap(matrix(1, ncol=3), tags=tags), "must be a numeric matrix with 4 column")
  expect_error(Gviz:::ImageMap(matrix(1, ncol=4), tags=tags), "Rownames must be set for the matrix in")
  expect_error(Gviz:::ImageMap(mat, tags=list(c(a="tag1"))), "must be a named list with character vector items.")
  expect_error(Gviz:::ImageMap(mat, tags=list(a=c(a="tag1", b="tag2"))), "following values in the")
  expect_error(Gviz:::ImageMap(mat, tags=list(a=c("tag1"))), "items in the 'tags' list must be named character")
})

test_that("DisplayPars works", {
  expect_s4_class(DisplayPars(), "DisplayPars")
  expect_identical(DisplayPars(a=1)@pars$a, 1)
  expect_error(DisplayPars(1), "All supplied arguments must be named.")
  expect_error(DisplayPars(1, b=2), "All supplied arguments must be named.")
})

# test_that("IdeogramTrack works", {
#   bands <- 
#   expect_s4_class(IdeogramTrack(genome="hg38", chromosome="chrX", bands=), "IdeogramTrack")
# })

test_that("GenomeAxisTrack works", {
  expect_s4_class(GenomeAxisTrack(), "GenomeAxisTrack")
})

test_that("DataTrack works", {
  expect_s4_class(DataTrack(), "DataTrack")
})

test_that("AnnnotationTrack works", {
  expect_s4_class(AnnotationTrack(), "AnnotationTrack")
})

test_that("GeneRegionTrack works", {
  expect_s4_class(GeneRegionTrack(), "GeneRegionTrack")
})

# test_that("BiomartGeneRegionTrack works", {
#   expect_s4_class(BiomartGeneRegionTrack(), "BiomartGeneRegionTrack")
# })

test_that("DetailsAnnotationTrack works", {
  expect_s4_class(DetailsAnnotationTrack(), "AnnotationTrack")
})

test_that("SequenceTrack works", {
  expect_s4_class(SequenceTrack(), "SequenceTrack")
})

test_that("AlignmentsTrack works", {
  expect_s4_class(AlignmentsTrack(), "AlignmentsTrack")
})
