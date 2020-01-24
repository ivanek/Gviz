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
  dp <- DisplayPars(a=1)
  
  expect_s4_class(dp, "DisplayPars")
  expect_identical(dp@pars$a, 1)
  expect_error(DisplayPars(1), "All supplied arguments must be named.")
  expect_error(DisplayPars(1, b=2), "All supplied arguments must be named.")
  # .updateDp
  expect_message(Gviz:::.updateDp(dp), "Note that the behaviour of the")
  expect_identical(Gviz:::.updateDp(dp, interactive=FALSE)@pars, list(a=1))
  dp@pars <- new.env()
  dp@pars[["a"]] <- 1
  expect_identical(Gviz:::.updateDp(dp, interactive=FALSE)@pars, list(a=1))
  expect_message(Gviz:::.updateDp(dp, interactive=FALSE), "The DisplayPars object has been updated")
  # getPar
  expect_identical(getPar(dp), list(a=1))
  expect_identical(getPar(dp, "a"), 1)
  # setPar
  expect_identical(setPar(dp, "b", 2, interactive=F)@pars$b, 2)
  expect_error(setPar(dp, "b", c(1,2), interactive=F), "equal length")
  # displayPars <-
  displayPars(dp) <- list(b=2)
  expect_identical(dp@pars, list(a=1, b=2))
  # getPar
  expect_identical(getPar(dp), list(a=1, b=2))
  expect_identical(getPar(dp, "a"), 1)
  # displayPars
  expect_identical(displayPars(dp), list(a=1, b=2))
  expect_identical(displayPars(dp, "a"), 1)
  # as.list
  expect_identical(as.list(dp), list(a=1, b=2))
})

# test_that("IdeogramTrack works", {
#   bands <- 
#   expect_s4_class(IdeogramTrack(genome="hg38", chromosome="chrX", bands=), "IdeogramTrack")
# })

test_that("GenomeAxisTrack works", {
  expect_s4_class(GenomeAxisTrack(), "GdObject")
  expect_s4_class(GenomeAxisTrack(), "GenomeAxisTrack")
})

test_that("DataTrack works", {
  expect_s4_class(DataTrack(), "GdObject")
  expect_s4_class(DataTrack(), "RangeTrack")
  expect_s4_class(DataTrack(), "NumericTrack")
  expect_s4_class(DataTrack(), "DataTrack")
})

test_that("AnnnotationTrack works", {
  expect_s4_class(AnnotationTrack(), "GdObject")
  expect_s4_class(AnnotationTrack(), "RangeTrack")
  expect_s4_class(AnnotationTrack(), "StackedTrack")
  expect_s4_class(AnnotationTrack(), "AnnotationTrack")
})

test_that("GeneRegionTrack works", {
  expect_s4_class(GeneRegionTrack(), "AnnotationTrack")
  expect_s4_class(GeneRegionTrack(), "GeneRegionTrack")
})

# test_that("BiomartGeneRegionTrack works", {
#   expect_s4_class(BiomartGeneRegionTrack(), "BiomartGeneRegionTrack")
# })

test_that("DetailsAnnotationTrack works", {
  expect_s4_class(DetailsAnnotationTrack(), "AnnotationTrack")
})

test_that("SequenceTrack works", {
  expect_s4_class(SequenceTrack(), "GdObject")
  expect_s4_class(SequenceTrack(), "SequenceTrack")
  expect_s4_class(SequenceTrack(), "SequenceDNAStringSetTrack")
})

test_that("AlignmentsTrack works", {
  expect_s4_class(AlignmentsTrack(), "GdObject")
  expect_s4_class(AlignmentsTrack(), "RangeTrack")
  expect_s4_class(AlignmentsTrack(), "AlignmentsTrack")
})

test_that("CustomTrack works", {
  expect_s4_class(CustomTrack(), "GdObject")
  expect_s4_class(CustomTrack(), "CustomTrack")
})
