# Gviz methods

ir <- IRanges(1L,5L)
gr <- GRanges("chr1", ir, score=1)
ir2 <- IRanges(2L,6L)
gr2 <- GRanges("chr1", ir2, score=2)

dt <- DataTrack(gr)
gt <- GenomeAxisTrack(gr)
at <- AnnotationTrack(gr)


test_that("General accessors work", {

  ranges(dt) <- gr2
  expect_identical(ranges(dt), gr2)
  expect_identical(range(dt), ir2)
  
  ranges(gt) <- gr2
  expect_identical(ranges(gt), gr2)
  expect_identical(range(gt), ir2)
  
  expect_identical(seqnames(at), "chr1")
  expect_identical(seqnames(dt), "chr1")
  
  expect_identical(min(dt), 2L)
  expect_identical(max(dt), 6L)
})


test_that("chromosome accessor works", {
  expect_identical(chromosome(at), "chr1")
  expect_identical(chromosome(dt), "chr1")
})