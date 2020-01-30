# Gviz methods

test_that("General accessors work", {
  ir <- IRanges(1,1)
  gr <- GRanges("chr1", ir, score=1)
  dt <- DataTrack()
  
  ranges(dt) <- gr
  expect_identical(ranges(dt), gr)
  expect_identical(range(dt), ir)
  
  gt <- GenomeAxisTrack()
  ranges(gt) <- gr
  expect_identical(ranges(gt), gr)
  expect_identical(range(gt), ir)
  
})
