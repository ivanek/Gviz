# Gviz methods

test_that("General accessors work", {

  expect_identical(ranges(dataTrack), granges(gr))
  expect_identical(range(dataTrack), ir)
  
  expect_identical(ranges(axisTrack), gr)
  expect_identical(range(axisTrack), ir)
  
  expect_identical(seqnames(annoTrack), "chr1")
  expect_identical(seqnames(dataTrack), "chr1")
  expect_identical(seqnames(seqTrack.dna), "chr1")
  expect_identical(seqnames(seqTrack.rna), "chr1")
  expect_identical(seqnames(seqTrack.bs), seqnames(BSgenome.Hsapiens.UCSC.hg19))

  expect_identical(seqlevels(annoTrack), "chr1")
  expect_identical(seqlevels(dataTrack), "chr1")
  expect_identical(seqlevels(seqTrack.dna), "chr1")
  expect_identical(seqlevels(seqTrack.rna), "chr1")
  expect_identical(seqlevels(seqTrack.bs), seqlevels(BSgenome.Hsapiens.UCSC.hg19))

  expect_identical(seqinfo(annoTrack), table("chr1"))
  
  expect_identical(min(dataTrack), 1L)
  expect_identical(max(dataTrack), 5L)
  
  expect_identical(start(dataTrack), 1L)
  expect_identical(end(dataTrack), 5L)
  expect_identical(width(dataTrack), 5L)
  expect_identical(length(dataTrack), 1L)
  
  expect_identical(start(ideoTrack), NULL)
  expect_identical(end(ideoTrack), NULL)
  expect_identical(width(ideoTrack), NULL)
  expect_identical(length(ideoTrack), 4L)
  
  expect_identical(start(axisTrack), 1L)
  expect_identical(end(axisTrack), 5L)
  expect_identical(width(axisTrack), 5L)
  
  expect_identical(start(seqTrack.dna), NULL)
  expect_identical(end(seqTrack.dna), NULL)
  expect_identical(width(seqTrack.dna), NULL)
  expect_identical(length(seqTrack.dna), 100L) # !
  
  expect_identical(length(highTrack), 1L)
  expect_identical(length(overTrack), 2L)
})

test_that("General replacement methods work", {
  ranges(dataTrack) <- gr2
  expect_identical(ranges(dataTrack), gr2)
  expect_identical(range(dataTrack), ir2)
  
  ranges(axisTrack) <- gr2
  expect_identical(ranges(axisTrack), gr2)
  expect_identical(range(axisTrack), ir2)
 
  ranges(annoTrack) <- gr2
  expect_identical(ranges(annoTrack), gr2)
  expect_identical(range(annoTrack), ir2)
  
  start(annoTrack) <- 3L
  expect_identical(start(annoTrack), 3L)
  start(axisTrack) <- 3L
  expect_identical(start(axisTrack), 3L)
  start(ideoTrack) <- 3L
  expect_identical(start(ideoTrack), NULL)
  
  end(annoTrack) <- 4L
  expect_identical(end(annoTrack), 4L)
  end(axisTrack) <- 4L
  expect_identical(end(axisTrack), 4L)
  end(ideoTrack) <- 4L
  expect_identical(end(ideoTrack), NULL)
  
  width(annoTrack) <- 10L
  expect_identical(width(annoTrack), 10L)
})


test_that("chromosome accessors and replacement methods work", {
  expect_identical(chromosome(annoTrack), "chr1")
  expect_identical(chromosome(dataTrack), "chr1")
})