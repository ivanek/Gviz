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
  
  width(ideoTrack) <- 10L
  expect_identical(ideoTrack, ideoTrack)
})

test_that("values accessors and replacement methods work", {
  expect_identical(values(annoTrack), data.frame(feature = "unknown", group = "1", 
                               id = "unknown", density = 1, stringsAsFactors=FALSE))
  expect_identical(values(axisTrack), data.frame(score = 1, stringsAsFactors=FALSE))
  expect_identical(values(dataTrack), as.matrix(c(score = 1)))
  expect_identical(values(DataTrack()), matrix(logical(), 0, 0))
  expect_identical(values(alnTrack), NULL)
  
  values(dataTrack) <- matrix(3, dimnames=list("score", NULL))
  expect_identical(values(dataTrack), as.matrix(c(score = 3)))
  values(dataTrack) <- 3
  expect_identical(values(dataTrack), as.matrix(c(3)))
  expect_error(values(dataTrack) <- matrix("3", dimnames=list("score", NULL)), "Not numeric or dimensions of replacement value do not match.")
  expect_error(values(dataTrack) <- c(3,3), "Not numeric or invalid length of replacement vector.")
  
})

test_that("subseq works", {
  expect_error(subseq(SequenceTrack(), start=NA, end=NA, width=3), "Two out of the three in")
  expect_error(subseq(SequenceTrack(), start=10, end=1), "'end' has to be bigger than 'start'")
  expect_error(subseq(seqTrack.dna, start=10, end=1), "'end' has to be bigger than 'start'")
  expect_error(subseq(SequenceTrack(), start=1, end=10e6+2), "Sequence is too big")
  expect_error(subseq(seqTrack.bs, start=1, end=10e6+2), "Sequence is too big")
  expect_warning(subseq(seqTrack.dna, start=1, end=10, width=10), "are provided, ignoring")
  expect_identical(as.character(subseq(SequenceTrack(), start=1, end=10)), as.character(DNAString("----------")))
  expect_identical(as.character(subseq(seqTrack.dna, start=1, end=10)), as.character(DNAString("ATTTCCCTGA")))
  expect_identical(as.character(subseq(seqTrack.dna, start=1, width=10)), as.character(DNAString("ATTTCCCTGA")))
  expect_identical(as.character(subseq(seqTrack.dna, end=10, width=10)), as.character(DNAString("ATTTCCCTGA")))
  expect_identical(as.character(subseq(seqTrack.dna, start=91, end=110)), as.character(DNAString("ACGTCTTCCA----------")))
  expect_identical(as.character(subseq(seqTrack.bs, start=1, width=10)), as.character(DNAString("NNNNNNNNNN")))
  expect_identical(as.character(subseq(seqTrack.dna, start=1)), as.character(dna.sq[[1]]))
  expect_identical(as.character(subseq(SequenceTrack(), start=1)), as.character(DNAString("-")))
  expect_identical(as.character(subseq(seqTrack.dna, end=1)), as.character(DNAString("A")))
  expect_identical(as.character(subseq(SequenceTrack(), end=1)), as.character(DNAString("-")))
  displayPars(seqTrack.dna)$complement <- TRUE
  expect_identical(as.character(subseq(seqTrack.dna, start=1, end=10)), as.character(DNAString("TAAAGGGACT")))
  
  
  expect_error(subseq(SequenceTrack(fastafile, chromosome="chr1")), "at least two out of ")
  expect_error(subseq(SequenceTrack(fastafile, chromosome="chr1"), start=1), "at least two out of ")
  expect_warning(as.character(subseq(SequenceTrack(fastafile, chromosome="chr1"), start=1, end=10, width=10)), "All ")
  expect_error(as.character(subseq(SequenceTrack(fastafile, chromosome="chr1"), start=NA, end=NA, width=10)), "Two ")
  expect_identical(as.character(subseq(SequenceTrack(fastafile, chromosome="chr1"), start=1, width=10)), as.character(DNAString("CTANGAGACG")))
  expect_identical(as.character(subseq(SequenceTrack(fastafile, chromosome="chr1"), end=10, width=10)), as.character(DNAString("CTANGAGACG")))
  ## expect_identical(as.character(subseq(SequenceTrack(fastafile, chromosome="chr1"), start=1, end=NA, width=10)), as.character(DNAString("CTANGAGACG")))
  expect_identical(as.character(subseq(SequenceTrack(fastafile, chromosome="chr1"), start=1, end=10)), as.character(DNAString("CTANGAGACG")))
})

test_that("chromosome accessors and replacement methods work", {
  expect_identical(chromosome(annoTrack), "chr1")
  chromosome(annoTrack) <- "chr2"
  expect_identical(chromosome(annoTrack), "chr2")

  expect_identical(chromosome(dataTrack), "chr1")

  expect_identical(chromosome(seqTrack.dna), "chr1")
  chromosome(seqTrack.dna) <- "chr2"
  expect_identical(chromosome(seqTrack.dna), "chr2")
  
  expect_identical(chromosome(ideoTrack), "chrI")
  chromosome(ideoTrack) <- "chrII"
  expect_identical(chromosome(ideoTrack), "chrII")

  expect_identical(chromosome(highTrack), "chr1")
  chromosome(highTrack) <- "chr2"
  expect_identical(chromosome(highTrack), "chr2")

  expect_identical(chromosome(overTrack), "chr1")
  chromosome(overTrack) <- "chr2"
  expect_identical(chromosome(overTrack), "chr2")
  
})


test_that("genome accessors and replacement methods work", {
  expect_identical(is.na(genome(annoTrack)), TRUE)
  genome(annoTrack) <- "hg38"
  expect_identical(genome(annoTrack), "hg38")

  expect_identical(is.na(genome(seqTrack.dna)), TRUE)
  ## genome(seqTrack.dna) <- "hg38"
  ## expect_identical(genome(seqTrack.dna), "hg38")
  
  expect_identical(genome(ideoTrack), "sacCer3")
  ## genome(ideoTrack) <- "sacCer2"
  ## expect_identical(genome(ideoTrack), "sacCer2")
  
  expect_identical(genome(overTrack), NULL)
  genome(overTrack) <- "hg38"
  expect_identical(genome(overTrack), NULL)
})

test_that("strand accessors and replacement methods work", {
  expect_identical(strand(dataTrack), "*")
  strand(dataTrack) <- "+"
  expect_identical(strand(dataTrack), "+")
  expect_error(strand(dataTrack) <- rep("+", 2), "Invalid")  
  
  expect_identical(strand(annoTrack), "*")
  strand(annoTrack) <- "+"
  expect_identical(strand(annoTrack), "+")
  
  expect_identical(strand(axisTrack), "*")
  
  annoTrack <- AnnotationTrack(c(gr,gr2))
  expect_identical(strand(annoTrack), rep("*", 2))
  strand(annoTrack) <- rep("+", 2)
  expect_identical(strand(annoTrack), rep("+", 2))
  expect_error(strand(annoTrack) <- rep("+", 3), "Invalid replacement value or length of replacement")
})

test_that("names accessors and replacement methods work", {
  expect_identical(names(annoTrack), "AnnotationTrack")
  names(annoTrack) <- "AnnoTrack"
  expect_identical(names(annoTrack), "AnnoTrack")
})


test_that("gene, symbol, transcript, feature, exon and group accessors and replacement methods work", {
  expect_identical(gene(geneTrack), as.character(geneModels$gene))
  gene(geneTrack) <- paste0(as.character(geneModels$gene), ".1")
  expect_identical(gene(geneTrack), paste0(as.character(geneModels$gene), ".1"))
  
  expect_identical(symbol(geneTrack), as.character(geneModels$symbol))
  symbol(geneTrack) <- paste0(as.character(geneModels$symbol), ".1")
  expect_identical(symbol(geneTrack), paste0(as.character(geneModels$symbol), ".1"))
  
  expect_identical(transcript(geneTrack), as.character(geneModels$transcript))
  transcript(geneTrack) <- paste0(as.character(geneModels$transcript), ".1")
  expect_identical(transcript(geneTrack), paste0(as.character(geneModels$transcript), ".1"))
  
  expect_identical(exon(geneTrack), as.character(geneModels$exon))
  exon(geneTrack) <- paste0(as.character(geneModels$exon), ".1")
  expect_identical(exon(geneTrack), paste0(as.character(geneModels$exon), ".1"))
  
  expect_identical(feature(geneTrack), as.character(geneModels$feature))
  feature(geneTrack) <- paste0(as.character(geneModels$feature), ".1")
  expect_identical(feature(geneTrack), paste0(as.character(geneModels$feature), ".1"))
  
  expect_identical(group(geneTrack),paste0(as.character(geneModels$transcript), ".1"))
  group(geneTrack) <- paste0(as.character(geneModels$transcript), ".2")
  expect_identical(group(geneTrack), paste0(as.character(geneModels$transcript), ".2"))
  
})

