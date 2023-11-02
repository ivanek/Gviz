## Gviz functions

test_that("checking the class and structure works", {
    expect_null(.checkClass(data.frame(), "data.frame"))
    expect_error(.checkClass(data.frame(), "matrix"), "of class")
    expect_error(.checkClass(class = "matrix"), "is missing with no default")
    expect_error(.checkClass(data.frame(), "data.frame", length = 2), "of length")
})

test_that("conversion of chromosome names works", {
    expect_identical(.chrName(character()), character())
    options("ucscChromosomeNames" = TRUE)
    expect_identical(.chrName(c("1", "chr1", "M", "MT", "X")), c("chr1", "chr1", "chrM", "chrM", "chrX"))
    expect_identical(.chrName("foo", force = TRUE), "chrfoo")
    expect_error(.chrName("foo"), "Please consider setting options\\(ucscChromosomeNames=FALSE\\)")
    options("ucscChromosomeNames" = FALSE)
    expect_identical(.chrName(c("1", "chr1", "M", "X")), c("1", "chr1", "M", "X"))
})


## not sure how to set up the test for the device resolution
## with opening PDF device?
test_that("finding location and pixel-size of current viewport works", {
    l <- list(
        location = c(x1 = 0, y1 = 0, x2 = 504, y2 = 504),
        size = c(width = 504, height = 504),
        ilocation = c(x1 = 0, y1 = 0, x2 = 7, y2 = 7),
        isize = c(width = 7, height = 7)
    )
    ## if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
    ## pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
    ## expect_identical(vpLocation(), l)
    expect_equal(vpLocation(), l)
    ## dev.off()
    ## unlink("Rplots.pdf")
})


test_that("estimating of device resolution works", {
    ## if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
    ## pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
    ## expect_identical(devRes(), c(xres=72, yres=72))
    expect_equal(devRes(), c(xres = 72, yres = 72))
    ## dev.off()
    ## unlink("Rplots.pdf")

    ## if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
    ## pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
    pushViewport(viewport())
    ## expect_identical(devRes(), c(xres=72, yres=72))
    expect_equal(devRes(), c(xres = 72, yres = 72))
    popViewport()
    ## dev.off()
    ## unlink("Rplots.pdf")
})

test_that("estimating of coordinates for an HTML image map works", {
    m <- matrix(rep(0:1, each = 4), ncol = 4, byrow = TRUE)
    mm <- matrix(rep(as.numeric(c(0, 1, 504, 503)), 2),
        ncol = 4,
        dimnames = list(NULL, c("x1", "y1", "x2", "y2"))
    )
    mm <- as.data.frame(mm)
    ## if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
    ## pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
    ## expect_identical(.getImageMap(m), mm)
    expect_equal(.getImageMap(m), mm)
    ## dev.off()
    ## unlink("Rplots.pdf")
})

test_that("deep copying of displaypars works", {
    expect_identical(.needsStacking(annoTrack), TRUE)
    expect_identical(.deepCopyPars(annoTrack), annoTrack)
    expect_identical(displayPars(.deepCopyPars(annoTrack)), displayPars(annoTrack))
})

test_that("check of stacking works", {
    expect_identical(.needsStacking(annoTrack), TRUE)
})

test_that("checking of the strand works", {
    displayPars(dataTrack)$reverseStrand <- TRUE
    displayPars(overTrack@trackList[[1]])$reverseStrand <- TRUE
    expect_identical(.whichStrand(annoTrack), "forward")
    expect_identical(
        .whichStrand(list(dataTrack, annoTrack, highTrack, overTrack)),
        c("reverse", "forward", "forward", "reverse", "forward")
    )
})

test_that("import of sequence from FASTA file works", {
    expect_identical(.import.fasta(fastafile, GRanges("chr3", IRanges(1, 4))), DNAStringSet())
    expect_identical(.import.fasta(fastafile, GRanges("chr1", IRanges(1, 4))), DNAStringSet(c("chr1" = "CTAN")))
})

test_that("import of sequence from 2bit file works", {
    twobitfile <- system.file("extdata/test.2bit", package = "Gviz")
    expect_identical(.import.2bit(twobitfile, GRanges("chr3", IRanges(1, 4))), DNAStringSet())
    expect_identical(.import.2bit(twobitfile, GRanges("chr1", IRanges(1, 4))), DNAStringSet(c("chr1" = "CTAN")))
})

test_that("import of alignments from BAM file works", {
    empty <- GRanges()
    mcols(empty) <- DataFrame(
        id = character(), cigar = character(), mapq = integer(), flag = integer(),
        md = character(), seq = DNAStringSet(), isize = integer(), groupid = integer(),
        status = factor(levels = c("mated", "ambiguous", "unmated"))
    )

    unlink(paste(bamfile, "bai", sep = "."))
    expect_error(.import.bam.alignments(bamfile, GRanges("chr1", IRanges(189891401, 189894000))), "Unable to find index for BAM file")
    bamfile <- Rsamtools::asBam(samfile, indexDestination = FALSE, overwrite = TRUE)
    Rsamtools::indexBam(bamfile)

    expect_identical(.import.bam.alignments(bamfile, GRanges("chr1", IRanges(189891401, 189894000))), bamgr)
    expect_identical(.import.bam.alignments(bamfile, GRanges("chr2", IRanges(1, 2))), empty)
    seqlevels(empty) <- "chr1"
    expect_identical(.import.bam.alignments(bamfile, GRanges("chr1", IRanges(1, 2))), empty)
})

test_that("estimate of required verticalSpace works", {
    expect_identical(.verticalSpace(dataTrack, 10), 5)
    expect_identical(.verticalSpace(axisTrack, 10), 1)
    expect_identical(.verticalSpace(ideoTrack, 10), 1)
    expect_identical(.verticalSpace(seqTrack.dna, 10), 1)
    expect_identical(.verticalSpace(setStacks(annoTrack), 10), 1) ## require stacks
    expect_identical(.verticalSpace(alnTrack, 10), 1)
})



test_that("conversion of ranges to summarizedJunctions works", {
    ## with + strand defined in readStrand column
    range <- GRanges("chr1", IRanges(1, 10), cigar = "1M8N1M", readStrand = Rle(factor("+", levels = c("+", "-", "*"))), entityId = 1)
    juns <- GRanges("chr1", IRanges(start = 2, end = 9), score = 1L, plus_score = 1L, minus_score = 0L)
    expect_identical(.create.summarizedJunctions.for.sashimi.junctions(range), juns)
    ## without strand definition
    range <- GRanges("chr1", IRanges(1, 10), cigar = "1M8N1M", entityId = 1)
    juns <- GRanges("chr1", IRanges(start = 2, end = 9), score = 1L, plus_score = 0L, minus_score = 0L)
    expect_identical(.create.summarizedJunctions.for.sashimi.junctions(range), juns)
})


test_that("conversion of junction to list for plotting works", {
    juns <- GRanges("chr1", IRanges(start = 2, end = 9), score = 1L, plus_score = 0L, minus_score = 0L)
    out <- list(
        x = c(2, 5, 9),
        y = c(0, 1, 0),
        id = c(1L, 1L, 1L),
        score = 1, scaled = 10
    )
    filt <- GRanges("chr1", IRanges(start = 2, end = 9))
    out2 <- list(
        x = c(3, 5, 8),
        y = c(0, 1, 0),
        id = c(1L, 1L, 1L),
        score = 1, scaled = 10
    )
    filt2 <- GRanges("chr1", IRanges(start = 3, end = 8))
    expect_identical(.convert.summarizedJunctions.to.sashimi.junctions(juns), out)
    ## if minimum score filter works
    expect_identical(
        .convert.summarizedJunctions.to.sashimi.junctions(juns, score = 2L),
        list(x = numeric(), y = numeric(), id = integer(), score = numeric(), scaled = numeric())
    )
    ## filter match
    expect_identical(.convert.summarizedJunctions.to.sashimi.junctions(juns, filter = filt), out)
    ## filter match (none)
    expect_identical(
        .convert.summarizedJunctions.to.sashimi.junctions(juns, filter = filt2),
        list(x = numeric(), y = numeric(), id = integer(), score = numeric(), scaled = numeric())
    )
    ## filter match (with filterTolerance)
    expect_identical(.convert.summarizedJunctions.to.sashimi.junctions(juns, filter = filt2, filterTolerance = 1), out2)
    ## negative filterTolerance
    expect_warning(
        .convert.summarizedJunctions.to.sashimi.junctions(juns, filter = filt, filterTolerance = -1),
        "can't be negative, taking absolute value of it"
    )
    ## transformation
    expect_identical(.convert.summarizedJunctions.to.sashimi.junctions(juns, trans = list(function(x) {
        x
    })), out)
    expect_error(.convert.summarizedJunctions.to.sashimi.junctions(juns, trans = 1), "must be a function with a single argument")
    expect_error(.convert.summarizedJunctions.to.sashimi.junctions(juns, trans = function(x) {
        x[-1]
    }), "invalid output")
})


## test_that("only first warning works", {
##   warn_once({
##     for (i in 1:10) { warn_if_first("foo", "oh, no! foo!") }
##     for (i in 1:10) { warn_if_first("bar", "oh, no! bar!") }
##     sapply(1:10, function(x) {
##       warn_if_first("foo", "oh, no! foo again! (not really)")
##       warn_if_first("foobar", "foobar, too!")
##     })
##     "DONE!"
##   })
## })
