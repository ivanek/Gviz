library(Gviz)

context("Gviz functions")

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

context("Gviz functions, graphical device")
## not sure how to set up the test for the device resolution
## with opening PDF device?
test_that("finding location and pixel-size of current viewport works", {
  l <- list(location = c(x1 = 0, y1 = 0, x2 = 504, y2 = 504), 
            size = c(width = 504, height = 504), 
            ilocation = c(x1 = 0, y1 = 0, x2 = 7, y2 = 7), 
            isize = c(width = 7, height = 7))
  #if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
  #pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
  expect_identical(Gviz:::vpLocation(), l)
  #dev.off()
  #unlink("Rplots.pdf")
})


test_that("estimating of device resolution works", {
  #if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
  #pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
  expect_identical(Gviz:::devRes(), c(xres=72, yres=72))
  #dev.off()
  #unlink("Rplots.pdf")
  
  #if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
  #pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
  pushViewport(viewport())
  expect_identical(Gviz:::devRes(), c(xres=72, yres=72))
  popViewport()
  #dev.off()
  #unlink("Rplots.pdf")
})

test_that("estimating of coordinates for an HTML image map works", {
  m <- matrix(rep(0:1, each=4), ncol=4, byrow=TRUE)
  mm <- matrix(rep(as.numeric(c(0,1,504,503)), 2), ncol=4,
               dimnames=list(NULL, c("x1","y1","x2","y2")))
  mm <- as.data.frame(mm)
  # if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
  # pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
  expect_identical(Gviz:::.getImageMap(m), mm)
  # dev.off()
  # unlink("Rplots.pdf")
})

context("Gviz functions, helper functions")

test_that("conversion of ranges to summarizedJunctions works", {
  ## with + strand defined in readStrand column
  range <- GRanges("chr1", IRanges(1,10), cigar="1M8N1M", readStrand=Rle(factor("+", levels=c("+","-","*"))), entityId=1)
  juns <- GRanges("chr1", IRanges(start=2,end=9), score=1L, plus_score=1L, minus_score=0L)
  expect_identical(Gviz:::.create.summarizedJunctions.for.sashimi.junctions(range), juns)
  ## without strand definition
  range <- GRanges("chr1", IRanges(1,10), cigar="1M8N1M", entityId=1)
  juns <- GRanges("chr1", IRanges(start=2,end=9), score=1L, plus_score=0L, minus_score=0L)
  expect_identical(Gviz:::.create.summarizedJunctions.for.sashimi.junctions(range), juns)
})


test_that("conversion of junction to list for plotting works", {
  #<- function(juns, score=1L, lwd.max=10, strand="*", filter=NULL, filterTolerance=0L) {
  juns <- GRanges("chr1", IRanges(start=2,end=9), score=1L, plus_score=0L, minus_score=0L)
  out <- list(x = c(2, 5, 9), 
              y = c(0, 1, 0), 
              id = c(1L, 1L, 1L), 
              score = 1, scaled = 10)
  expect_identical(Gviz:::.convert.summarizedJunctions.to.sashimi.junctions(juns), out)
  # if minimum score filter works
  expect_identical(Gviz:::.convert.summarizedJunctions.to.sashimi.junctions(juns, score=2L), 
                   list(x=numeric(), y=numeric(), id=integer(), score=numeric(), scaled=numeric()))
})
