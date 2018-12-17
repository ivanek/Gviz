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

## not sure how to set up the test for the device resolution
## with opening PDF device?
test_that("estimating of device resolution works", {
  #if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
  #pdf(file="Rplots.pdf", width=7, height=5, pointsize=12)
  expect_identical(Gviz:::devRes(), c(xres=72, yres=72))
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

