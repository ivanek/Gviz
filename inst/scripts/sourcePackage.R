##library(methods)
library(grid)
library(IRanges)
library(rtracklayer)
library(GenomicRanges)
library(lattice)
library(biomaRt)
library(RColorBrewer)
library(Biobase)
library(grid)
library(AnnotationDbi)

## options(error=recover)
removeAndSource <- function(file, remove=FALSE){
    if(remove){
        exp <- parse(file)
        sapply(exp, function(x){
            x <- as.character(x)
            if(x[[1]]=="setMethod")
                suppressWarnings(try(removeMethod(x[[2]], ifelse(length(grep("\\(", x[3]))>0, eval(parse(text=x[3])), x[3]))))
            if(x[[1]]=="setReplaceMethod")
                suppressWarnings(try(removeMethod(paste(x[[2]], "<-", sep=""), ifelse(length(grep("\\(", x[3]))>0, eval(parse(text=x[3])), x[3]))))
            if(x[[1]]=="setAs")
                suppressWarnings(try(removeMethod("coerce", x[2:3])))
        })
    }
    source(file)
}

path <- "~/Rpacks/Gviz/R"
files <- c("Gviz.R", "AllGenerics.R", "AllClasses.R", "Gviz-methods.R")
sapply(file.path(path, files), removeAndSource)


## dtTrack <- DataTrack(start=seq(1,1000, len=100), width=10, data=matrix(runif(400), nrow=4), chromosome=1, genome="mm9", name="random data")
