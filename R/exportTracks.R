## Write all tracks in a list of tracks into a single BED file.
# FIXME: Need to support WIG exports as well...

#' Export Gviz tracks into an annotation file representation.
#'
#' This function is still a bit experimental. Write all tracks provided as a
#' list `tracks` into a single BED file. So far only BED export is supported.
#'
#' @name exportTracks
#'
#' @param tracks A list of annotation track objects to be exported into a
#' single BED file.
#' @param range A numeric vector or length 2. The genomic range to display when
#' opening the file in a browser.
#' @param chromosome The chromosome to display when opening the file in a browser.
#' @param file Character, the path to the file to write into.
#'
#' @return The function is called for its side effect of writing to a file.
#'
#' @author Florian Hahne
#'
#' @examples
#' ## export AnnotationTrack to BED file
#' at <- AnnotationTrack(start = seq(1, 10), width = 1, chromosome = "chr1")
#' exportTracks(list(at),
#'     range = c(1, 10), chromosome = "chr1",
#'     file = paste0(tempfile(), ".bed")
#' )
#' @importClassesFrom rtracklayer UCSCData
#' @importFrom utils write.table
#' @importFrom grDevices col2rgb
#'
#' @export exportTracks
exportTracks <- function(tracks, range, chromosome, file) {
    if (missing(file)) {
        file <- "customTracks.bed"
    }
    con <- file(file, open = "wt")
    writeLines(
        sprintf("browser position %s:%i-%i", chromosome, range[1], range[2]),
        con
    )
    writeLines("browser hide all", con)
    for (t in seq_along(tracks))
    {
        track <- tracks[[t]]
        if (length(track) > 0 && (is(track, "AnnotationTrack") || is(track, "GeneRegion"))) {
            track <- as(track, "UCSCData")
            writeLines(as(track@trackLine, "character"), con)
            ## nextMet <- selectMethod("export.bed", c("RangedData", "characterORconnection"))
            ## nextMet(as(track, "GRanhes"), con)
            .expBed(as(track, "GRanges"), con)
        }
    }
    close(con)
}

## This function is broken in the rtracklayer package
.expBed <- function(object, con, variant = c("base", "bedGraph", "bed15"), color, append) {
    variant <- match.arg(variant)
    name <- strand <- thickStart <- thickEnd <- color <- NULL
    blockCount <- blockSizes <- blockStarts <- NULL
    df <- data.frame(as.character(seqnames(object)), start(object) - 1, end(object))
    score <- score(object)
    if (!is.null(score)) {
        if (!is.numeric(score) || any(is.na(score))) {
            stop("Scores must be non-NA numeric values")
        }
    }
    if (variant == "bedGraph") {
        if (is.null(score)) {
            score <- 0
        }
        df$score <- score
    } else {
        blockSizes <- object$blockSizes
        blockStarts <- object$blockStarts
        if (variant == "bed15" && is.null(blockSizes)) {
            blockStarts <- blockSizes <- ""
        }
        if (!is.null(blockSizes) || !is.null(blockStarts)) {
            if (is.null(blockSizes)) {
                stop("'blockStarts' specified without 'blockSizes'")
            }
            if (is.null(blockStarts)) {
                stop("'blockSizes' specified without 'blockStarts'")
            }
            lastBlock <- function(x) sub(".*,", "", x)
            lastSize <- lastBlock(blockSizes)
            lastStart <- lastBlock(blockStarts)
            if (any(df[[2]] + as.integer(lastSize) + as.integer(lastStart) != df[[3]]) ||
                any(sub(",.*", "", blockStarts) != 0)) {
                stop("blocks must span entire feature")
            }
            blockCount <- vapply(strsplit(blockSizes, ","), length, FUN.VALUE = numeric(1L))
        }
        if (is.null(color)) {
            color <- object$itemRgb
        }
        if (is.null(color) && !is.null(blockCount)) {
            color <- "0"
        } else if (!is.null(color)) {
            nacol <- is.na(color)
            colmat <- col2rgb(color)
            color <- paste(colmat[1, ], colmat[2, ], colmat[3, ], sep = ",")
            color[nacol] <- "0"
        }
        thickStart <- object$thickStart
        thickEnd <- object$thickEnd
        if (is.null(thickStart) && !is.null(color)) {
            thickStart <- start(object)
            thickEnd <- end(object)
        }
        strand <- object$strand
        if (!is.null(thickStart) && is.null(strand)) {
            strand <- rep(NA, length(object))
        }
        if (!is.null(strand) && is.null(score)) {
            score <- 0
        }
        name <- object$name
        if (is.null(name)) {
            name <- rownames(object)
        }
        if (!is.null(score) && is.null(name)) {
            name <- rep(NA, length(object))
        }
        df$name <- name
        df$score <- score
        df$strand <- strand
        df$thickStart <- thickStart
        df$thickEnd <- thickEnd
        df$itemRgb <- color
        df$blockCount <- blockCount
        df$blockSizes <- blockSizes
        df$blockStarts <- blockStarts
        if (variant == "bed15") {
            df$expCount <- object$expCount
            df$expIds <- object$expIds
            df$expScores <- object$expScores
        }
    }
    scipen <- getOption("scipen")
    options(scipen = 100)
    on.exit(options(scipen = scipen))
    write.table(df, con,
        sep = "\t", col.names = FALSE, row.names = FALSE,
        quote = FALSE, na = ".", append = append
    )
}
