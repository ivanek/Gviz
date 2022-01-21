#' @include RangeTrack-class.R
NULL

## HighlightTrack Class ------------------------------------------------------

#' HighlightTrack class and methods
#'
#'
#' A container for other track objects from the Gviz package that allows for
#' the addition of a common highlighting area across tracks.
#'
#'
#' A track to conceptionally group other Gviz track objects into a meta track
#' for the sole purpose of overlaying all the contained tracks with the same
#' highlighting region as defined by the objects genomic ranges. During
#' rendering the contained tracks will be treated as if they had been provided
#' to the \code{plotTracks} function as individual objects.
#'
#' @template HighlightTrack-class_param
#'
#' @name HighlightTrack-class
#'
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{HighlightTrack}.
#' @section Objects from the Class:
#'
#' Objects can be created using the constructor function \code{HighlightTrack}.
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @examples
#' ## Object construction:
#' set.seed(123)
#' dat <- runif(100, min = -2, max = 22)
#' gt <- GenomeAxisTrack()
#' dt <- DataTrack(data = dat, start = sort(sample(200, 100)), width = 1, genome = "hg19")
#'
#' ht <- HighlightTrack(trackList = list(gt, dt))
#' @exportClass HighlightTrack
setClass("HighlightTrack",
    representation = representation(trackList = "list"),
    contains = c("RangeTrack"),
    prototype = prototype(dp = DisplayPars(
        col = "red",
        fill = "#FFE3E6",
        inBackground = TRUE
    ))
)

## Initialize ----------------------------------------------------------------

#' @describeIn HighlightTrack-class Initialize.
#' @export
setMethod("initialize", "HighlightTrack", function(.Object, trackList, ...) {
    .Object <- .updatePars(.Object, "HighlightTrack")
    .Object@trackList <- trackList
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

#' @describeIn HighlightTrack-class Constructor function for
#' `HighlightTrack-class`.
#' @export
HighlightTrack <- function(trackList = list(), range = NULL, start = NULL, end = NULL, width = NULL, chromosome, genome,
                           name = "HighlightTrack", ...) {
    ## Some defaults
    covars <- .getCovars(range)
    n <- max(c(length(start), length(end), length(width)), nrow(covars))
    ## Build a GRanges object from the inputs
    .missingToNull(c("strand", "chromosome", "genome"))
    args <- list(chromosome = chromosome, genome = genome)
    defs <- list(strand = "*", density = 1, chromosome = "chrNA", genome = NA)
    range <- .buildRange(
        range = range, start = start, end = end, width = width,
        args = args, defaults = defs, chromosome = chromosome, trackType = "HighlightTrack"
    )
    if (is.list(range)) {
        range <- GRanges()
    }
    if (!is.list(trackList)) {
        trackList <- list(trackList)
    }
    if (!all(vapply(trackList, is, class2 = "GdObject", FUN.VALUE = logical(1)))) {
        stop("All elements in 'trackList' must inherit from 'GdObject'")
    }
    ## If no chromosome was explicitly asked for we just take the first one in the GRanges object
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- if (length(range) > 0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    }
    ## And finally the object instantiation
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    return(new("HighlightTrack", trackList = trackList, chromosome = chromosome[1], range = range, name = name, genome = genome, ...))
}

## General accessors ---------------------------------------------------------

#' @describeIn HighlightTrack-class set display parameters using the values of
#' the named list in value. See \code{\link{settings}} for details on display
#' parameters and customization.
#' @export
setReplaceMethod("displayPars", signature("HighlightTrack", "list"), function(x, recursive = FALSE, value) {
    x <- setPar(x, value, interactive = FALSE)
    if (recursive) {
        x@trackList <- lapply(x@trackList, function(y) {
            displayPars(y) <- value
            return(y)
        })
    }
    return(x)
})

#' @describeIn HighlightTrack-class return the number of subtracks.
#' @export
setMethod("length", "HighlightTrack", function(x) length(x@trackList))

#' @describeIn HighlightTrack-class replace the value of the track's chromosome.
#' This has to be a valid UCSC chromosome identifier or an integer or character
#' scalar that can be reasonably coerced into one.
#' @export
setReplaceMethod("chromosome", "HighlightTrack", function(GdObject, value) {
    GdObject@trackList <- lapply(GdObject@trackList, function(x) {
        chromosome(x) <- value[1]
        x
    })
    GdObject@chromosome <- .chrName(value[1])
    return(GdObject)
})


## Annotation Accessors ------------------------------------------------------
## Stacking ------------------------------------------------------------------

#' @describeIn HighlightTrack-class Rrecompute the stacks based on the available
#' space and on the object's track items and stacking settings.
#' This really just calls the `setStacks` methods for the contained tracks and
#' only exists for dispatching reasons.
#' @export
setMethod("setStacks", "HighlightTrack", function(GdObject, ...) {
    GdObject@trackList <- lapply(GdObject@trackList, setStacks, ...)
    return(GdObject)
})

## Consolidate ---------------------------------------------------------------

#' @describeIn HighlightTrack-class Consolidate
#' For a `HighlightTrack` apply the method on each of the subtracks in
#' the `trackList` slot
# #' @keywords internal
#' @export
setMethod("consolidateTrack", signature(GdObject = "HighlightTrack"), function(GdObject, chromosome, ...) {
    GdObject@trackList <- lapply(GdObject@trackList, consolidateTrack, chromosome = chromosome, ...)
    return(GdObject)
})

## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------

#' @describeIn HighlightTrack-class subset all the contained tracks in an HighlightTrack by coordinates and sort if necessary.
#' @export
setMethod("subset", signature(x = "HighlightTrack"), function(x, ...) {
    x@trackList <- lapply(x@trackList, subset, ...)
    x <- callNextMethod()
    return(x)
})

## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------
## DrawAxis ------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------
## SetAs ---------------------------------------------------------------------
## Show ----------------------------------------------------------------------

#' @describeIn HighlightTrack-class Show method.
#' @export
setMethod("show", signature(object = "HighlightTrack"), function(object) {
    cat(sprintf(
        "HighlightTrack '%s' containing %i subtrack%s\n%s\n", names(object), length(object),
        ifelse(length(object) == 1, "", "s"), .annotationTrackInfo(object)
    ))
})
