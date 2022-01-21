#' @include RangeTrack-class.R
NULL

## OverlayTrack Class --------------------------------------------------------

#' OverlayTrack class and methods
#'
#'
#' A container for other track objects from the Gviz package that allows for
#' overlays of their content on the same region of the plot.
#'
#'
#' A track to conceptionally group other Gviz track objects into a meta track
#' in order to merge them into a single overlay visualization. Only the first
#' track in the supplied list will be inferred when setting up the track title
#' and axis, for all the other tracks only the panel content is plotted.
#'
#' @name OverlayTrack-class
#'
#' @param trackList A list of Gviz track objects that all have to inherit from
#' class \code{GdObject}.
#' @param name Character scalar of the track's name. This is not really used
#' and only exists fro completeness.
#' @param \dots All additional parameters are ignored.
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{OverlayTrack}.
#' @section Objects from the Class:
#'
#' Objects can be created using the constructor function \code{OverlayTrack}.
#'
#' @author Florian Hahne
#' @inherit GdObject-class seealso
#' @examples
#'
#'
#' ## Object construction:
#' set.seed(123)
#' dat <- runif(100, min = -2, max = 22)
#' dt1 <- DataTrack(data = dat, start = sort(sample(200, 100)), width = 1, genome = "hg19")
#' dt2 <- DataTrack(data = dat, start = sort(sample(200, 100)), width = 1, genome = "hg19")
#' ot <- OverlayTrack(trackList = list(dt1, dt2))
#' @exportClass OverlayTrack
setClass("OverlayTrack",
    representation = representation(trackList = "list"),
    contains = c("GdObject"),
    prototype = prototype(dp = DisplayPars())
)

## Initialize ----------------------------------------------------------------

#' @describeIn OverlayTrack-class Initialize.
#' @export
setMethod("initialize", "OverlayTrack", function(.Object, trackList, ...) {
    .Object <- .updatePars(.Object, "OverlayTrack")
    .Object@trackList <- trackList
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

#' @describeIn OverlayTrack-class Constructor function for `OverlayTrack-class`.
#' @export
OverlayTrack <- function(trackList = list(), name = "OverlayTrack", ...) {
    if (!is.list(trackList)) {
        trackList <- list(trackList)
    }
    if (!all(vapply(trackList, is, class2 = "GdObject", FUN.VALUE = logical(1)))) {
        stop("All elements in 'trackList' must inherit from 'GdObject'")
    }
    return(new("OverlayTrack", trackList = trackList, name = name, ...))
}

## General accessors ---------------------------------------------------------

#' @describeIn OverlayTrack-class set display parameters using the values of
#' the named list in value. See \code{\link{settings}} for details on
#' display parameters and customization.
#' @export
setReplaceMethod("displayPars", signature("OverlayTrack", "list"), function(x, recursive = FALSE, value) {
    x <- setPar(x, value, interactive = FALSE)
    if (recursive) {
        x@trackList <- lapply(x@trackList, function(y) {
            displayPars(y) <- value
            return(y)
        })
    }
    return(x)
})

#' @describeIn OverlayTrack-class return the number of subtracks.
#' @export
setMethod("length", "OverlayTrack", function(x) length(x@trackList))

#' @describeIn OverlayTrack-class return the chromosome for which the track
#' is defined.
#' @export
setMethod("chromosome", "OverlayTrack", function(GdObject) {
    unique(unlist(lapply(GdObject@trackList, chromosome)))
})

#' @describeIn OverlayTrack-class replace the value of the track's chromosome.
#' This has to be a valid UCSC chromosome identifier or an integer or character
#' scalar that can be reasonably coerced into one.
#' @export
setReplaceMethod("chromosome", "OverlayTrack", function(GdObject, value) {
    GdObject@trackList <- lapply(GdObject@trackList, function(x) {
        chromosome(x) <- value
        x
    })
    return(GdObject)
})


## Annotation Accessors ------------------------------------------------------
## Stacking ------------------------------------------------------------------

#' @describeIn OverlayTrack-class recompute the stacks based on the available
#' space and on the object's track items and stacking settings.
#' This really just calls the setStacks methods for the contained
#' tracks and only exists for dispatching reasons.
#' @export
setMethod("setStacks", "OverlayTrack", function(GdObject, ...) {
    GdObject@trackList <- lapply(GdObject@trackList, setStacks, ...)
    return(GdObject)
})

## Consolidate ---------------------------------------------------------------

#' @describeIn OverlayTrack-class #' For a `OverlayTrack` apply the method on
#' each of the subtracks in the `trackList` slot
#' @keywords internal
#' @export
setMethod("consolidateTrack", signature(GdObject = "OverlayTrack"), function(GdObject, chromosome, ...) {
    GdObject@trackList <- lapply(GdObject@trackList, consolidateTrack, chromosome = chromosome, ...)
    return(GdObject)
})

## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------

#' @describeIn OverlayTrack-class plot subset all the contained tracks in an
#' `OverlayTrack` by coordinates and sort if necessary.
#' @export
setMethod("subset", signature(x = "OverlayTrack"), function(x, ...) {
    x@trackList <- lapply(x@trackList, subset, ...)
    return(x)
})

## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------
## DrawAxis ------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------

#' @describeIn OverlayTrack-class plot the object to a graphics device.
#' The return value of this method is the input object, potentially updated
#' during the plotting operation. Internally, there are two modes in which the
#' method can be called. Either in 'prepare' mode, in which case no plotting is
#' done but the object is preprocessed based on the available space, or in
#' 'plotting' mode, in which case the actual graphical output is created.
#' Since subsetting of the object can be potentially costly, this can be
#' switched off in case subsetting has already been performed before or
#' is not necessary.
#'
#' @export
setMethod("drawGD", signature("OverlayTrack"), function(GdObject, ...) {
    GdObject@trackList <- lapply(GdObject@trackList, drawGD, ...)
    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------
## Show ----------------------------------------------------------------------

#' @describeIn OverlayTrack-class Show method.
#' @export
setMethod("show", signature(object = "OverlayTrack"), function(object) {
    cat(sprintf(
        "OverlayTrack '%s' containing %i subtrack%s\n", names(object), length(object),
        ifelse(length(object) == 1, "", "s")
    ))
})
