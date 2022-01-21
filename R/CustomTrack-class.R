#' @include GdObject-class.R
NULL

#' CustomTrack class and methods
#'
#' A fully customizable track object to be populated via a user-defined
#' plotting function.
#'
#' A track to allow for any sort of plotting, with the currently displayed
#' genomic location set. Essentially this acts as a simple callback into the
#' `Gviz` plotting machinery after all the track panels and coordinates
#' have been set up. It is entirely up to the user what to plot in the track,
#' or even to use the predefined coordinate system. The only prerequisite is
#' that all plotting operations need to utilize Grid graphics.
#'
#' @template CustomTrack-class_param
#'
#' @name CustomTrack-class
#'
#' @return
#' The return value of the constructor function is a new object of class
#' `CustomTrack`.
#'
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#'
#' @examples
#' ## Object construction:
#'
#' ## An empty object
#' CustomTrack()
#' @export
setClass("CustomTrack",
    contains = c("GdObject"),
    representation = representation(
        plottingFunction = "function",
        variables = "list"
    ),
    prototype = prototype(dp = DisplayPars())
)

## Initialize ----------------------------------------------------------------

#' @describeIn CustomTrack-class Initialize.
#' @export
setMethod("initialize", "CustomTrack", function(.Object, plottingFunction, variables, ...) {
    .Object <- .updatePars(.Object, "CustomTrack")
    .Object@plottingFunction <- plottingFunction
    .Object@variables <- variables
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

#' @param plottingFunction A user-defined function to be executed once the
#' track coordinates have been properly set up. The function needs to accept
#' two mandatory arguments: `GdObject`, the `CustomTrack` object to
#' be plotted, and `prepare`, a logical flag indicating whether the
#' function has been called in preparation mode or in drawing mode. It also
#' needs to return the input `GdObject`, potentially with modifications.
#' @param variables A list of additional variables for the user-defined
#' plotting function.
#' @param name Character scalar of the track's name.
#' @param \dots Additional items which will all be interpreted as further
#' display parameters. See `settings` and the "Display Parameters"
#' section below for details.
#' @describeIn CustomTrack-class Objects can be created using the constructor function.
#' @export
CustomTrack <- function(plottingFunction = function(GdObject, prepare = FALSE, ...) {}, variables = list(), name = "CustomTrack", ...) {
    return(new("CustomTrack", plottingFunction = plottingFunction, variables = variables, name = name, ...))
}

## General accessors ---------------------------------------------------------
## Annotation Accessors ------------------------------------------------------
## Stacking ------------------------------------------------------------------
## Consolidate ---------------------------------------------------------------
## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------
## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------
## DrawAxis ------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------

#' @describeIn CustomTrack-class plot the object to a graphics device.
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
setMethod("drawGD", signature("CustomTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, ...) {
    rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
    xscale <- if (!rev) c(minBase, maxBase) else c(maxBase, minBase)
    pushViewport(viewport(xscale = xscale, clip = TRUE))
    tmp <- GdObject@plottingFunction(GdObject, prepare = prepare)
    if (!is(tmp, "CustomTrack")) {
        warning("The plotting function of a CustomTrack has to return the input object. Using the original CustomTrack object now.")
    } else {
        GdObject <- tmp
    }
    popViewport(1)
    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------
## Show ----------------------------------------------------------------------

#' @describeIn CustomTrack-class Show method.
#' @export
setMethod("show", signature(object = "CustomTrack"), function(object) {
    cat(sprintf("CustomTrack '%s'\n", names(object)))
})
