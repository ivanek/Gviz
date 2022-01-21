#' @include DisplayPars-class.R
NULL

###  GdObject-class ----------------------------------------------------------

#' GdObject class and methods
#'
#' The virtual parent class for all track items in the Gviz package.  This
#' class definition contains all the common entities that are needed for a
#' track to be plotted. During object instantiation for any of the sub-classes
#' inheriting from `GdObject`, this class' global initializer has to be
#' called in order to assure that all necessary settings are present.
#'
#' @template DisplayParams_list
#'
#' @template ImageMap-class_param
#' @template GdObject-class_param
#'
#' @template GdObject-class_slot
#'
#' @name GdObject-class
#'
#' @return A virtual class: No objects may be created from it.
#'
#' @author Florian Hahne
#'
#' @seealso
#' \code{\linkS4class{DisplayPars}}
#'
#' \code{\linkS4class{GdObject}}
#'
#' \code{\linkS4class{GRanges}}
#'
#' \code{\linkS4class{HighlightTrack}}
#'
#' \code{\linkS4class{ImageMap}}
#'
#' \code{\linkS4class{IRanges}}
#'
#' \code{\linkS4class{RangeTrack}}
#'
#' \code{\linkS4class{DataTrack}}
#'
#' \code{\link{collapsing}}
#'
#' \code{\link{grouping}}
#'
#' \code{\link{panel.grid}}
#'
#' \code{\link{plotTracks}}
#'
#' \code{\link{settings}}
#'
#' @examples
#' ## This is a reference class therefore we show below
#' ## an example from AnnotationTrack:
#'
#' ## An empty object
#' AnnotationTrack()
#'
#' ## Construct from individual arguments
#' st <- c(2000000, 2070000, 2100000, 2160000)
#' ed <- c(2050000, 2130000, 2150000, 2170000)
#' str <- c("-", "+", "-", "-")
#' gr <- c("Group1", "Group2", "Group1", "Group3")
#'
#' annTrack <- AnnotationTrack(
#'     start = st, end = ed, strand = str, chromosome = 7,
#'     genome = "hg19", feature = "test", group = gr,
#'     id = paste("annTrack item", 1:4),
#'     name = "generic annotation", stacking = "squish"
#' )
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(annTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#' ## Plotting
#' plotTracks(annTrack)
#' @exportClass GdObject
setClass("GdObject",
    representation = representation("VIRTUAL",
        dp = "DisplayPars",
        name = "character",
        imageMap = "ImageMapOrNULL"
    ),
    prototype = prototype(
        dp = DisplayPars(
            alpha = 1,
            alpha.title = NULL,
            background.panel = "transparent",
            background.title = "lightgray",
            background.legend = "transparent",
            cex.axis = NULL,
            cex.title = NULL,
            cex = 1,
            col.axis = "white",
            col.border.title = "white",
            col.frame = "lightgray",
            col.grid = .DEFAULT_SHADED_COL,
            col.line = NULL,
            col.symbol = NULL,
            col.title = "white",
            col = .DEFAULT_SYMBOL_COL,
            collapse = TRUE,
            fill = .DEFAULT_FILL_COL,
            fontcolour = "black",
            fontface.title = 2,
            fontface = 1,
            fontfamily.title = "sans",
            fontfamily = "sans",
            fontsize = 12,
            frame = FALSE,
            grid = FALSE,
            h = -1,
            lineheight = 1,
            lty.grid = "solid",
            lty = "solid",
            lwd.border.title = 1,
            lwd.title = 1,
            lwd.grid = 1,
            lwd = 1,
            min.distance = 1,
            min.height = 3,
            min.width = 1,
            reverseStrand = FALSE,
            rotation.title = 90,
            rotation = 0,
            showAxis = TRUE,
            showTitle = TRUE,
            size = 1,
            v = -1
        ),
        name = "GdObject",
        imageMap = NULL
    )
)

##  GdObject Constructor -----------------------------------------------------

## We add everything that hasn't been clobbered up so far as
## additional DisplayParameters.  Also, the `dp` slot must be
## re-initiated here in order to get a fresh environment for each
## instance of the class.

#' @describeIn GdObject-class Initialize the object. This involves setting up a
#' new environment for the display parameters and filling it up with the current
#' settings. All arguments that have not been clobbered up by one of the
#' sub-class initializers are considered to be additional display parameters
#' and are also added to the environment. See `settings` for details on setting
#' graphical parameters for tracks.
#' @keywords internal
setMethod("initialize", "GdObject", function(.Object, name, ...) {
    ## update the default parameters first
    .makeParMapping()
    .Object <- .updatePars(.Object, "GdObject")
    ## now rebuild the slot to get a new environment
    pars <- getPar(.Object, hideInternal = FALSE)
    .Object@dp <- DisplayPars()
    .Object <- setPar(.Object, pars, interactive = FALSE)
    if (!missing(name)) {
        .Object@name <- if (is.null(name)) "" else name
    }
    ## Finally clobber up everything that's left
    .Object <- setPar(.Object, list(...), interactive = FALSE)
    return(.Object)
})

##  GdObject Methods Setters -------------------------------------------------

## We need to set and query DisplayPars in the the initializer, hence
## the appropriate methods have to be defined here first.
#' @describeIn GdObject-class set the single display parameter name to value.
#' Note that display parameters in the `GdObject-class` are pass-by-reference,
#' so no re-assignment to the symbol `obj` is necessary. See settings for
#' details on display parameters and customization.
setMethod("setPar", signature("GdObject", "character"), function(x, name, value, interactive = TRUE) {
    newDp <- setPar(x@dp, name, value, interactive = interactive)
    x@dp <- newDp
    return(x)
})

#' @describeIn GdObject-class set display parameters by the values of the named
#' list in value. Note that display parameters in the `GdObject-class` are
#' pass-by-reference, so no re-assignment to the symbol `obj` is necessary.
#' See settings for details on display parameters and customization.
setMethod("setPar", signature("GdObject", "list"), function(x, value, interactive = TRUE) {
    newDp <- setPar(x@dp, value, interactive = interactive)
    x@dp <- newDp
    return(x)
})

#' @describeIn GdObject-class set display parameters using the values of the
#' named list in `value`. See `settings` for details on display parameters
#' and customization.
setReplaceMethod("displayPars", signature("GdObject", "list"), function(x, recursive = FALSE, value) {
    x <- setPar(x, value, interactive = FALSE)
    return(x)
})

##  GdObject Methods Getters -------------------------------------------------

#' @describeIn GdObject-class alias for the `displayPars` method.
#' See `settings` for details on display parameters and customization.
setMethod("getPar", c("GdObject", "character"), function(x, name, asIs = FALSE) getPar(x@dp, name, asIs = asIs))

#' @describeIn GdObject-class alias for the `displayPars` method.
#' See `settings` for details on display parameters and customization.
setMethod("getPar", c("GdObject", "missing"), function(x, hideInternal = TRUE) getPar(x@dp, hideInternal = hideInternal))

#' @describeIn GdObject-class list the value of the display parameter name.
#' See `settings` for details on display parameters and customization.
setMethod("displayPars", c("GdObject", "character"), function(x, name) getPar(x, name))

#' @describeIn GdObject-class list the value of all available display parameters.
#' See `settings` for details on display parameters and customization.
setMethod("displayPars", c("GdObject", "missing"), function(x, hideInternal = TRUE) getPar(x, hideInternal = hideInternal))

##  GdObject Methods Coord and Tags  -----------------------------------------

#' @describeIn GdObject-class return the coordinates from the internal image map.
setMethod("coords", "GdObject", function(ImageMap) coords(imageMap(ImageMap)))

#' @describeIn GdObject-class return the tags from the internal image map.
setMethod("tags", "GdObject", function(ImageMap) tags(imageMap(ImageMap)))

##  GdObject Methods Subset --------------------------------------------------

#' @describeIn GdObject-class subset a `GdObject` by coordinates.
#' Most of the respective sub-classes inheriting from `GdObject` overwrite this
#' method, the default is to return the unaltered input object.
setMethod("subset", signature(x = "GdObject"), function(x, ...) x)

##  GdObject Methods Name ----------------------------------------------------

#' @describeIn GdObject-class return the value of the `name` slot.
setMethod("names", "GdObject", function(x) x@name)

#' @describeIn GdObject-class set the value of the `name` slot.
setReplaceMethod(
    "names", signature("GdObject", "character"),
    function(x, value) {
        x@name <- value[1]
        return(x)
    }
)

##  GdObject Methods Group ---------------------------------------------------

#' @describeIn GdObject-class Generics for `group`.
#' @exportMethod group
#' @keywords internal
setGeneric("group", function(GdObject, ...) standardGeneric("group"))

#' @exportMethod "group<-"
#' @describeIn GdObject-class Generics for `group<-`.
#' @keywords internal
setGeneric("group<-", function(GdObject, value) standardGeneric("group<-"))

#' @describeIn GdObject-class return grouping information for the individual
#' items in the track. Unless overwritten in one of the sub-classes,
#' this usually returns `NULL`.
setMethod("group", "GdObject", function(GdObject) NULL)

##  GdObject Methods ImageMap ------------------------------------------------

#' @exportMethod imageMap
#' @describeIn GdObject-class Generics for `imageMap`.
#' @keywords internal
setGeneric("imageMap", function(GdObject, ...) standardGeneric("imageMap"))

#' @describeIn  GdObject-class Extract the content of the `imageMap` slot.
#' @keywords internal
#' @export
setMethod("imageMap", "GdObject", function(GdObject) GdObject@imageMap)

#' @exportMethod "imageMap<-"
#' @describeIn GdObject-class Generics for `imageMap<-`.
#' @keywords internal
setGeneric("imageMap<-", function(GdObject, value) standardGeneric("imageMap<-"))

#' @describeIn  GdObject-class Replace the content of the `imageMap` slot.
#' @keywords internal
#' @export
setReplaceMethod(
    "imageMap", signature("GdObject", "ImageMapOrNULL"),
    function(GdObject, value) {
        GdObject@imageMap <- value
        return(GdObject)
    }
)

##  GdObject Methods drawAxis ------------------------------------------------

#' @describeIn GdObject-class Generics for `drawAxis`.
#' @export
setGeneric("drawAxis", function(GdObject, ...) standardGeneric("drawAxis"))

#' @describeIn GdObject-class add a y-axis to the title panel of a track if
#' necessary. Unless overwritten in one of the sub-classes this usually
#' does not plot anything and returns `NULL`.
#' @export
setMethod("drawAxis", signature(GdObject = "GdObject"), function(GdObject, ...) {
    return(NULL)
})

##  GdObject Methods drawGrid ------------------------------------------------

#' @describeIn GdObject-class Generics for `drawGrid`.
setGeneric("drawGrid", function(GdObject, ...) standardGeneric("drawGrid"))

# #' @describeIn GdObject-class superpose a grid on top of a track if necessary.
# #' Unless overwritten in one of the sub-classes this usually does not plot
# #' anything and returns `NULL`.
#' @keywords internal
setMethod("drawGrid", signature(GdObject = "GdObject"), function(GdObject, ...) {
    return(NULL)
})

##  GdObject Methods drawGd --------------------------------------------------

#' @describeIn GdObject-class Generics for `drawGD`.
#' @export
setGeneric("drawGD", function(GdObject, ...) standardGeneric("drawGD"))

##  GdObject Methods Annotation accessors ------------------------------------

.getAnn <- function(GdObject, type) {
    return(as.character(values(GdObject)[[type]]))
}
.setAnn <- function(GdObject, value, type) {
    v <- values(GdObject)
    if (length(value) > 1 && length(value) != nrow(v)) {
        stop(
            "The length of the replacement value for the '", type, "' annotation does not match the number ",
            "of features in the track."
        )
    }
    v[[type]] <- value
    mcols(GdObject@range) <- v
    return(GdObject)
}

#' @exportMethod gene
#' @describeIn GdObject-class Generics for `gene`.
#' @keywords internal
setGeneric("gene", function(GdObject, ...) standardGeneric("gene"))

#' @exportMethod "gene<-"
#' @describeIn GdObject-class Generics for `gene<-`.
#' @keywords internal
setGeneric("gene<-", function(GdObject, value) standardGeneric("gene<-"))

#' @exportMethod symbol
#' @describeIn GdObject-class Generics for `symbol`.
#' @keywords internal
setGeneric("symbol", function(GdObject, ...) standardGeneric("symbol"))

#' @exportMethod "symbol<-"
#' @describeIn GdObject-class Generics for `symbol<-`.
#' @keywords internal
setGeneric("symbol<-", function(GdObject, value) standardGeneric("symbol<-"))

#' @exportMethod transcript
#' @describeIn GdObject-class Generics for `transcript`.
#' @keywords internal
setGeneric("transcript", function(GdObject, ...) standardGeneric("transcript"))

#' @exportMethod "transcript<-"
#' @describeIn GdObject-class Generics for `transcript<-`.
#' @keywords internal
setGeneric("transcript<-", function(GdObject, value) standardGeneric("transcript<-"))

#' @exportMethod exon
#' @describeIn GdObject-class Generics for `exon`.
#' @keywords internal
setGeneric("exon", function(GdObject, ...) standardGeneric("exon"))

#' @exportMethod "exon<-"
#' @describeIn GdObject-class Generics for `exon<-`.
#' @keywords internal
setGeneric("exon<-", function(GdObject, value) standardGeneric("exon<-"))

#' @exportMethod feature
#' @describeIn GdObject-class Generics for `feature`.
#' @keywords internal
setGeneric("feature", function(GdObject, ...) standardGeneric("feature"))

#' @exportMethod "feature<-"
#' @describeIn GdObject-class Generics for `feature<-`.
#' @keywords internal
setGeneric("feature<-", function(GdObject, value) standardGeneric("feature<-"))

#' @exportMethod identifier
#' @describeIn GdObject-class Generics for `identifier`.
#' @keywords internal
setGeneric("identifier", function(GdObject, ...) standardGeneric("identifier"))

#' @exportMethod "identifier<-"
#' @describeIn GdObject-class Generics for `identifier<-`.
#' @keywords internal
setGeneric("identifier<-", function(GdObject, value) standardGeneric("identifier<-"))


##  GdObject Methods  General accessors --------------------------------------

#' @exportMethod chromosome
#' @describeIn GdObject-class Generics for `chromosome`.
#' @keywords internal
setGeneric("chromosome", function(GdObject, ...) standardGeneric("chromosome"))

#' @describeIn GdObject-class return the chromosome for which the track is defined.
#' @export
setMethod("chromosome", "GdObject", function(GdObject) {
    return(NULL)
})

#' @exportMethod "chromosome<-"
#' @describeIn GdObject-class Generics for `chromosome`.
#' @keywords internal
setGeneric("chromosome<-", function(GdObject, value) standardGeneric("chromosome<-"))

#' @describeIn GdObject-class replace the value of the track's chromosome. This
#' has to be a valid UCSC chromosome identifier or an integer or character
#' scalar that can be reasonably coerced into one.
#' @export
setReplaceMethod("chromosome", "GdObject", function(GdObject, value) {
    return(GdObject)
})

#' @exportMethod position
#' @describeIn GdObject-class Generics for `position`.
#' @keywords internal
setGeneric("position", function(GdObject, ...) standardGeneric("position"))

#' @describeIn GdObject-class return the track's genome.
#' @export
setMethod("genome", "GdObject", function(x) NULL)

#' @describeIn GdObject-class set the track's genome. Usually this has to be a
#' valid UCSC identifier, however this is not formally enforced here.
#' @export
setReplaceMethod("genome", "GdObject", function(x, value) {
    return(x)
})

##  GdObject Methods  Prepare tracks for plotting ----------------------------

#' @exportMethod consolidateTrack
#' @describeIn GdObject-class Generics for `consolidateTrack`.
#' @keywords internal
setGeneric("consolidateTrack", function(GdObject, ...) standardGeneric("consolidateTrack"))

#' @describeIn GdObject-class Consolidate.
#' Determine whether there is `alpha` settings or not, and add this information
#' as the internal display parameter `.__hasAlphaSupport`.
#' @export
# #' @keywords internal
setMethod("consolidateTrack", signature(GdObject = "GdObject"), function(GdObject, alpha, ...) {
    pars <- list(...)
    pars <- pars[names(pars) != ""]
    pars[[".__hasAlphaSupport"]] <- alpha
    displayPars(GdObject) <- pars
    return(GdObject)
})

#' @noRd
#' @keywords internal
setGeneric("collapseTrack", function(GdObject, ...) standardGeneric("collapseTrack"))

#' @exportMethod stacking
#' @describeIn GdObject-class Generics for `stacking`.
#' @keywords internal
setGeneric("stacking", function(GdObject, ...) standardGeneric("stacking"))

#' @exportMethod "stacking<-"
#' @describeIn GdObject-class Generics for `stacking<-`.
#' @keywords internal
setGeneric("stacking<-", function(GdObject, value) standardGeneric("stacking<-"))

#' @exportMethod stacks
#' @describeIn GdObject-class Generics for `stacks`.
#' @keywords internal
setGeneric("stacks", function(GdObject, ...) standardGeneric("stacks"))

#' @exportMethod setStacks
#' @describeIn GdObject-class Generics for ``.
#' @keywords internal
setGeneric("setStacks", function(GdObject, ...) standardGeneric("setStacks"))

#' @describeIn GdObject-class set stacks.
#' @export
setMethod("setStacks", "GdObject", function(GdObject, ...) GdObject)

#' @exportMethod setCoverage
#' @describeIn GdObject-class Generics for ``.
#' @keywords internal
setGeneric("setCoverage", function(GdObject, ...) standardGeneric("setCoverage"))

##  GdObject Methods Internal methods ----------------------------------------

#' @noRd
#' @keywords internal
setGeneric(".buildRange", function(range, start, end, width, ...) standardGeneric(".buildRange"))
