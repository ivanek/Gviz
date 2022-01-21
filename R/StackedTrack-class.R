#' @include RangeTrack-class.R
NULL

## StackedTrack --------------------------------------------------------------

#' StackedTrack class and methods
#'
#'
#' The virtual parent class for all track types in the Gviz package which
#' contain potentially overlapping annotation items that have to be stacked
#' when plotted.
#'
#' @template GdObject-class_slot
#' @template RangeTrack-class_slot
#' @template StackedTrack_slot
#' @template StackedTrack_param
#'
#' @name StackedTrack-class
#'
#' @return A virtual Class: No objects may be created from it.

#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @examples
#' ## This is a reference class therefore we show below
#' ## an example from AnnotationTrack
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
#'
#' ## Stacking
#' stacking(annTrack)
#' stacking(annTrack) <- "dense"
#' plotTracks(annTrack)
#' @exportClass StackedTrack
setClass("StackedTrack",
    representation = representation("VIRTUAL",
        stacking = "character",
        stacks = "numeric"
    ),
    prototype = prototype(
        name = "StackedTrack",
        stacking = "squish",
        stackingValues = c("hide", "dense", "squish", "pack", "full"),
        dp = DisplayPars(
            stackHeight = 0.75,
            reverseStacking = FALSE
        )
    ),
    contains = "RangeTrack"
)

## Initialize ----------------------------------------------------------------

## Need to fill the stacks slot here, don't want to recompute all the time
#' @describeIn StackedTrack-class Initialize.
#' @export
setMethod("initialize", "StackedTrack", function(.Object, stacking, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "StackedTrack")
    pt <- getClass("StackedTrack")@prototype
    if (!missing(stacking)) {
        if (!all(stacking %in% pt@stackingValues)) {
            stop(
                sprintf("Problem initializing %s, accepts the following values for 'stacking': ", class(.Object)),
                paste(pt@stackingValues, collapse = ", "), "\n"
            )
        }
        .Object@stacking <- stacking
        r <- list(...)$range
        ## stacks <- if(length(r)>0) disjointBins(ranges(r)) else 0
        ## .Object@stacks <- stacks
        .Object@stacks <- numeric()
    }
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## General accessors ---------------------------------------------------------
## Annotation Accessors ------------------------------------------------------
## Stacking ------------------------------------------------------------------

## Stacking controls what to do with overlapping annotation regions.

#' @describeIn StackedTrack-class  return the current stacking type.
#' @export
setMethod("stacking", "StackedTrack", function(GdObject) GdObject@stacking)


#' @describeIn StackedTrack-class  set the object's stacking type to one in
#' `c(hide, dense, squish, pack,full)`.
#' @export
setReplaceMethod(
    "stacking", c("StackedTrack", "character"),
    function(GdObject, value) {
        pt <- getClass("StackedTrack")@prototype
        if (!all(value %in% pt@stackingValues)) {
            stop(
                "Problem initializing StackedTrack,  need the following values for 'stacking':",
                paste(pt@stackingValues, collapse = ", "), "\n"
            )
        }
        GdObject@stacking <- value
        displayPars(GdObject) <- list(stacking = value)
        return(GdObject)
    }
)

## Recompute or return the stacking information for different types of StackedTrack objects. Stacking in needed when
## annotation regions in the objects overlap and when the stacking type is set to squish, full or pack. Since stacking
## can be dependent on the available space (if feature annotation is added) we need to be able to recompute this before
## we start the actual plotting. For the different sub-classes of StackedTracks we need different behaviour of setStacks:
##   o StackedTrack: there are no groups, so each feature can be treated separately
##   o AnnotationTrack and GeneRegionTrack: features can be grouped, in which case we have to avoid overlapping of the whole group region,
##      i.e, from the start of the first group item to the end of the last. In addition to grouping we have to factor in additional space
##      before each group if group/transcript annotation is enabled (gpar showId==TRUE). To do so we need to figure out the current
##      fontsize on the device, which means that a device already has to be open and the appropriate viewport has been
##      pushed to the stack. Hence we have to call setStacks immediately before the actual plotting.
## 'stacks' should return a vector of stacks, where each factor level of the vector indicates membeship
## to a particular stacking level. 'setStacks' returns the updated GdObject.

#' @describeIn StackedTrack-class return the stack indices for each track item.
#' @export
setMethod("stacks", "StackedTrack",
    function(GdObject) if (length(GdObject@stacks)) GdObject@stacks else NULL
)


#' @describeIn StackedTrack-class  recompute the stacks based on the available
#' space and on the object's track items and stacking settings.
#' @export
setMethod("setStacks", "StackedTrack", function(GdObject, ...) {
    bins <- if (!.needsStacking(GdObject)) rep(1, length(GdObject)) else disjointBins(range(GdObject))
    GdObject@stacks <- bins
    return(GdObject)
})

## Consolidate ---------------------------------------------------------------

#' @describeIn StackedTrack-class Consolidate.
# 'For `StackedTrack`s set the stacking (which could have been passed in as
#' a display parameter)
# #' @keywords internal
#' @export
setMethod("consolidateTrack", signature(GdObject = "StackedTrack"), function(GdObject, ...) {
    GdObject <- callNextMethod()
    st <- .dpOrDefault(GdObject, "stacking")
    if (!is.null(st)) {
        stacking(GdObject) <- st
    }
    return(GdObject)
})

## Collapse  -----------------------------------------------------------------
## There is a natural limit of what can be plotted as individual features caused by the maximum resolution of the device.
## Essentially no object can be smaller than the equivalent of a single pixel on the screen or whatever a pixel corresponds
## to on other devices.
## Thus we need to adapt the all objects in the track to the current resolution. Things that are closer than a certain limit
## distance can be collapsed (if that is desired), things that are smaller than 1 pixel can be blown up to a minimum size.
## All collapseTrack methods should always return a GdObject instance of similar type as the input to allow us to keep the
## collpasing step optional and for all the downstream plotting operations to still work. The internal (mostly the GRanges)
## objects can be modified to whatever is desired. Please note that this method is called after the drawing viewport has been
## pushed, and hence all the coordinate systems are already in place.
## Available arguments are:
##    o GdObject: the input AnnotationTrack or GeneRegionTrack track object
##    o min.width: the minimum width in pixel, everything that's smaller will be expanded to this size.
##    o min.distance: the minimum distance between two features in pixels below which to start collapsing track items.
##    o collapse: logical, collapse overlapping items into a single meta-item.
##    o diff: the equivalent of 1 pixel in the native coordinate system.
##    o xrange: the data range on the x axis. Can be used for some preliminary subsetting to speed things up
## Subset --------------------------------------------------------------------

#' @describeIn StackedTrack-class subset the items in the `StackedTrack` object.
#' This is essentially similar to subsetting of the `GRanges` object in the
#' range slot. For most applications, the subset method may be more appropriate.
#' @export
setMethod("[", signature(x = "StackedTrack"), function(x, i, j, ..., drop = TRUE) {
    x <- callNextMethod(x, i)
    x@stacks <- x@stacks[i]
    return(x)
})

#' @describeIn StackedTrack-class subset a `StackedTrack` by coordinates and sort if necessary.
#' @export
setMethod("subset", signature(x = "StackedTrack"), function(x, from = NULL, to = NULL, sort = FALSE, stacks = FALSE, ...) {
    x <- callNextMethod(x = x, from = from, to = to, sort = sort)
    if (stacks) {
        x <- setStacks(x)
    }
    return(x)
})


## Position ------------------------------------------------------------------
## In some cases we don't need a range but rather a single position. Most of the time this is simply taking the geometric
## mean of the range. If numeric values are associated to positions we have to be able to extract those as well.

## DrawGrid ------------------------------------------------------------------
## Draw a grid in the background of a GdObject. For some subclasses this is meaningless, and the default function will
## return NULL without plotting anything.

## DrawAxis ------------------------------------------------------------------
## The individual bits and pieces of a Gviz plot are all drawn by separate renderers. Currently, those are a y-axis,
## a grid, and the actual track panel.
##
## For certain GdObject subclasses we may want to draw a y-axis. For others an axis is meaningless, and the default function
## will return NULL without plotting anything.

## DrawGD --------------------------------------------------------------------

## All the drawGD methods should support two modes, triggered by the boolean argument 'prepare':
##    In prepare mode: nothing is plotted but the object is prepared for plotting bases on the available space. The return
##       value of the method in this mode should always be the updated object. If nothing needs to be prepared, i.e., if the
##       plotting is independent from the available space, simply return the original object
##    In plotting mode: the object is plotted. Return value is the object with optional HTML image map information
##    added to the imageMap slot
## Since subsetting can be potentially expensive when the data are large we want to minimize this operation. Essentially it
## should be done only once before any other plotting or computation starts, hence we expect the GdObject in the drawGD
## methods to already be trimmed to the correct size

## The default method for all `StackedTrack` types which always should be called
## (this has to be done explicitly using `callNextMethod`)

## Although the stacking type is not stored as a displayParameter we still want to check whether it is
## included there and set the actual stacking of the object accordingly

#' @describeIn StackedTrack-class plot the object to a graphics device.
#' The return value of this method is the input object, potentially updated
#' during the plotting operation. Internally, there are two modes in which
#' the method can be called. Either in 'prepare' mode, in which case no
#' plotting is done but the stacking information is updated based on the
#' available space, or in 'plotting' mode, in which case the actual graphical
#' output is created. Note that the method for this particular subclass is
#' usually called through inheritance and not particularly useful on its own.
#' @export
setMethod("drawGD", signature("StackedTrack"), function(GdObject, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if ((is.logical(debug) && debug)) {
        browser()
    }
    st <- .dpOrDefault(GdObject, "stacking")
    if (!is.null(st)) {
        stacking(GdObject) <- st
    }
    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------
## Show ----------------------------------------------------------------------
