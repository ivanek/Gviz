#' @include DataTrack-class.R
#' @include ReferenceTrack-class.R

NULL

## AnnotationTrack Class -----------------------------------------------------

#' AnnotationTrack class and methods
#'
#' A fairly generic track object for arbitrary genomic range annotations, with
#' the option of grouped track items. The extended `DetailsAnnotationTrack`
#' provides a more flexible interface to add user-defined custom
#' information for each range.
#'
#' @name AnnotationTrack-class
#'
#' @template GdObject-class_slot
#' @template RangeTrack-class_slot
#' @template StackedTrack_slot
#' @template AnnotationTrack_slot
#'
#' @template AnnotationTrack_param
#'
#' @author Florian Hahne, Arne Mueller
#'
#' @inherit GdObject-class seealso
#'
#' @examples
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
#'
#' ## Or from a data.frame
#' df <- data.frame(
#'     start = st, end = ed, strand = str, id = paste("annTrack item", 1:4),
#'     feature = "test", group = gr
#' )
#' annTrack <- AnnotationTrack(
#'     range = df, genome = "hg19", chromosome = 7,
#'     name = "generic annotation", stacking = "squish"
#' )
#'
#' ## Or from a GRanges object
#' gr <- GRanges(
#'     seqnames = "chr7", range = IRanges(start = df$start, end = df$end),
#'     strand = str
#' )
#' genome(gr) <- "hg19"
#' mcols(gr) <- df[, -(1:3)]
#' annTrack <- AnnotationTrack(
#'     range = gr, name = "generic annotation",
#'     stacking = "squish"
#' )
#'
#' ## Finally from a GRangesList
#' grl <- split(gr, values(gr)$group)
#' AnnotationTrack(grl)
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
#' ## Track names
#' names(annTrack)
#' names(annTrack) <- "foo"
#' plotTracks(annTrack)
#'
#' ## Subsetting and splitting
#' subTrack <- subset(annTrack, to = 2155000)
#' length(subTrack)
#' subTrack[1:2]
#' split(annTrack, c(1, 2, 1, 2))
#'
#' ## Accessors
#' start(annTrack)
#' end(annTrack)
#' width(annTrack)
#' position(annTrack)
#' width(subTrack) <- width(subTrack) + 1000
#'
#' strand(annTrack)
#' strand(subTrack) <- "-"
#'
#' chromosome(annTrack)
#' chromosome(subTrack) <- "chrX"
#'
#' genome(annTrack)
#' genome(subTrack) <- "mm9"
#'
#' range(annTrack)
#' ranges(annTrack)
#'
#' ## Annotation
#' identifier(annTrack)
#' identifier(annTrack, "lowest")
#' identifier(subTrack) <- "bar"
#'
#' feature(annTrack)
#' feature(subTrack) <- "foo"
#'
#' values(annTrack)
#'
#' ## Grouping
#' group(annTrack)
#' group(subTrack) <- "Group 1"
#' chromosome(subTrack) <- "chr7"
#' plotTracks(subTrack)
#'
#' ## Stacking
#' stacking(annTrack)
#' stacking(annTrack) <- "dense"
#' plotTracks(annTrack)
#'
#' ## coercion
#' as(annTrack, "data.frame")
#' as(annTrack, "UCSCData")
#'
#' ## HTML image map
#' coords(annTrack)
#' tags(annTrack)
#' annTrack <- plotTracks(annTrack)$foo
#' coords(annTrack)
#' tags(annTrack)
#'
#' ## DetailsAnnotationTrack
#' library(lattice) # need to use grid grapics
#'
#' ## generate two random distributions per row (probe/feature)
#' ## the difference between the distributions increases from probe 1 to 4
#' m <- matrix(c(rgamma(400, 1)), ncol = 100)
#' m[, 51:100] <- m[, 51:100] + 0:3
#' ## rownames must be accessible by AnnotationTrack element identifier
#' rownames(m) <- identifier(annTrack, "lowest")
#'
#' ## create a lattice density plot for the values (signals) of the two groups
#' ## as the chart must be placed into a pre-set grid view port we have to use
#' ## print without calling plot.new! Note, use a common prefix for all lattice.
#' ## Avoid wasting space by removing y-axis decorations.
#'
#' ## Note, in this example 'm' will be found in the environment the 'details'
#' ## function is defined in. To avoid overwriting 'm' you should use a closure
#' ## or environment to access 'm'.
#' details <- function(identifier, ...) {
#'     d <- data.frame(signal = m[identifier, ], group = rep(c("grp1", "grp2"), each = 50))
#'     print(densityplot(~signal,
#'         group = group, data = d, main = identifier,
#'         scales = list(draw = FALSE, x = list(draw = TRUE)), ylab = "", xlab = "",
#'     ), newpage = FALSE, prefix = "plot")
#' }
#'
#' deTrack <- AnnotationTrack(
#'     range = gr, genome = "hg19", chromosome = 7,
#'     name = "generic annotation with details per entry", stacking = "squish",
#'     fun = details, details.ratio = 1
#' )
#'
#' plotTracks(deTrack)
#'
#' set.seed(1234)
#' deTrack <- AnnotationTrack(
#'     range = gr, genome = "hg19", chromosome = 7,
#'     name = "generic annotation with details per entry",
#'     stacking = "squish", fun = details,
#'     details.ratio = 1, selectFun = function(...) {
#'         sample(c(FALSE, TRUE), 1)
#'     }
#' )
#'
#' plotTracks(deTrack)
#' @importClassesFrom rtracklayer UCSCData
#'
#' @exportClass AnnotationTrack
setClass("AnnotationTrack",
    contains = "StackedTrack",
    prototype = prototype(
        columns = c("feature", "group", "id"),
        stacking = "squish",
        name = "AnnotationTrack",
        dp = DisplayPars(
            arrowHeadWidth = 30,
            arrowHeadMaxWidth = 40,
            cex.group = 0.6,
            cex = 1,
            col.line = "darkgray",
            col = "transparent",
            featureAnnotation = NULL,
            fill = "lightblue",
            fontfamily.group = "sans",
            fontcolor.group = .DEFAULT_SHADED_COL,
            fontcolor.item = "white",
            fontface.group = 2,
            fontsize.group = 12,
            groupAnnotation = NULL,
            just.group = "left",
            lex = 1,
            lineheight = 1,
            lty = "solid",
            lwd = 1,
            mergeGroups = FALSE,
            min.height = 3,
            min.width = 1,
            rotation = 0,
            rotation.group = 0,
            rotation.item = 0,
            shape = "arrow",
            showFeatureId = FALSE,
            showId = FALSE,
            showOverplotting = FALSE,
            size = 1
        )
    )
)

## Initialize ----------------------------------------------------------------

## Essentially we just check for the correct GRanges columns here
#' @describeIn AnnotationTrack-class Show method.
#' @export
setMethod("initialize", "AnnotationTrack", function(.Object, ...) {
    if (is.null(list(...)$range) && is.null(list(...)$genome) && is.null(list(...)$chromosome)) {
        return(.Object)
    }
    ## the display parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "AnnotationTrack")
    range <- list(...)$range
    if (!is.null(range) && length(.Object)) {
        if (!all(.Object@columns %in% colnames(values(range)))) {
            stop(
                "Problem initializing AnnotationTrack, need the following columns:",
                toString(.Object@columns, collapse = ", ")
            )
        }
        grp <- if (is(.Object, "GeneRegionTrack")) values(range)$transcript else values(range)$group
        if (any(vapply(split(as.character(strand(range)), grp), function(x) length(unique(x)), numeric(1)) != 1)) {
            stop("Grouped elments of a RangeTrack can not be on opposing strands")
        }
    }
    .Object <- callNextMethod()
    return(.Object)
})

## The file-based version of the AnnotationTrack class. This will mainly provide a means to dispatch to
## a special 'subset' method which should stream the necessary data from disk.

#' @describeIn AnnotationTrack-class The file-based version of the `AnnotationTrack-class`.
#' @exportClass ReferenceAnnotationTrack
setClass("ReferenceAnnotationTrack", contains = c("AnnotationTrack", "ReferenceTrack"))

## Initialize ----------------------------------------------------------------

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritance has some strange features with regards to method selection
#' @describeIn AnnotationTrack-class Initialize.
#' @export
setMethod("initialize", "ReferenceAnnotationTrack", function(.Object, stream, reference, mapping = list(),
                                                             args = list(), defaults = list(), ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object = .Object, reference = reference, stream = stream,
        mapping = mapping, args = args, defaults = defaults)
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

## Constructor. The following arguments are supported:
##    o range: a data.frame or a GRanges object containing the information
##       about the track items. If a data.frame, it needs to be coerceable
##     to a GRanges object, i.e., it needs at least the mandatory 'start', 'stop' and
##     'strand' columns. Additional optional columns are:
##        - feature: the type of the item. Can be mapped to colors via the DisplayPars.
##        - group: a grouping factor to connect track items.
##        - id: a unique identifier for a feature. This will be plotted if showId==TRUE
##             Note that internally we use the value of ID as the seqnames slot in the
##             internal GRanges object. Defaults for all missing columns are generated.
##    Instead of using the 'range' parameter, all these values can also be passed as
##    individual vectors, in which case they need to be of similar length.
##    o start, end, width: numeric vectors of the item start and end coordinates
##    o strand: the strand information may be provided in the form '+' for
##       the Watson strand, '-' for the Crick strand or '*' for any of the two.
##    o feature, group, id: individual vectors of equal length as described above.
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o stacking: character controlling the stacking of overlapping items. One in 'hide',
##       'dense', 'squish', 'pack' or 'full'.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as further DisplayParameters

#' @describeIn AnnotationTrack-class Constructor function for
#' `AnnotationTrack-class`

#' @return
#' The return value of the constructor function is a new object of class
#' \code{AnnotationTrack} or of class \code{DetailsAnnotationTrack}, depending
#' on the constructor arguments. Typically the user will not have to be
#' troubled with this distinction and can rely on the constructor to make the
#' right choice.
#'
#' @section Objects from the class:
#'
#' Objects can be created using the constructor function
#' \code{AnnotationTrack}.
#' @export
AnnotationTrack <- function(range = NULL, start = NULL, end = NULL, width = NULL, feature, group, id, strand, chromosome,
                            genome, stacking = "squish", name = "AnnotationTrack", fun, selectFun, importFunction,
                            stream = FALSE, ...) {
    ## Some defaults
    covars <- .getCovars(range)
    isStream <- FALSE
    if (!is.character(range)) {
        n <- max(c(length(start), length(end), length(width)), nrow(covars))
        if (is.null(covars[["feature"]]) && missing(feature)) {
            feature <- rep("unknown", n)
        }
        if (is.null(covars[["id"]]) && missing(id)) {
            id <- make.unique(rep(if (!is.null(feature)) as.character(feature) else covars[["feature"]], n)[seq_len(n)])
        }
        if (is.null(covars[["group"]]) && missing(group)) {
            group <- seq_len(n)
        }
    }
    ## Build a GRanges object from the inputs
    .missingToNull(c("feature", "group", "id", "strand", "chromosome", "importFunction", "genome"))
    args <- list(feature = feature, group = group, id = id, strand = strand, chromosome = chromosome, genome = genome)
    defs <- list(feature = "unknown", group = "unknown", id = "unknown", strand = "*", density = 1, chromosome = "chrNA", genome = NA)
    range <- .buildRange(
        range = range, groupId = "group", start = start, end = end, width = width,
        args = args, defaults = defs, chromosome = chromosome, trackType = "AnnotationTrack",
        importFun = importFunction, stream = stream
    )
    if (is.list(range)) {
        isStream <- TRUE
        slist <- range
        range <- GRanges()
    }
    ## Pipes have a special meaning for merged groups, so we can't have them in the initial group vector
    mcols(range)[["group"]] <- gsub("|", "", mcols(range)[["group"]], fixed = TRUE)
    ## If no chromosome was explicitely asked for we just take the first one in the GRanges object
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- if (length(range) > 0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    }
    ## And finally the object instantiation, we have to distinguish between DetailsAnnotationTracks and normal ones
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    if (missing(fun)) {
        if (!isStream) {
            return(new("AnnotationTrack",
                chromosome = as.character(chromosome[1]), range = range,
                name = name, genome = genome, stacking = stacking, ...
            ))
        } else {
            ## A bit hackish but for some functions we may want to know which track type we need but at the
            ## same time we do not want to enforce this as an additional argument
            e <- new.env()
            e[["._trackType"]] <- "AnnotationTrack"
            environment(slist[["stream"]]) <- e
            return(new("ReferenceAnnotationTrack",
                chromosome = as.character(chromosome[1]), range = range,
                name = name, genome = genome, stacking = stacking, stream = slist[["stream"]], reference = slist[["reference"]],
                mapping = slist[["mapping"]], args = args, defaults = defs, ...
            ))
        }
    } else {
        if (!is.function(fun)) {
            stop("'fun' must be a function")
        }
        if (missing(selectFun)) {
            selectFun <- function(...) {
                return(TRUE)
            }
        }
        if (!is.function(selectFun)) {
            stop("'selectFun' must be a function")
        }
        return(new("DetailsAnnotationTrack",
            chromosome = as.character(chromosome[1]), range = range,
            name = name, genome = genome, stacking = stacking, fun = fun, selectFun = selectFun, ...
        ))
    }
}

#' @describeIn AnnotationTrack-class directly extends `AnnotationTrack.`
#' @exportClass DetailsAnnotationTrack
setClass("DetailsAnnotationTrack",
    contains = "AnnotationTrack",
    representation = representation(fun = "function", selectFun = "function"),
    prototype = prototype(
        fun = function(...) {},
        selectFun = function(...) {
            return(TRUE)
        },
        dp = DisplayPars(
            details.minWidth = 100,
            details.ratio = Inf,
            details.size = 0.5,
            detailsBorder.col = "darkgray",
            detailsBorder.fill = "transparent",
            detailsBorder.lty = "solid",
            detailsBorder.lwd = 1,
            detailsConnector.cex = 1,
            detailsConnector.col = "darkgray",
            detailsConnector.lty = "dashed",
            detailsConnector.lwd = 1,
            detailsConnector.pch = 20,
            detailsFunArgs = list(),
            groupDetails = FALSE
        )
    )
)

#' @describeIn AnnotationTrack-class Constructor function for
#' `DetailsAnnotationTrack-class`
#'
#' The `DetailsAnnotationTrack` class directly extends `AnnotationTrack.`
#' The purpose of this track type is to add an arbitrarily detailed plot
#' section (typically consisting of additional quantitative data) for each
#' range element of an `AnnotationTrack.` This allows a locus wide view of
#' annotation elements together with any kind of details per feature or element
#'  that may for instance provide insight on how some complex quantitative
#'  measurements change according to their position in a locus. If the
#'  quantitative data is too complex for a `DataTrack` e.g. because it requires
#'  extra space or a trellis-like representation, a `DetailsAnnotationTrack` can
#'  be used instead. Example: An `AnnotationTrack` shows the positions of a
#'  number of probes from a microarray, and you want a histogram of the signal
#'  intensity distribution derived from all samples at each of these probe
#'  location. Another example usage would be to show for each element of an
#'  `AnnotationTrack` an xy-plot of the signal against some clinical measurement
#'  such as blood pressure. The limitation for applications of this type of
#'  track is basically only the available space of the device you are
#'  plotting to.
#'
#' This flexibility is possible by utilizing a simple function model
#' to perform all the detailed plotting. The functionality of this plotting
#' function fun is totally up to the user, and the function environment is
#' prepared in a way that all necessary information about the plotted
#' annotation feature is available. To restrict the details section to only
#' selected number of annotation features one can supply another function
#' `selectFun`, which decides for each feature separately whether details are
#' available or not. Finally, an arbitrary number of additional arguments can
#' be passed on to these two function by means of the `detailsFunArgs` display
#' parameter. This is expected to be a named list, and all list elements are
#' passed along to the plotting function fun and to the selector function
#' `selectFun` as additional named arguments. Please note that some argument
#' names like `start`, `end` or `identifier` are reserved and can not be used
#' in the `detailsFunArgs` list. For examples of plotting functions,
#' see the 'Examples' section.
#' @export
DetailsAnnotationTrack <- function(...) AnnotationTrack(...)

## Initialize ----------------------------------------------------------------

#' @describeIn AnnotationTrack-class Initialize.
#' @export
setMethod("initialize", "DetailsAnnotationTrack", function(.Object, fun, selectFun, ...) {
    ## the diplay parameter defaults
    .Object <- .updatePars(.Object, "DetailsAnnotationTrack")
    .makeParMapping()
    .Object@fun <- fun
    .Object@selectFun <- selectFun
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## General accessors ---------------------------------------------------------
## Annotation Accessors ------------------------------------------------------

#' @describeIn AnnotationTrack-class extract the group membership for all track items.
#' @export
setMethod("group", "AnnotationTrack", function(GdObject) .getAnn(GdObject, "group"))

#' @describeIn AnnotationTrack-class replace the grouping information for track items.
#' The replacement value must be a factor of appropriate length or another
#' vector that can be coerced into such.
#' @export
setReplaceMethod("group", signature("AnnotationTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "group"))

## Context-dependent meta-accessors to the identifier data of an AnnotationTrack and a GeneRegionTrack. For the former, those will
## be the content of the 'id', the 'group' or the 'feature' metadata column, for the latter, one in
## the selection of 'symbol', 'gene', 'transcript', 'feature' or 'exon'. For historical reasons, logical values are also allowed, where TRUE
## is equivalent to 'symbol' and FALSE is equivalent to 'gene'. If not provided explicitely in the 'type' argument, the identifier
## type is inferred from the 'transcriptAnnotation' display parameter (a.k.a. 'geneSymbols') for a GeneRegionTrack, and from the
## 'groupAnnotation' parameter for an AnnotationTrack. The special value 'lowest' should work across all classes and provided the most
## detailed annotation ('id' for AnnotationTrack and 'exon' for GeneRegionTrack).

#' @describeIn AnnotationTrack-class return track item identifiers.
#' Depending on the setting of the optional argument lowest, these are either
#' the group identifiers or the individual item identifiers.
#' @export
setMethod("identifier", "AnnotationTrack", function(GdObject, type = .dpOrDefault(GdObject, "groupAnnotation", "group")) {
    if (is.null(type)) {
        type <- "group"
    }
    id <- switch(as.character(type),
        "group" = group(GdObject),
        "id" = .getAnn(GdObject, "id"),
        "feature" = feature(GdObject),
        "lowest" = .getAnn(GdObject, "id"),
        group(GdObject)
    )
    id[is.na(id)] <- "NA"
    return(id)
})

#' @describeIn AnnotationTrack-class Set the track item identifiers.
#' The replacement value has to be a character vector of appropriate length.
#' This always replaces the group-level identifiers, so essentially it is
#' similar to `groups<-`.
#' @export
setReplaceMethod("identifier", c("AnnotationTrack", "character"), function(GdObject, value) {
    type <- .dpOrDefault(GdObject, "groupAnnotation", "group")
    switch(as.character(type),
        "group" = group(GdObject) <- value,
        "id" = ranges(GdObject)$id <- value,
        "feature" = feature(GdObject) <- value,
        "lowest" = ranges(GdObject)$id <- value,
        group(GdObject) <- value
    )
    return(GdObject)
})
## Stacking ------------------------------------------------------------------

#' @describeIn AnnotationTrack-class Recompute the stacks based on the available
#' space and on the object's track items and stacking settings.
#' @export
setMethod("setStacks", "AnnotationTrack", function(GdObject, recomputeRanges = TRUE) {
    if (!.needsStacking(GdObject) || length(GdObject) == 0) {
        ## No stacks needed, so there is just a single bin
        bins <- rep(1, length(GdObject))
    } else {
        gp <- group(GdObject)
        if (!is.factor(gp)) {
            gp <- factor(gp, levels=unique(gp))
        }
        needsGrp <- any(duplicated(gp))
        lranges <- .dpOrDefault(GdObject, ".__groupRanges")
        gpt <- if (needsGrp) table(gp) else rep(1, length(GdObject))
        if (recomputeRanges || is.null(lranges) || length(lranges) != length(gpt)) {
            GdObject <- .computeGroupRange(GdObject)
            lranges <- .dpOrDefault(GdObject, ".__groupRanges")
        }
        uid <- if (.transcriptsAreCollapsed(GdObject)) {
            sprintf("uid%i", seq_along(identifier(GdObject)))
        } else {
            make.unique(identifier(GdObject, type = "lowest"))
        }

        bins <- rep(disjointBins(lranges), gpt)
        names(bins) <- if (needsGrp) unlist(split(uid, gp)) else uid
        bins <- bins[uid]
    }
    bins <- if (length(bins)) bins else 0
    GdObject@stacks <- bins
    return(GdObject)
})
## Consolidate ---------------------------------------------------------------

#' @describeIn AnnotationTrack-class Consolidate.
#' Determine whether there is group label annotation or not, and add this
#' information as the internal display parameter `.__hasAnno`. Precompute
#' the grouped ranges together with optional labels in order to determine
#' the correct plotting range later.
#' @keywords internal
setMethod("consolidateTrack", signature(GdObject = "AnnotationTrack"), function(GdObject, hasAxis = FALSE,
                                                                                hasTitle = .dpOrDefault(GdObject, "showTitle", TRUE),
                                                                                title.width = NULL, ...) {
    GdObject <- callNextMethod(GdObject, ...)
    if (length(GdObject)) {
        ## ids are shown if either set by the showId parameter, or through the use of the transcriptAnnotation or
        ## groupAnnotation parameters
        ids <- identifier(GdObject)
        sid <- .dpOrDefault(GdObject, "showId")
        ta <- if (is(GdObject, "GeneRegionTrack")) .dpOrDefault(GdObject, "transcriptAnnotation") else .dpOrDefault(GdObject, "groupAnnotation")
        sid <- if (is.null(sid)) !is.null(ta) && ta != "none" else sid
        hasAnno <- sid && !all(ids == "")
        ## element ids are shown if either the showFeatureId or showExonId parameters are TRUE, or if featureAnnotation
        ## or exonAnnotation is set
        sfid <- if (is(GdObject, "GeneRegionTrack")) .dpOrDefault(GdObject, "showExonId") else .dpOrDefault(GdObject, "showFeatureId")
        fa <- if (is(GdObject, "GeneRegionTrack")) .dpOrDefault(GdObject, "exonAnnotation") else .dpOrDefault(GdObject, "featureAnnotation")
        sfid <- if (is.null(sfid)) !is.null(fa) && fa != "none" else sfid
        displayPars(GdObject) <- list(".__hasAnno" = hasAnno, showId = sid, showFeatureId = sfid, showExonId = sfid)
        GdObject <- .computeGroupRange(GdObject, hasAxis = hasAxis, hasTitle = hasTitle, title.width = title.width)
    } else {
        displayPars(GdObject) <- list(".__hasAnno" = FALSE, showId = FALSE, showFeatureId = FALSE, showExonId = FALSE)
    }
    return(GdObject)
})

## Collapse  -----------------------------------------------------------------

## A slightly quicker function to compute overlaps between two GRanges objects
.myFindOverlaps <- function(gr1, gr2) {
    gr1 <- sort(gr1)
    gr2 <- sort(gr2)
    gr1 <- split(gr1, seqnames(gr1))
    gr2 <- split(gr2, seqnames(gr2))
    queryHits(findOverlaps(ranges(gr1), ranges(gr2)))
}
## Find all elements in a GRanges object 'grange' with distance smaller than 'minXDist' and merge them along with their additional
## metadata columns. 'elements' is a frequency table of items per group, and it is needed to figure out whether all items of a given
## group have been merged. 'GdObject' is the input track object from which certain information has to be extracted. The output of this
## function is a list with elements
##   o range: the updated GRanges object
##   o needsRestacking: logical flag indicating whether stacks have to be recomputed
##   o split: the original merged annotation (need to work on this, currently not used)
##   o merged: logical vector indicating which of the elements in 'range' constitute fully merged groups
.collapseAnnotation <- function(grange, minXDist, elements, GdObject, offset = 0) {
    needsRestacking <- TRUE
    annoSplit <- merged <- NULL
    anno <- as.data.frame(grange)
    for (i in colnames(anno)) {
        if (is.factor(anno[, i])) {
            anno[, i] <- as.character(anno[, i])
        }
    }
    cols <- c("strand", "density", "gdensity", "feature", "id", "start", "end", if (is(GdObject, "GeneRegionTrack")) {
        c("gene", "exon", "transcript", "symbol", "rank")
    } else {
        "group"
    })
    missing <- which(!cols %in% colnames(anno))
    for (i in missing) {
        anno[, cols[missing]] <- if (cols[i] == "density") 1 else NA
    }
    rRed <- if (length(grange) > 1) reduce(grange, min.gapwidth = minXDist, with.revmap = TRUE) else grange
    if (length(rRed) < length(grange)) {
        ## Some of the items have to be merged and we need to make sure that the additional annotation data that comes with it
        ## is processed in a sane way.
        needsRestacking <- TRUE
        mapping <- rep(seq_along(rRed$revmap), elementNROWS(rRed$revmap))
        ## We start by finding the items that have not been reduced
        identical <- mapping %in% which(table(mapping) == 1)
        newVals <- anno[identical, cols]
        ## Here we hijack the seqnames column to indicate whether the whole group has been merged
        if (nrow(newVals)) {
            newVals$seqnames <- elements[as.character(anno[identical, "seqnames"])] == 1
            newVals$gdensity <- ifelse(elements[as.character(anno[identical, "seqnames"])] == 1, 1, NA)
        }
        ## Now find out which original items have been merged
        grange <- grange[!identical]
        rRed <- rRed[-(mapping[identical])]
        index <- mapping[!identical]
        annoSplit <- split(anno[!identical, ], index)
        cid <- function(j) sprintf("[Cluster_%i]  ", j + offset)
        ## FIXME: We could speed this up by running it in C
        tmpVals <- lapply(seq_along(annoSplit), function(i) {
            x <- annoSplit[[i]]
            if (is(GdObject, "GeneRegionTrack")) {
                c(
                    strand = ifelse(length(unique(x[, "strand"])) == 1, as.character(x[1, "strand"]), "*"),
                    density = sum(as.integer(x[, "density"])),
                    gdensity = ifelse(is.na(head(x[, "gdensity"], 1)), 1, sum(as.integer(x[, "gdensity"]))),
                    feature = ifelse(length(unique(x[, "feature"])) == 1, as.character(x[1, "feature"]), "composite"),
                    id = ifelse(length(unique(x[, "id"])) == 1, as.character(x[1, "id"]), cid(i)),
                    start = min(x[, "start"]),
                    end = max(x[, "end"]),
                    gene = ifelse(length(unique(x[, "gene"])) == 1, as.character(x[1, "gene"]), cid(i)),
                    exon = ifelse(length(unique(x[, "exon"])) == 1, as.character(x[1, "exon"]), cid(i)),
                    transcript = ifelse(length(unique(x[, "transcript"])) == 1, as.character(x[1, "transcript"]), cid(i)),
                    symbol = ifelse(length(unique(x[, "symbol"])) == 1, as.character(x[1, "symbol"]), cid(i)),
                    rank = min(as.integer(x[, "rank"])), seqnames = as.vector(nrow(x) == elements[x[1, "seqnames"]])
                )
            } else {
                c(
                    strand = ifelse(length(unique(x[, "strand"])) == 1, as.character(x[1, "strand"]), "*"),
                    density = sum(as.integer(x[, "density"])),
                    gdensity = ifelse(is.na(head(x[, "gdensity"], 1)), 1, sum(as.integer(x[, "gdensity"]))),
                    feature = ifelse(length(unique(x[, "feature"])) == 1, as.character(x[1, "feature"]), "composite"),
                    id = ifelse(length(unique(x[, "id"])) == 1, as.character(x[1, "id"]), cid(i)),
                    start = min(x[, "start"]),
                    end = max(x[, "end"]),
                    group = ifelse(length(unique(x[, "group"])) == 1, as.character(x[1, "group"]), cid(i)),
                    seqnames = as.vector(nrow(x) == elements[x[1, "seqnames"]])
                )
            }
        })
        newVals <- rbind(newVals, as.data.frame(do.call(rbind, tmpVals), stringsAsFactors = FALSE))
        merged <- as.logical(newVals$seqnames)
        grange <- GRanges(
            seqnames = chromosome(GdObject), strand = newVals[, "strand"],
            ranges = IRanges(start = as.integer(newVals[, "start"]), end = as.integer(newVals[, "end"]))
        )
        cnMatch <- match(c(colnames(values(GdObject)), "gdensity"), colnames(newVals))
        mcols(grange) <-
            if (any(is.na(cnMatch))) newVals[, setdiff(colnames(newVals), c("strand", "start", "end", "seqnames"))] else newVals[, cnMatch]
    } else {
        grange2 <- GRanges(seqnames = chromosome(GdObject), strand = strand(grange), ranges = ranges(grange))
        mcols(grange2) <- mcols(grange)
        grange <- grange2
    }
    return(list(range = grange, needsRestacking = needsRestacking, split = annoSplit, merged = merged, offset = length(annoSplit)))
}

## For AnnotationTracks we need to collapse the all regions along with the additional annotation.
## For GeneRegionTracks we essentially need to do the same thing as for AnnotationTracks, however the additional annotation columns
## are quite different. We do this in multiple turn with increasing levels of complexity:
##    1.) merge all individual items within a group that can no longer be separated
##    2.) merge overlapping groups with just a single remaining item (optional, if mergeGroups==TRUE)

#' @describeIn AnnotationTrack-class preprocess the track before plotting.
#' This will collapse overlapping track items based on the available resolution
#' and increase the width and height of all track objects to a minimum value
#' to avoid rendering issues. See collapsing for details.
#' @keywords internal
setMethod(
    "collapseTrack", signature(GdObject = "AnnotationTrack"),
    function(GdObject, diff = .pxResolution(coord = "x"), xrange) {
        ## We first add the original unmodified GdObject as a display parameter to be able to reference back if we ever need to
        displayPars(GdObject) <- list(".__OriginalGdObject" = .deepCopyPars(GdObject))
        collapse <- .dpOrDefault(GdObject, "collapse", TRUE)
        min.width <- .dpOrDefault(GdObject, "min.width", 2)
        min.distance <- .dpOrDefault(GdObject, "min.distance", 2)
        minXDist <- max(0, ceiling(min.distance * diff))
        r <- ranges(GdObject)
        ## Compute native coordinate equivalent to 1 pixel and resize
        rNew <- .resize(r, min.width, diff)
        needsRestacking <- any(r != rNew)
        r <- rNew
        ## Collapse all items within a group to a single meta-item (if collapseTranscripts is set)
        if (is(GdObject, "GeneRegionTrack")) {
            ctrans <- .dpOrDefault(GdObject, "collapseTranscripts", FALSE)
            if (is.logical(ctrans) && ctrans) {
                ctrans <- "gene"
            }
            switch(ctrans,
                "gene" = {
                    newVals <- unlist(endoapply(split(values(r), paste(gene(GdObject), strand(GdObject))), head, 1))
                    newVals$exon <- NA
                    newVals$feature <- "merged"
                    newVals$transcript <- newVals$gene
                    r <- unlist(range(split(r, gene(GdObject))))
                    mcols(r) <- newVals
                    GdObject@range <- r
                },
                "longest" = {
                    r <- unlist(endoapply(split(r, gene(GdObject)), function(x) {
                        xs <- split(x, x$transcript)
                        xs[[which.max(vapply(xs, function(y) abs(diff(as.numeric(as.data.frame(ranges(range(y)))[, c("start", "end")]))), FUN.VALUE = numeric(1L)))]]
                    }))
                    GdObject@range <- r
                },
                "shortest" = {
                    r <- unlist(endoapply(split(r, gene(GdObject)), function(x) {
                        xs <- split(x, x$transcript)
                        xs[[which.min(vapply(xs, function(y) abs(diff(as.numeric(as.data.frame(ranges(range(y)))[, c("start", "end")]))), FUN.VALUE = numeric(1L)))]]
                    }))
                    GdObject@range <- r
                },
                "meta" = {
                    newVals <- unlist(endoapply(split(values(r), paste(gene(GdObject), strand(GdObject))), head, 1))
                    newVals$feature <- "merged"
                    newVals$transcript <- newVals$gene
                    rtmp <- reduce(split(r, paste(gene(GdObject), strand(GdObject))))
                    newVals <- newVals[rep(seq_along(rtmp), elementNROWS(rtmp)), ]
                    newVals$exon <- paste("merged_exon_", unlist(lapply(elementNROWS(rtmp), function(x) seq(1, x)), use.names = FALSE), sep = "")
                    r <- unlist(rtmp)
                    mcols(r) <- newVals
                    GdObject@range <- r
                }
            )
        }
        ## Collapse overlapping ranges (less than minXDist space between them) and process the annotation data
        if (collapse) {
            ## Merge all items in those groups for which no individual items can be separated
            elements <- table(group(GdObject))
            rr <- GRanges(seqnames = as.character(group(GdObject)), ranges = IRanges(start = start(r), end = end(r)), strand = strand(r))
            mcols(rr) <- mcols(r)
            rr <- sort(unique(rr))
            mergedAnn <- .collapseAnnotation(rr, minXDist, elements, GdObject)
            needsRestacking <- needsRestacking || mergedAnn$needsRestacking
            ## Now we take a look whether there are any groups that could be merged (if mergeGroups is TRUE)
            if (.dpOrDefault(GdObject, "mergeGroups", FALSE) && any(mergedAnn$merged)) {
                rr <- sort(mergedAnn$range[mergedAnn$merged])
                strand(rr) <- "*"
                mergedAnn2 <- .collapseAnnotation(rr, minXDist, elements, GdObject, mergedAnn$offset)
                needsRestacking <- needsRestacking || mergedAnn2$needsRestacking
                mergedAnn$range <- c(mergedAnn$range[!mergedAnn$merged], mergedAnn2$range)
            }
            r <- mergedAnn$range
        }
        ## Reconstuct the track object and return
        GdObject@range <- r
        ## if(needsRestacking)
        GdObject <- setStacks(GdObject)
        return(GdObject)
    }
)

## Subset --------------------------------------------------------------------

## In order to keep the grouping information for track regions in the clipped areas we have to
## keep all group elements that overlap with the range. We still want to record the requested
## ranges in the internal '.__plottingRange' display parameter.

#' @describeIn AnnotationTrack-class subset a `AnnotationTrack` by coordinates
#' and sort if necessary.
#' @export
setMethod("subset", signature(x = "AnnotationTrack"), function(x, from = NULL, to = NULL, sort = FALSE, stacks = FALSE, use.defaults = TRUE, ...) {
    ## Subset to a single chromosome first
    lx <- length(x)
    csel <- seqnames(x) != chromosome(x)
    if (any(csel)) {
        x <- x[!csel]
    }
    if (length(x)) {
        ## Nothing to do if everything is within the range
        granges <- unlist(range(split(ranges(x), group(x))))
        ranges <- if (use.defaults) {
            .defaultRange(x, from = from, to = to)
        } else {
            c(
                from = ifelse(is.null(from), min(start(granges)) - 1, from),
                to = ifelse(is.null(to), max(end(granges)) + 1, to)
            )
        }
        if (!(any(end(x) < ranges["from"] | start(x) > ranges["to"]))) {
            if (stacks) {
                x <- setStacks(x)
            }
            return(x)
        }
        ## Now remove everything except for the overlapping groups by first subselecting all groups in the range...
        gsel <- names(granges)[subjectHits(findOverlaps(GRanges(seqnames = chromosome(x), ranges = IRanges(min(ranges), max(ranges))), granges))]
        x <- x[group(x) %in% gsel]
        if (sort) {
            x <- x[order(range(x)), ]
        }
        if (stacks) {
            x <- setStacks(x)
        }
        displayPars(x) <- list(".__plottingRange" = ranges)
    }
    if (length(x) != lx) {
        x <- .computeGroupRange(x)
    }
    return(x)
})

## ReferenceDataTracks need to stream the data from file and then pass the results on to the next method

#' @describeIn AnnotationTrack-class subset a `ReferenceAnnotationTrack` by
#' coordinates and sort if necessary.
#' @export
setMethod("subset", signature(x = "ReferenceAnnotationTrack"), function(x, from, to, chromosome, ...) {
    ## We only need to reach out into the referenced file once if the range is already contained in the object
    if (missing(from) || is.null(from) || missing(to) || is.null(to)) {
        stop("Need both start and end location to subset a ReferenceAnnotationTrack")
    }
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- Gviz::chromosome(x)
    }
    subRegion <- GRanges(seqnames = chromosome[1], ranges = IRanges(start = from, end = to))
    if (length(ranges(x)) == 0 || all(overlapsAny(ranges(x), subRegion))) {
        cMap <- .resolveColMapping(x@stream(x@reference, subRegion), x@args, x@mapping)
        x@range <- .buildRange(cMap$data, args = cMap$args, defaults = x@defaults, trackType = "AnnotationTrack")
        chromosome(x) <- chromosome[1]
    }
    return(callNextMethod(x = x, from = from, to = to, drop = FALSE, ...))
})

## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------
setMethod("drawGrid", signature(GdObject = "AnnotationTrack"), function(GdObject, from, to) {
    if (.dpOrDefault(GdObject, "grid", FALSE)) {
        pushViewport(dataViewport(xData = c(from, to), extension = c(0, 0), yData = 0:1, clip = TRUE))
        panel.grid(
            h = 0, v = .dpOrDefault(GdObject, "v", -1),
            col = .dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty = .dpOrDefault(GdObject, "lty.grid", 1),
            lwd = .dpOrDefault(GdObject, "lwd.grid", 1)
        )
        popViewport(1)
    }
})
## DrawAxis ------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------

## Draw gene models as found in AnnotationTracks or GeneRegionTracks
##
## Calculate all coordinates and values for the individual stacks first, and append to the
## input list 'toDraw'. This allows us to plot the whole track at once, making use of grid's
## vectorization. Without this tweak, every stack would have to be drawn individually, which
## can be painfully slow...
## Because all of the y coordinates are calculated in the coordinate system of the current
## viewport for the respective stack we have to add offset values.

## Compute the coordinates and colors for all track items (e.g. exons in a GeneRegionTrack)
.boxes <- function(GdObject, offsets) {
    ylim <- c(0, 1)
    h <- diff(ylim)
    middle <- mean(ylim)
    sh <- max(0, min(h, .dpOrDefault(GdObject, "stackHeight", 0.75)))
    space <- (h - (h * sh)) / 2
    shift <- switch(.dpOrDefault(GdObject, "stackJust", "middle"),
        "top" = space,
        "bottom" = -space,
        0
    )
    if (inherits(GdObject, "GeneRegionTrack")) {
        thinBox <- .dpOrDefault(GdObject, "thinBoxFeature", .THIN_BOX_FEATURES)
        space <- ifelse(feature(GdObject) %in% thinBox, space + ((middle -
            space) / 2), space)
    }
    shape <- .dpOrDefault(GdObject, "shape", "arrow")
    color <- .getBiotypeColor(GdObject)
    id <- identifier(GdObject, type = .dpOrDefault(
        GdObject, ifelse(is(GdObject, "GeneRegionTrack"), "exonAnnotation", "featureAnnotation"),
        "lowest"
    ))
    sel <- grepl("\\[Cluster_[0-9]*\\]", id)
    id[sel] <- sprintf(
        "%i merged\n%s", as.integer(.getAnn(GdObject, "density")[sel]),
        ifelse(class(GdObject) %in% c("AnnotationTrack", "DetailsAnnotationTrack"), "features", "exons")
    )
    yy1 <- ylim[1] + space + offsets + shift
    yy2 <- ylim[2] - space + offsets + shift
    boxes <- data.frame(
        cx1 = start(GdObject), cy1 = yy1, cx2 = start(GdObject) + width(GdObject), cy2 = yy2,
        fill = color, strand = strand(GdObject), text = id, textX = start(GdObject) + (width(GdObject) / 2), textY = middle + offsets,
        .getImageMap(cbind(start(GdObject), yy1, end(GdObject), yy2)),
        start = start(GdObject), end = end(GdObject), values(GdObject), exonId = id, origExonId = identifier(GdObject, type = "lowest"),
        stringsAsFactors = FALSE
    )
    rownames(boxes) <- if (.transcriptsAreCollapsed(GdObject)) {
        sprintf("uid%i", seq_along(identifier(GdObject)))
    } else {
        make.unique(identifier(GdObject, type = "lowest"))
    }
    return(boxes)
}

## Compute the coordinates for the bars connecting grouped items and the group labels

#' @importFrom grDevices col2rgb hsv rgb2hsv
.barsAndLabels <- function(GdObject) {
    bins <- stacks(GdObject)
    stacks <- max(bins)
    res <- .pxResolution(coord = "x")
    gp <- group(GdObject)
    if (!is.factor(gp)) {
        gp <- factor(gp, levels=unique(gp))
    }
    grpSplit <- split(range(GdObject), gp)
    grpRanges <- unlist(range(grpSplit))
    needBar <- vapply(grpSplit, length, FUN.VALUE = numeric(1L)) > 1 & width(grpRanges) > res
    ## If we draw the bar from start to end of the range we sometimes see little overlaps that extend beyond the first or last item.
    ## In order to fix this, we just subtract the equivalent of min.width pixels from both ends of each group range
    min.swidth <- res * .dpOrDefault(GdObject, "min.width", 2)
    nstart <- start(grpRanges[needBar]) + min.swidth
    nend <- end(grpRanges[needBar]) - min.swidth
    sel <- (nend - nstart) > 0
    start(grpRanges[needBar][sel]) <- nstart[sel]
    end(grpRanges[needBar][sel]) <- nend[sel]
    strand <- vapply(split(strand(GdObject), gp), function(x) {
        tmp <- unique(x)
        if (length(tmp) > 1) "*" else tmp
    }, FUN.VALUE = character(1L))
    yloc <- vapply(split((stacks - bins) + 1, gp), function(x) unique(x), FUN.VALUE = numeric(1L)) + 0.5
    color <- if (length(grep("__desatCol", values(GdObject)$feature[1]))) {
        .dpOrDefault(GdObject, "fill", .DEFAULT_FILL_COL)
    } else {
        vapply(split(.getBiotypeColor(GdObject), gp), head, 1, FUN.VALUE = character(1L))
    }
    bars <- data.frame(
        sx1 = start(grpRanges)[needBar], sx2 = end(grpRanges)[needBar], y = yloc[needBar], strand = strand[needBar],
        col = color[needBar], stringsAsFactors = FALSE
    )
    labs <- .dpOrDefault(GdObject, ".__groupLabels")[names(grpRanges)]
    if (!is.null(labs)) {
        lsel <- grepl("\\[Cluster_[0-9]*\\]", labs)
        if (any(lsel)) {
            gdens <- as.integer(vapply(split(.getAnn(GdObject, "gdensity"), gp), head, 1, FUN.VALUE = character(1L)))
            labs[lsel] <- sprintf(
                "%i merged %s  ", gdens[lsel],
                ifelse(class(GdObject) %in% c("AnnotationTrack", "DetailsAnnotationTrack"), "groups", "transcript models")
            )
        }
        just <- .dpOrDefault(GdObject, "just.group", "left")
        rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
        sizes <- .dpOrDefault(GdObject, ".__groupLabelWidths")[names(grpRanges), , drop = FALSE]
        pr <- .dpOrDefault(GdObject, ".__plottingRange", data.frame(from = min(start(GdObject)), to = max(end(GdObject))))
        ## grpRangesCut <- restrict(grpRanges, start=as.integer(pr["from"]), end=as.integer(pr["to"]))
        switch(just,
            "left" = {
                cx <- if (!rev) start(grpRanges) - sizes$after else end(grpRanges) + sizes$after
                cy <- yloc
                algn <- c("right", "center")
            },
            "right" = {
                cx <- if (!rev) end(grpRanges) + sizes$after else start(grpRanges) - sizes$after
                cy <- yloc
                algn <- c("left", "center")
            },
            "above" = {
                ## cx <- start(grpRangesCut) + width(grpRangesCut)/2
                ## indx <- which(seq_len(length(grpRanges)) %in% queryHits(findOverlaps(grpRanges, grpRangesCut)))
                ## cx <- cx[indx] # needs revision
                ## cy <- yloc[indx] + 0.5
                ## labs <- labs[indx]
                cx <- start(grpRanges) + width(grpRanges) / 2
                cy <- yloc + 0.5
                algn <- c("center", "top")
            },
            "below" = {
                ## cx <- start(grpRangesCut) + width(grpRangesCut)/2
                ## indx <- which(seq_len(length(grpRanges)) %in% queryHits(findOverlaps(grpRanges, grpRangesCut)))
                ## cx <- cx[indx] # needs revision
                ## cy <- yloc[indx] - 0.5
                ## labs <- labs[indx]
                cx <- start(grpRanges) + width(grpRanges) / 2
                cy <- yloc - 0.5
                algn <- c("center", "bottom")
            },
            stop(sprintf("Unknown label justification '%s'", just))
        )

        labels <- data.frame(txt = labs, x = cx, y = cy, stringsAsFactors = FALSE)
    } else {
        labels <- algn <- NA
    }
    return(list(bars = bars, labels = labels, align = algn))
}

#' @describeIn AnnotationTrack-class plot the object to a graphics device.
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
setMethod("drawGD", signature("AnnotationTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, subset = TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if ((is.logical(debug) && debug) || debug == "prepare") {
        browser()
    }
    imageMap(GdObject) <- NULL
    if (!length(GdObject)) {
        return(invisible(GdObject))
    }
    ## In prepare mode we need to make sure that the stacking information is updated from the optional display parameter (by calling
    ## the StackedTrack drawGD method) and also perform the collapsing of track items which could potentially lead to re-stacking.
    if (prepare) {
        GdObject <- callNextMethod(GdObject, ...)
        bins <- stacks(GdObject)
        stacks <- max(bins)
        ## We need to collapse the track object based on the current screen resolution (note that this may trigger re-stacking)
        pushViewport(dataViewport(xData = c(minBase, maxBase), extension = 0, yscale = c(1, stacks + 1), clip = TRUE))
        GdObject <- collapseTrack(GdObject, diff = .pxResolution(coord = "x"), xrange = c(minBase, maxBase))
        popViewport(1)
        return(invisible(GdObject))
    }
    if ((is.logical(debug) && debug) || debug == "draw") {
        browser()
    }
    ## If there are too many stacks for the available device resolution we cast an error, otherwise we set the viewport
    bins <- stacks(GdObject)
    stacks <- max(bins)
    rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
    xscale <- if (!rev) c(minBase, maxBase) else c(maxBase, minBase)
    yscale <- if (!.dpOrDefault(GdObject, "reverseStacking", FALSE)) c(1, stacks + 1) else c(stacks + 1, 1)
    pushViewport(dataViewport(xscale = xscale, extension = 0, yscale = yscale, clip = TRUE))
    res <- .pxResolution(coord = "x")
    curVp <- vpLocation()
    if (curVp$size["height"] / stacks < .dpOrDefault(GdObject, "min.height", 3)) {
        stop("Too many stacks to draw. Either increase the device size or limit the drawing to a smaller region.")
    }
    ## We adjust the color saturation to indicate overplotting if necessary
    if (.dpOrDefault(GdObject, "showOverplotting", FALSE)) {
        dens <- as.numeric(values(GdObject)$density)
        if (length(unique(dens)) != 1) {
            minSat <- max(0.25, 1 / max(dens))
            minDens <- min(dens)
            rDens <- diff(range(dens))
            saturation <- minSat + ((dens - minDens) / rDens / (1 / (1 - minSat)))
            bc <- unique(.getBiotypeColor(GdObject))
            baseCol <- rgb2hsv(col2rgb(bc))
            desatCols <- unlist(lapply(saturation, function(x) hsv(baseCol[1, ], x, baseCol[3, ])))
            names(desatCols) <- paste(unique(feature(GdObject)), rep(dens, each = length(bc)), sep = "_")
            feature(GdObject) <- paste(feature(GdObject), dens, sep = "_")
            desatCols <- desatCols[unique(names(desatCols))]
            displayPars(GdObject) <- as.list(desatCols)
        }
    }
    ## Now we can pre-compute all the coordinates and settings for the elements to be drawn, ...
    box <- .boxes(GdObject, (stacks - bins) + 1)
    barsAndLab <- .barsAndLabels(GdObject)
    bar <- barsAndLab$bars
    bartext <- barsAndLab$labels
    ## ... get all the necessary display parameters
    shape <- .dpOrDefault(GdObject, "shape", "arrow")
    col.line <- .dpOrDefault(GdObject, "col.line")[1]
    border <- .dpOrDefault(GdObject, "col")[1]
    if (is.null(border)) {
        border <- ifelse(is(GdObject, "GeneRegionTrack"), NA, "transparent")
    }
    lwd <- .dpOrDefault(GdObject, "lwd", 2)
    lty <- .dpOrDefault(GdObject, "lty", 1)
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
    rotation <- .dpOrDefaultFont(GdObject, "rotation", "item", 0)
    rotation.group <- .dpOrDefaultFont(GdObject, "rotation", "group", 0)
    just <- .dpOrDefault(GdObject, "just.group", "left")
    ## ... and finally draw whatever is needed
    if (nrow(box) > 0) {
        ## If we want to place a label on top or below the ranges we need to know how much size that will take up
        ## and adjust the plotting shapes accordingly
        drawLabel <- .dpOrDefault(GdObject, ".__hasAnno", FALSE) && !is.null(bartext) && !anyNA(bartext) &&
            nrow(bartext) > 0 && stacking(GdObject) != "dense"
        bs <- as.vector((stacks - bins) + 1)
        if (drawLabel && just %in% c("above", "below")) {
            labelHeights <- max(.getStringDims(GdObject, bartext$txt, subtype = "group")$height)
            labelSpace <- as.numeric(convertHeight(unit(2, "points"), "native"))
            avSpace <- (1 - max(box$cy2 - box$cy1)) / 2
            vadjust <- max(0, (labelHeights + (2 * labelSpace)) - avSpace)
            bh <- 1 - (avSpace * 2)
            bhNew <- bh - vadjust
            sfac <- bhNew / bh
            if (just == "above") {
                bar$y <- bar$y - (vadjust / 2)
                box$textY <- box$textY - (vadjust / 2)
                bartext$y <- bartext$y - labelSpace
                box$cy1 <- bs + ((box$cy1 %% 1 - avSpace) * sfac) + avSpace
                box$cy2 <- bs + ((box$cy2 %% 1 - avSpace) * sfac) + avSpace
            } else {
                bar$y <- bar$y + (vadjust / 2)
                box$textY <- box$textY + (vadjust / 2)
                bartext$y <- bartext$y + labelSpace
                box$cy1 <- bs + ((box$cy1 %% 1 - avSpace) * sfac) + avSpace + vadjust
                box$cy2 <- bs + ((box$cy2 %% 1 - avSpace) * sfac) + avSpace + vadjust
            }
        }
        ## Plotting of the (arrow)bar
        if (nrow(bar) > 0) {
            .arrowBar(bar$sx1, bar$sx2,
                y = bar$y, bar$strand, box[, seq_len(4), drop = FALSE],
                W = .dpOrDefault(GdObject, "arrowFeatherWidth", 3), D = .dpOrDefault(GdObject, "arrowFeatherDistance", 20),
                col = if (is.null(col.line)) bar$col else rep(col.line, length(bar$col)), lwd = lwd, lty = lty,
                alpha = alpha, barOnly = (!"smallArrow" %in% .dpOrDefault(GdObject, "shape", "box") || stacking(GdObject) == "dense"),
                diff = res, min.height = .dpOrDefault(GdObject, "min.height", 3)
            )
        }
        ## Plotting of the boxes
        box$col <- if (is.na(border)) box$fill else border
        if ("box" %in% shape || ("smallArrow" %in% shape && !("arrow" %in% shape || "ellipse" %in% shape))) {
            .filledBoxes(box, lwd = lwd, lty = lty, alpha = alpha)
        }
        ## Plotting of the ellipses
        if ("ellipse" %in% shape) {
            ellCoords <- .box2Ellipse(box)
            grid.polygon(
                x = ellCoords$x1, y = ellCoords$y1, id = ellCoords$id,
                gp = gpar(col = if (is.na(border)) box$fill else border, fill = box$fill, lwd = lwd, lty = lty, alpha = alpha),
                default.units = "native"
            )
        }
        ## Plotting of the filled arrows
        if ("arrow" %in% shape && !"box" %in% shape) {
            .filledArrow(box, lwd = lwd, lty = lty, alpha = alpha, min.width = 4 * res, max.width = .dpOrDefault(GdObject, "arrowHeadMaxWidth", 40) * res)
        }
        ## Plotting of the filled arrows with fixed head size
        if ("fixedArrow" %in% shape && !"box" %in% shape) {
            .filledArrow(box,
                lwd = lwd, lty = lty, alpha = alpha, min.width = 4 * res, absoluteWidth = TRUE,
                W = .dpOrDefault(GdObject, "arrowHeadWidth", 30) * res
            )
        }
        ## Plotting of the item labels
        if (.dpOrDefault(GdObject, "showFeatureId", FALSE)) {
            # grid.text(str2expression(box$text), box$textX, box$textY,
            grid.text(box$text, box$textX, box$textY,
                rot = rotation, gp = .fontGp(GdObject, subtype = "item"),
                default.units = "native", just = c("center", "center")
            )
        }
        ## Plotting of the group labels
        if (drawLabel) {
            grid.text(bartext$txt, bartext$x, bartext$y,
                rot = rotation.group, gp = .fontGp(GdObject, subtype = "group"),
                default.units = "native", just = barsAndLab$align
            )
        }
    }
    popViewport(1)
    ## Finally we set up the image map
    ## FIXME: we may want to record the merging information here
    im <- if (!is.null(box)) {
        coords <- as.matrix(box[, c("x1", "y1", "x2", "y2"), drop = FALSE])
        restCols <- setdiff(colnames(box), c("x1", "x2", "y1", "y2", "cx1", "cx2", "cy1", "cy2", "textX", "textY"))
        tags <- lapply(restCols, function(x) {
            tmp <- as.character(box[, x])
            names(tmp) <- rownames(coords)
            tmp
        })
        names(tags) <- restCols
        tags$title <- identifier(GdObject)
        ImageMap(coords = coords, tags = tags)
    } else {
        NULL
    }
    imageMap(GdObject) <- im
    return(invisible(GdObject))
})

## DrawGD DetailsAnnotationTrack ---------------------------------------------

## Create a data.frame with the distinct details function arguments (like start, end, ...)
.buildArgsDf <- function(GdObject) {
    groupDetails <- .dpOrDefault(GdObject, "groupDetails", FALSE)
    rr <- if (groupDetails) unlist(range(split(ranges(GdObject), group(GdObject)))) else ranges(GdObject)
    args <- data.frame(
        start = as.integer(start(rr)), end = as.integer(end(rr)), strand = as.character(strand(rr)),
        chromosome = as.character(seqnames(rr)),
        identifier = as.character(if (groupDetails) names(rr) else identifier(GdObject, type = "lowest")),
        stringsAsFactors = FALSE
    )
    return(args)
}

#' @describeIn AnnotationTrack-class plot the object to a graphics device.
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
setMethod(
    "drawGD", signature("DetailsAnnotationTrack"),
    function(GdObject, minBase, maxBase, prepare = FALSE, ...) {
        debug <- .dpOrDefault(GdObject, "debug", FALSE)
        if ((is.logical(debug) && debug) || debug == "prepare") {
            browser()
        }
        adf <- .buildArgsDf(GdObject)
        args <- .dpOrDefault(GdObject, "detailsFunArgs", fromPrototype = TRUE)
        groupDetails <- .dpOrDefault(GdObject, "groupDetails", FALSE)
        if (prepare) {
            GdObject <- callNextMethod()
            GdObject <- GdObject[order(start(GdObject))]
            indices <- if (groupDetails) seq_len(length(unique(group(GdObject)))) else seq_len(length(GdObject))
            xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
            pushViewport(viewport(xscale = xscale))
            select <- lapply(indices, function(i) {
                warn <- FALSE
                iargs <- as.list(adf[i, ])
                iargs$index <- i
                iargs$GdObject <- GdObject
                iargs$GdObject.original <- .dpOrDefault(GdObject, ".__OriginalGdObject", GdObject)
                args <- c(args[setdiff(names(args), names(iargs))], iargs)
                res <- do.call(GdObject@selectFun, args)
                if (length(res) != 1 || !is.logical(res) || is.na(res)) {
                    res <- TRUE
                    warn <- TRUE
                }
                c(res = res, warn = warn)
            })
            select <- do.call(rbind, select)
            if (any(select[, "warn"])) {
                warning("The result of function 'selectFun' has to be a single logical value. Forcing the value to 'TRUE'")
            }
            select <- select[, "res"]
            popViewport(1)
            displayPars(GdObject) <- list(".__select" = select)
            return(invisible(GdObject))
        }
        if ((is.logical(debug) && debug) || debug == "draw") {
            browser()
        }
        n <- length(GdObject)
        col <- rep(.dpOrDefault(GdObject, "detailsConnector.col", fromPrototype = TRUE), n)[seq_len(n)]
        lty <- rep(.dpOrDefault(GdObject, "detailsConnector.lty", fromPrototype = TRUE), n)[seq_len(n)]
        lwd <- rep(.dpOrDefault(GdObject, "detailsConnector.lwd", fromPrototype = TRUE), n)[seq_len(n)]
        pch <- rep(.dpOrDefault(GdObject, "detailsConnector.pch", fromPrototype = TRUE), n)[seq_len(n)]
        cex <- rep(.dpOrDefault(GdObject, "detailsConnector.cex", fromPrototype = TRUE), n)[seq_len(n)]
        border.lty <- rep(.dpOrDefault(GdObject, "detailsBorder.lty", fromPrototype = TRUE), n)[seq_len(n)]
        border.lwd <- rep(.dpOrDefault(GdObject, "detailsBorder.lwd", fromPrototype = TRUE), n)[seq_len(n)]
        border.col <- rep(.dpOrDefault(GdObject, "detailsBorder.col", fromPrototype = TRUE), n)[seq_len(n)]
        border.fill <- rep(.dpOrDefault(GdObject, "detailsBorder.fill", fromPrototype = TRUE), n)[seq_len(n)]
        minwidth <- .dpOrDefault(GdObject, "details.minWidth", fromPrototype = TRUE)
        size <- .dpOrDefault(GdObject, "details.size", fromPrototype = TRUE)
        xyratio <- .dpOrDefault(GdObject, "details.ratio", fromPrototype = TRUE)
        if (0 >= size || size > 1) {
            warning("details.size must be >0 and <1 - reset to 0.5")
            size <- 0.5
        }
        selection <- .dpOrDefault(GdObject, ".__select", rep(TRUE, length(GdObject)))
        len <- sum(selection)
        bins <- if (!groupDetails) stacks(GdObject) else unlist(lapply(split(stacks(GdObject), group(GdObject)), unique))
        stacks <- max(bins)
        if (len > 0) {
            if (((maxBase - minBase) / len) / .pxResolution(coord = "x") < minwidth) {
                warning("too much detail for available space (plot fewer annotation or increase details.minWidth)!")
                popViewport(1)
                GdObject <- callNextMethod()
                return(GdObject)
            }
            rr <- if (groupDetails) unlist(range(split(ranges(GdObject), group(GdObject)))) else ranges(GdObject)
            xloc1 <- (end(rr) - start(rr)) / 2 + start(rr)
            yloc1 <- (stacks - (bins - 0.5) + 1)
            xloc2 <- ((1 / len * seq_len(len)) - 1 / len + (1 / len * 0.5))
            yloc2 <- rep(1, len)
            ## draw details plots (via user supplied function 'fun')
            pushViewport(viewport(height = size, y = 1 - size, just = c(0.5, 0)))
            w <- 1
            v <- 0
            vpl <- vpLocation()
            r <- vpl$size["width"] / len / vpl$size["height"]
            if (r > xyratio) {
                w <- xyratio / r
                v <- ((1 / len) - (1 / len * w)) / 2
            }
            indices <- if (groupDetails) seq_len(length(unique(group(GdObject)))) else seq_len(length(GdObject))
            j <- 1
            pres <- list()
            hasError <- FALSE
            for (i in indices[selection]) {
                pushViewport(viewport(width = 1 / len * w, x = ((1 / len * j) - 1 / len) + (v), just = c(0, 0.5)))
                grid.rect(gp = gpar(col = border.col[i], lwd = border.lwd[i], lty = border.lty[i], fill = border.fill[i]))
                iargs <- as.list(adf[i, ])
                iargs$index <- i
                iargs$GdObject <- GdObject
                iargs$GdObject.original <- .dpOrDefault(GdObject, ".__OriginalGdObject", GdObject)
                args <- c(args[setdiff(names(args), names(iargs))], iargs)
                pres[[as.character(j)]] <- try(do.call(GdObject@fun, args), silent = TRUE)
                if (!is.null(pres) && is(pres[[as.character(j)]], "try-error")) {
                    hasError <- TRUE
                    grid.segments(x0 = c(0.1, 0.1), x1 = c(0.9, 0.9), y0 = c(0.9, 0.1), y1 = c(0.1, 0.9), gp = gpar(col = "red", lwd = 3))
                }
                popViewport(1)
                j <- j + 1
            }
            if (hasError) {
                warning("There have been errors in the detail plotting function:\n", paste(pres, collapse = "\n"))
            }
            popViewport(1)
            ## plot AnnotationTrack and connectors to details
            xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
            pushViewport(viewport(
                xscale = xscale,
                yscale = c(1, stacks + 1), clip = FALSE,
                height = 1 - size, y = 0, just = c(.5, 0)
            ))
            GdObject <- callNextMethod()
            grid.segments(
                x0 = unit(xloc1[selection], "native"), x1 = xloc2, y0 = unit(yloc1[selection], "native"),
                y1 = yloc2, gp = gpar(col = col, lwd = lwd, lty = lty, cex = cex)
            )
            grid.points(x = unit(xloc2, "npc"), y = unit(yloc2, "npc"), gp = gpar(col = col, cex = cex), pch = pch)
            grid.points(x = unit(xloc1[selection], "native"), y = unit(yloc1[selection], "native"), gp = gpar(col = col, cex = cex), pch = pch)
            popViewport(1)
        } else {
            GdObject <- callNextMethod()
        }
        return(invisible(GdObject))
    }
)

## SetAs ---------------------------------------------------------------------

#' @importFrom rtracklayer GenomicData
setAs(
    "AnnotationTrack", "UCSCData",
    function(from, to) {
        ranges <- range(from)
        dcolor <- as.integer(col2rgb(.dpOrDefault(from, "col")))
        line <- new("BasicTrackLine",
            name = names(from),
            description = names(from),
            visibility = stacking(from), color = dcolor, itemRgb = TRUE
        )
        vals <- values(from)
        color <- .getBiotypeColor(from)
        strand <- as.character(strand(from))
        strand[strand == "*"] <- "+"
        new("UCSCData", rtracklayer::GenomicData(ranges,
            chrom = chromosome(from),
            id = gsub(" ", "_", vals$id),
            name = gsub(" ", "_", as.character(vals$id)), itemRgb = color,
            strand = strand
        ),
        trackLine = line
        )
    }
)

setAs("GRanges", "AnnotationTrack", function(from, to) AnnotationTrack(range = from))

setAs("GRangesList", "AnnotationTrack", function(from, to) AnnotationTrack(range = from))

## Show ----------------------------------------------------------------------

## A helper function to plot general information about an AnnotationTrack
.annotationTrackInfo <- function(object) {
    msg <- sprintf(
        paste("| genome: %s\n| active chromosome: %s\n",
            "| annotation features: %s",
            sep = ""
        ),
        genome(object),
        chromosome(object),
        length(object)
    )
    addfeat <- length(object@range) - length(object)
    if (addfeat > 0) {
        msg <- c(msg, .addFeatInfo(object, addfeat), "Call chromosome(obj) <- 'chrId' to change the active chromosome")
    }
    return(paste(msg, collapse = "\n"))
}

## We have to show the name, genome and currently active chromosome, and, if more ranges are available on additional
## chromosomes some information about that

#' @describeIn AnnotationTrack-class Show method.
#' @export
setMethod("show", signature(object = "AnnotationTrack"), function(object) {
    cat(sprintf("AnnotationTrack '%s'\n%s\n", names(object), .annotationTrackInfo(object)))
})

#' @describeIn AnnotationTrack-class Show method.
#' @export
setMethod("show", signature(object = "ReferenceAnnotationTrack"), function(object) {
    .referenceTrackInfo(object, "ReferenceAnnotationTrack")
})
