#' @include GdObject-class.R
NULL

## RangeTrack Class ----------------------------------------------------------

#' RangeTrack class and methods
#'
#'
#' The virtual parent class for all track items in the Gviz package that
#' contain some form of genomic ranges (start, end, strand, chromosome
#' and the associated genome.)
#'
#'
#' @name RangeTrack-class
#'
#' @template GdObject-class_slot
#' @template RangeTrack-class_slot
#'
#' @template RangeTrack-class_param
#'
#' @return A virtual class: No objects may be created from it.
#'
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
#' @exportClass RangeTrack
setClass("RangeTrack",
    representation = representation("VIRTUAL",
        range = "GRangesOrIRanges",
        chromosome = "character",
        genome = "character"
    ),
    contains = "GdObject",
    prototype = prototype(
        chromosome = "chr1",
        dp = DisplayPars(),
        genome = "ANY",
        name = "RangeTrack",
        range = GRanges()
    )
)

## Coercing all input to the appropriate form
#' @describeIn RangeTrack-class Initialize.
#' @export
setMethod("initialize", "RangeTrack", function(.Object, range, chromosome, genome, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "RangeTrack")
    if (!missing(chromosome) && !is.null(chromosome)) {
        .Object@chromosome <- .chrName(chromosome)[1]
    }
    if (!missing(genome) && !is.null(genome)) {
        .Object@genome <- genome
    }
    if (!missing(range) && is(range, "GRanges")) {
        .Object@range <- range
    }
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})


## .buildRange ---------------------------------------------------------------
##
## Helper methods to build a GRanges object from the input arguments.
##
## Coordinates for grouped elements may be passed in as comma-separated values
## (e.g. "1,5,9"), in which case we need to split and convert to numeric.
## This also implies that the additional arguments (like feature, group, etc.)
## have to be replicated accordingly. We handle this by passing along the repeat
## vector 'by' to the numeric method below.

#' @importFrom Biobase rowMax
setMethod(
    ".buildRange", signature(
        "NULLOrMissing", "FactorOrCharacterOrNULL", "FactorOrCharacterOrNULL",
        "FactorOrCharacterOrNULL"
    ),
    function(range, start, end, width, asIRanges = FALSE, ...) {
        ## The inputs coordinates are all empty
        if (!length(start) && !length(end)) {
            return(if (asIRanges) IRanges() else GRanges())
        }
        delim <- ","
        coords <- c("start", "end", "width")
        lengths <- vapply(coords, function(x) length(get(x)), FUN.VALUE = numeric(1L))
        items <- structure(as.list(rep(0, 3)), names = coords)
        by <- NULL
        for (i in coords) {
            val <- get(i)
            if (!is.null(val)) {
                val <- strsplit(as.character(val), delim)
                lengths[i] <- length(val)
                items[[i]] <- listLen(val)
                val <- unlist(val)
                assign(i, as.numeric(val))
            }
        }
        len <- max(lengths)
        if (!all(unlist(lapply(items, function(x) length(x) == 1 && x == 0)))) {
            by <- Biobase::rowMax(do.call(cbind, lapply(items, function(x) {
                if (length(x) == 1) rep(x, len) else if (length(x)) x else rep(0, len)
            })))
        }
        return(.buildRange(start = start, end = end, width = width, asIRanges = asIRanges, by = by, len = len, ...))
    }
)

## Helper function to handle settings and defaults
.fillWithDefaults <- function(range = as.data.frame(matrix(ncol = 0, nrow = len)),
                              defaults, args, len, by = NULL, ignore = NULL) {
    for (a in setdiff(names(defaults), ignore)) {
        range[[a]] <- if (is.null(args[[a]])) {
            if (is.null(defaults[[a]])) {
                stop("The mandatory argument '", a, "' is missing with no default")
            }
            val <- defaults[[a]]
            if (length(val) == 1) {
                val <- rep(val, len)
            }
            if (!length(val) %in% c(len, sum(by))) {
                stop("Number of elements in '", a, "' is invalid")
            }
            if (!is.null(by) && length(val) != sum(by)) rep(val, by) else val
        } else {
            val <- args[[a]]
            if (length(val) == 1) {
                val <- rep(val, len)
            }
            if (!length(val) %in% c(len, sum(by))) {
                stop("Number of elements in argument '", a, "' is invalid")
            }
            if (!is.null(by) && length(val) != sum(by)) rep(val, by) else val
        }
    }
    return(range)
}

## For numeric vectors we can immediately create a data frame after some sanity checking and pass that on to the next method.
setMethod(
    ".buildRange", signature("NULLOrMissing", "NumericOrNULL", "NumericOrNULL", "NumericOrNULL"),
    function(range, start, end, width, asIRanges = FALSE, by = NULL, len, args, defaults, ...) {
        ## Some of the arguments are mutually exclusive and we want to catch this here.
        if (is.null(width)) {
            ## The inputs coordinates are all empty
            ## CHECK if this might solve the issue on windows
            if (is.null(start) && is.null(end)) {
                return(if (asIRanges) IRanges() else GRanges())
            }
            if (is.null(start) || is.null(end)) {
                stop("Must specify either start and end or width")
            }
        } else {
            if (is.null(end) && !is.null(start)) {
                end <- as.integer(start) + as.integer(width)
            } else if (is.null(start) && !is.null(end)) {
                start <- as.integer(end) - as.integer(width)
            } else {
                stop("Can't pass all three of 'start', 'end' and 'width'")
            }
        }
        if (length(start) != 1 && length(end) != 1 && length(start) != length(end)) {
            stop("Start and end must be vectors of the same length")
        }
        if (missing(len)) {
            len <- length(start)
        }
        range <- data.frame()
        if (length(start) > 0) {
            if (asIRanges) {
                return(IRanges(start = as.integer(start), end = as.integer(end)))
            }
            range <- .fillWithDefaults(data.frame(start = as.integer(start), end = as.integer(end)), defaults, args, len, by)
        }
        return(.buildRange(range = range, asIRanges = asIRanges, args = args["genome"], defaults = defaults, ...))
    }
)

## For data.frames we need to check for additional arguments
## (like feature, group, etc.), the chromosome information
## and create the final GRanges object
setMethod(
    ".buildRange", signature("data.frame"),
    function(range, asIRanges = FALSE, args = list(), defaults = list(), chromosome = NULL, trackType, ...) {
        if (asIRanges) {
            range <- .fillWithDefaults(range, defaults, args, len = nrow(range), ignore = setdiff(names(defaults), c("start", "end", "genome")))
            return(IRanges(start = as.integer(range$start), end = as.integer(range$end)))
        }
        mandArgs <- c("start", "end", "genome", names(defaults))
        ## Not quite sure how whether existing chromosome information in a GRanges object should generally have precedence over the
        ## chromosome constructor, but probably that should be the case
        if ("chromosome" %in% colnames(range)) {
            args$chromosome <- NULL
        }
        missing <- setdiff(union(setdiff(mandArgs, c(colnames(range))), names(which(!vapply(args, is.null, FUN.VALUE = logical(1L))))), "genome")
        range <- .fillWithDefaults(range, defaults[missing], args[missing], len = nrow(range))
        range$chromosome <- .chrName(as.character(range$chromosome))
        grange <- GRanges(ranges = IRanges(start = range$start, end = range$end), strand = range$strand, seqnames = range$chromosome)
        mcols(grange) <- range[, setdiff(colnames(range), c(
            "start", "end", "strand", "width", "chromosome", "genome", "seqnames",
            "ranges", "seqlevels", "seqlengths", "isCircular", "element"
        ))]
        if (trackType != "DataTrack") {
            mcols(grange) <- mcols(grange)[, intersect(names(defaults), colnames(mcols(grange)))]
        }
        suppressWarnings(genome(grange) <- unname(if (is.null(args[["genome"]])) defaults[["genome"]] else as.character(args[["genome"]])[[1]]))
        return(grange)
    }
)


## For GRanges we just need to check for the existence of additional
## arguments (like feature, group, etc.)
setMethod(
    ".buildRange", signature("GRanges"),
    function(range, asIRanges = FALSE, args = list(), defaults = list(), trackType = NULL, ...) {
        if (asIRanges) {
            return(ranges(range))
        }
        if (length(range)) {
            mandArgs <- names(defaults)
            ## Not quite sure how whether existing chromosome information in a GRanges object should generally have precedence over the
            ## chromosome constructor, but probably that should be the case
            args$chromosome <- NULL
            range <- renameSeqlevels(range, setNames(.chrName(seqlevels(range)), seqlevels(range)))
            missing <- setdiff(union(setdiff(mandArgs, c("chromosome", "strand", colnames(mcols(range)))), names(which(!vapply(args, is.null, FUN.VALUE = logical(1L))))), "genome")
            newVars <- .fillWithDefaults(DataFrame(chromosome = as.character(seqnames(range)), strand = as.character(strand(range)), mcols(range), check.names = FALSE),
                defaults[missing], args[missing],
                len = length(range), ignore = c("flag")
            )
            if (any(c("start", "end", "strand", "chromosome") %in% colnames(newVars))) {
                gen <- genome(range)
                range <- GRanges(
                    seqnames = if (is.null(newVars[["chromosome"]])) seqnames(range) else (newVars[["chromosome"]]),
                    strand = if (is.null(newVars[["strand"]])) strand(range) else (newVars[["strand"]]),
                    ranges = IRanges(
                        start = if (is.null(newVars[["start"]])) start(range) else (newVars[["start"]]),
                        end = if (is.null(newVars[["end"]])) end(range) else (newVars[["end"]])
                    )
                )
                if (length(unique(gen)) != 1) {
                    warning("Tracks can only be defined for a single genome. Forcing all reads to belong to genome '", gen[1], "'")
                }
                defaults[["genome"]] <- as.character(gen)[1]
            }
            mcols(range) <- newVars[, setdiff(colnames(newVars), c(
                "start", "end", "strand", "width", "chromosome", "genome", "seqnames",
                "ranges", "seqlevels", "seqlengths", "isCircular", "element"
            )), drop = FALSE]
        }
        if (trackType != "DataTrack") {
            mcols(range) <- mcols(range)[, intersect(names(defaults), colnames(mcols(range)))]
        }
        ## The genome information may or may not be encoded in the GRanges object at this time but we want it in there for sure
        genome <- if (!is.null(args[["genome"]])) args[["genome"]] else .getGenomeFromGRange(range, defaults[["genome"]])
        suppressWarnings(genome(range) <- unname(genome))[1]
        return(range)
    }
)

## For IRanges we need to deal with additional arguments
## (like feature, group, etc.) and create the final GRanges object
setMethod(
    ".buildRange", signature("IRanges"),
    function(range, asIRanges = FALSE, args = list(), defaults = list(), chromosome = NULL, strand, ...) {
        if (asIRanges) {
            return(range)
        }
        if (missing(chromosome) || is.null(chromosome)) {
            stop("Unable to find chromosome information in any of the arguments")
        }
        range <- GRanges(seqnames = .chrName(chromosome), ranges = range, strand = if (!is.null(args$strand)) args$strand else "*")
        if (length(range)) {
            vals <- .fillWithDefaults(defaults = defaults, args = args, len = (length(range)), by = NULL, ignore = "strand")
            mcols(range) <- vals
        }
        return(range)
    }
)

## For GRangesLists we capture the grouping information from the list
## structure, `unlist` and use the `GRanges` method
setMethod(
    ".buildRange", signature("GRangesList"),
    function(range, groupId = "group", ...) {
        grps <- rep(names(range), elementNROWS(range))
        range <- unlist(range)
        names(range) <- NULL
        mcols(range)[[groupId]] <- grps
        return(.buildRange(range = range, ...))
    }
)



## RangeTrack Methods ranges, range ------------------------------------------

## Extract the full GRanges object from the range slot
## of an object inheriting from RangeTrack

#' @describeIn RangeTrack-class return the genomic coordinates for the track
#' along with all additional annotation information as an object of
#' class `GRanges.`
setMethod("ranges", "RangeTrack", function(x) x@range)
setReplaceMethod("ranges", "RangeTrack", function(x, value) {
    x@range <- value
    return(x)
})

## Extract the IRanges part of the GRanges object from the range slot
## of an object inheriting from RangeTrack

#' @describeIn RangeTrack-class return the genomic coordinates for the
#' track as an object of class IRanges.
#' @export
setMethod("range", "RangeTrack", function(x) ranges(x@range))

## RangeTrack Methods seqnames, seqlevels and seqinfo, genome ----------------

#' @describeIn RangeTrack-class return the track's seqnames.
#' @export
setMethod("seqnames", "RangeTrack", function(x) as.character(seqnames(ranges(x))))

#' @describeIn RangeTrack-class return the track's seqlevels.
#' @export
setMethod("seqlevels", "RangeTrack", function(x) unique(seqnames(x)))

#' @describeIn RangeTrack-class return the track's seqinfo.
#' @export
setMethod("seqinfo", "RangeTrack", function(x) table(seqnames(x)))

#' @describeIn RangeTrack-class return the track's genome.
#' @export
setMethod("genome", "RangeTrack", function(x) x@genome)

#' @describeIn RangeTrack-class set the track's genome. Usually this has to
#' be a valid UCSC identifier, however this is not formally enforced here.
#' @export
setReplaceMethod("genome", "RangeTrack", function(x, value) {
    x@genome <- value[1]
    genome(ranges(x)) <- as.vector(value[1])
    return(x)
})

#' @describeIn RangeTrack-class return the chromosome for which the track is defined.
#' @export
setMethod("chromosome", "RangeTrack", function(GdObject) GdObject@chromosome)

#' @describeIn RangeTrack-class replace the value of the track's chromosome.
#' This has to be a valid UCSC chromosome identifier or an integer or character
#' scalar that can be reasonably coerced into one.
#' @export
setReplaceMethod("chromosome", "RangeTrack", function(GdObject, value) {
    GdObject@chromosome <- .chrName(value[1])
    return(GdObject)
})
## RangeTrack Methods start, end, min and max ranges -------------------------

#' @describeIn RangeTrack-class the start of the track items in genomic
#' coordinates.
#' @export
setMethod("start", "RangeTrack", function(x) if (length(x)) as.integer(start(range(x))) else NULL)

#' @describeIn RangeTrack-class replace the start of the track items in
#' genomic coordinates.
#' @export
setReplaceMethod("start", "RangeTrack", function(x, value) {
    start(x@range) <- value
    return(x)
})

#' @describeIn RangeTrack-class the end of the track items in genomic
#' coordinates.
#' @export
setMethod("end", "RangeTrack", function(x) if (length(x)) as.integer(end(range(x))) else NULL)

#' @describeIn RangeTrack-class replace the end of the track items in
#' genomic coordinates.
#' @export
setReplaceMethod("end", "RangeTrack", function(x, value) {
    end(x@range) <- value
    return(x)
})

#' @describeIn RangeTrack-class the width of the track items in genomic
#' coordinates.
#' @export
setMethod("width", "RangeTrack", function(x) if (length(x)) as.integer(width(range(x))) else NULL)

#' @describeIn RangeTrack-class replace the width of the track items in genomic
#' coordinates.
#' @export
setReplaceMethod("width", "RangeTrack", function(x, value) {
    width(x@range) <- value
    return(x)
})

#' @describeIn RangeTrack-class return the start position for the leftmost range item.
#' @export
setMethod("min", "RangeTrack", function(x) min(start(x)))

#' @describeIn RangeTrack-class return the end position for the rightmost range item.
#' @export
setMethod("max", "RangeTrack", function(x) max(end(x)))

#' @describeIn RangeTrack-class return the number of items in the track.
#' @export
setMethod("length", "RangeTrack", function(x) sum(seqnames(x) == chromosome(x)))

#' @describeIn RangeTrack-class  return a vector of strand specifiers for all
#' track items, in the form '+' for the Watson strand, '-' for the Crick
#' strand or '*' for either of the two.
#' @export
setMethod("strand", "RangeTrack", function(x) as.character(strand(ranges(x))))

#' @describeIn RangeTrack-class replace the strand information for the track
#' items. The replacement value needs to be an appropriate scalar or vector
#' of strand values.
#' @export
setReplaceMethod("strand", "RangeTrack", function(x, value) {
    if ((length(value) != 1 && length(value) != length(x)) || !all(value %in% c("+", "-", "*"))) {
        stop(
            "Invalid replacement value or length of replacement value ",
            "for the strand information does not match the ",
            "number of items in the track"
        )
    }
    r <- ranges(x)
    strand(r) <- value
    x@range <- r
    return(x)
})

#' @describeIn RangeTrack-class the arithmetic mean of the track item's
#' coordionates, i.e., `(end(obj)-start(obj))/2`.
#' @export
setMethod("position", signature("RangeTrack"), definition = function(GdObject, from = NULL, to = NULL, sort = FALSE, ...) {
    if (!is.null(from) && !is.null(to)) {
        GdObject <- subset(GdObject, from = from, to = to, sort = sort, ...)
    }
    pos <- if (length(GdObject)) rowMeans(cbind(start(GdObject), end(GdObject))) else numeric()
    return(pos)
})

## RangeTrack Methods subsetting, split---------------------------------------

#' @describeIn RangeTrack-class subset the items in the `RangeTrack` object.
#' This is essentially similar to subsetting of the `GRanges` object in the
#' `range` slot. For most applications, the subset method may be more appropriate.
#' @export
setMethod("[", signature(x = "RangeTrack"), function(x, i, j, ..., drop = TRUE) {
    x <- .deepCopyPars(x)
    x@range <- x@range[i, , drop = drop]
    return(x)
})

#' @describeIn RangeTrack-class subset a `RangeTrack` by coordinates and
#' sort if necessary.
#' @export
setMethod("subset", signature(x = "RangeTrack"), function(x, from = NULL, to = NULL, sort = FALSE, drop = TRUE, use.defaults = TRUE, ...) {
    ## Not needed anymore...
    ## Subset to a single chromosome first
    if (drop) {
        csel <- seqnames(x) != chromosome(x)
        if (any(csel)) {
            x <- x[, !csel]
        }
    }
    if (!length(x)) {
        return(x)
    }
    ranges <- if (use.defaults) .defaultRange(x, from = from, to = to) else c(from = ifelse(is.null(from), -Inf, from), to = ifelse(is.null(to), Inf, to))
    lsel <- end(x) < ranges["from"]
    if (any(lsel)) {
        lsel[max(0, max(which(lsel)) - 1)] <- FALSE
    }
    rsel <- start(x) > ranges["to"]
    if (any(rsel)) {
        rsel[min(length(x), min(which(rsel)) + 1)] <- FALSE
    }
    if (any(lsel) || any(rsel)) {
        x <- x[!(lsel | rsel), ]
    }
    if (sort) {
        x <- x[order(range(x)), ]
    }
    return(x)
})



#' @describeIn RangeTrack-class split a `RangeTrack` object by an appropriate
#' factor vector (or another vector that can be coerced into one). The output
#' of this operation is a list of objects of the same class as the input
#' object, all inheriting from class `RangeTrack.`
#' @export
setMethod("split", signature("RangeTrack"),
    definition = function(x, f, ...) {
        rs <- split(ranges(x), factor(f))
        lapply(rs, function(y) {
            x@range <- y
            return(x)
        })
    }
)

## RangeTrack Methods value, feature -----------------------------------------

#' @describeIn RangeTrack-class return all additional annotation information
#' except for the genomic coordinates for the track items as a `data.frame`.
#' @export
setMethod("values", "RangeTrack", function(x) as.data.frame(values(ranges(x))))

#' @describeIn RangeTrack-class return the grouping information for track
#' items. For certain sub-classes, groups may be indicated by different colour
#'  schemes when plotting. See grouping or `AnnotationTrack` and
#'  `GeneRegionTrack` for details.
#'  @export
setMethod("feature", signature(GdObject = "RangeTrack"), function(GdObject) .getAnn(GdObject, "feature"))

#' @describeIn RangeTrack-class set the grouping information for track items.
#' This has to be a factor vector (or another type of vector that can be
#' coerced into one) of the same length as the number of items in the
#' `RangeTrack.` See grouping or `AnnotationTrack` and `GeneRegionTrack` for
#'  details.
#'  @export
setReplaceMethod("feature", signature("RangeTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "feature"))


## RangeTrack Methods consolidate --------------------------------------------

#' @describeIn RangeTrack-class Consolidate.
#' @param GdObject the input track object
#' @param chromosome the currently active chromosome which may have to be set
#' for a `RangeTrack` or a `SequenceTrack` object
#' parameters
# #' @keywords internal
#' @export
setMethod("consolidateTrack", signature(GdObject = "RangeTrack"), function(GdObject, chromosome, ...) {
    if (!is.null(chromosome)) {
        chromosome(GdObject) <- chromosome
    }
    GdObject <- callNextMethod(GdObject, ...)
    return(GdObject)
})

## RangeTrack SetAs ----------------------------------------------------------


#' @noRd
#' @keywords internal
setAs(
    "RangeTrack", "data.frame",
    function(from, to) as(as(ranges(from), "DataFrame"), "data.frame")
)
