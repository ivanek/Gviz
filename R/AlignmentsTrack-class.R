#' @include DataTrack-class.R
#' @include ReferenceTrack-class.R
#' @include StackedTrack-class.R
#' @include SequenceTrack-class.R
NULL


setClassUnion("SequenceTrackOrNULL", c("SequenceTrack", "NULL"))

## AlignmentsTrack Class -----------------------------------------------------

#' AlignmentsTrack class and methods
#'
#'
#' A class to represent short sequences that have been aligned to a reference
#' genome as they are typically generated in next generation sequencing
#' experiments.
#'
#' @name AlignmentsTrack-class
#'
#' @param range An optional meta argument to handle the different input types.
#' If the \code{range} argument is missing, all the relevant information to
#' create the object has to be provided as individual function arguments
#' (see below).
#'
#' The different input options for \code{range} are:
#'
#' \describe{
#'
#' \item{A \code{character} string:}{ the path to a \code{BAM} file containing
#' the read alignments. To be precise, this will result in the instantiation of
#' a \code{ReferenceAlignmentsTrack} object, but for the user this
#' implementation detail should be of no concern.}
#'
#' \item{A \code{GRanges} object:}{ the genomic ranges of the individual reads
#' as well as the optional additional metadata columns \code{id}, \code{cigar},
#' \code{mapq}, \code{flag}, \code{isize}, \code{groupid}, \code{status},
#' \code{md} and \code{seqs} (see description of the individual function
#' parameters below for details). Calling the constructor on a \code{GRanges}
#' object without further arguments, e.g. \code{AlignmentsTrack(range=obj)} is
#' equivalent to calling the coerce method \code{as(obj, "AlignmentsTrack")}.}
#'
#' \item{An \code{\linkS4class{IRanges}} object:}{ almost identical to the
#' \code{GRanges} case, except that the chromosome and strand information as
#' well as all additional metadata has to be provided in the separate
#' \code{chromosome}, \code{strand}, \code{feature}, \code{group} or \code{id}
#' arguments, because it can not be directly encoded in an \code{IRanges}
#' object. Note that none of those inputs are mandatory, and if not provided
#' explicitely the more or less reasonable default values \code{chromosome=NA}
#' and \code{strand="*"} are used. }
#'
#' \item{A \code{data.frame} object:}{ the \code{data.frame} needs to contain
#' at least the two mandatory columns \code{start} and \code{end} with the
#' range coordinates. It may also contain a \code{chromosome} and a
#' \code{strand} column with the chromosome and strand information for each
#' range. If missing it will be drawn from the separate \code{chromosome} or
#' \code{strand} arguments. In addition, the \code{id}, \code{cigar},
#' \code{mapq}, \code{flag}, \code{isize}, \code{groupid}, \code{status},
#' \code{md} and \code{seqs} data can be provided as additional columns. The
#' above comments about potential default values also apply here.}
#'
#' }
#'
#' @template AlignmentsTrack-class_param
#'
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{AlignmentsTrack} or \code{ReferenceAlignmentsTrack}.
#' @section Objects from the Class:
#'
#' Objects can be created using the constructor function
#' \code{AlignmentsTrack}.
#'
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @examples
#' ## Creating objects
#' afrom <- 2960000
#' ato <- 3160000
#' alTrack <- AlignmentsTrack(system.file(
#'     package = "Gviz", "extdata",
#'     "gapped.bam"
#' ), isPaired = TRUE)
#' plotTracks(alTrack, from = afrom, to = ato, chromosome = "chr12")
#'
#' ## Omit the coverage or the pile-ups part
#' plotTracks(alTrack,
#'     from = afrom, to = ato, chromosome = "chr12",
#'     type = "coverage"
#' )
#' plotTracks(alTrack,
#'     from = afrom, to = ato, chromosome = "chr12",
#'     type = "pileup"
#' )
#'
#' ## Including sequence information with the constructor
#' if (require(BSgenome.Hsapiens.UCSC.hg19)) {
#'     strack <- SequenceTrack(Hsapiens, chromosome = "chr21")
#'     afrom <- 44945200
#'     ato <- 44947200
#'     alTrack <- AlignmentsTrack(system.file(
#'         package = "Gviz", "extdata",
#'         "snps.bam"
#'     ), isPaired = TRUE, referenceSequence = strack)
#'     plotTracks(alTrack, chromosome = "chr21", from = afrom, to = ato)
#'
#'     ## Including sequence information in the track list
#'     alTrack <- AlignmentsTrack(system.file(
#'         package = "Gviz", "extdata",
#'         "snps.bam"
#'     ), isPaired = TRUE)
#'     plotTracks(c(alTrack, strack),
#'         chromosome = "chr21", from = 44946590,
#'         to = 44946660
#'     )
#' }
#' @importFrom Rsamtools scanBamFlag scanBamHeader scanBam ScanBamParam scanFaIndex scanFa BamFile scanBamWhat bamWhich
#'
#' @exportClass AlignmentsTrack
setClass("AlignmentsTrack",
    representation = representation(
        stackRanges = "GRanges",
        sequences = "DNAStringSet",
        referenceSequence = "SequenceTrackOrNULL"
    ),
    contains = "StackedTrack",
    prototype = prototype(
        stacking = "squish",
        name = "AlignmentsTrack",
        coverageOnly = FALSE,
        stackRanges = GRanges(),
        sequences = Biostrings::DNAStringSet(),
        referenceSequence = NULL,
        dp = DisplayPars(
            alpha.reads = 0.5,
            alpha.mismatch = 1,
            cex = 0.7,
            cex.mismatch = NULL,
            col.coverage = NULL,
            col.gap = .DEFAULT_SHADED_COL,
            col.mates = .DEFAULT_BRIGHT_SHADED_COL,
            col.deletion = "#000000",
            col.insertion = "#984EA3",
            col.mismatch = .DEFAULT_SHADED_COL,
            col.reads = NULL,
            col.sashimi = NULL,
            col = .DEFAULT_SHADED_COL,
            collapse = FALSE,
            coverageHeight = 0.1,
            fill.coverage = NULL,
            fill.reads = NULL,
            fill = "#BABABA",
            fontface.mismatch = 2,
            lty.coverage = NULL,
            lty.gap = NULL,
            lty.mates = NULL,
            lty.deletion = NULL,
            lty.insertion = NULL,
            lty.mismatch = NULL,
            lty.reads = NULL,
            lty = 1,
            lwd.coverage = NULL,
            lwd.gap = NULL,
            lwd.mates = NULL,
            lwd.deletion = NULL,
            lwd.insertion = NULL,
            lwd.mismatch = NULL,
            lwd.reads = NULL,
            lwd.sashimiMax = 10,
            lwd = 1,
            max.height = 10,
            min.height = 5,
            minCoverageHeight = 50,
            minSashimiHeight = 50,
            noLetters = FALSE,
            sashimiFilter = NULL,
            sashimiFilterTolerance = 0L,
            sashimiHeight = 0.1,
            sashimiScore = 1,
            sashimiStrand = "*",
            sashimiTransformation = NULL,
            showIndels = FALSE,
            showMismatches = TRUE,
            size = NULL,
            transformation = NULL,
            type = c("coverage", "pileup")
        )
    )
)

## Initialize ----------------------------------------------------------------

#' @describeIn AlignmentsTrack-class Initialize.
#' @export
setMethod("initialize", "AlignmentsTrack", function(.Object, stackRanges = GRanges(), stacks = numeric(), sequences = DNAStringSet(),
                                                    referenceSequence = NULL, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "AlignmentsTrack")
    .Object@stackRanges <- stackRanges
    .Object <- callNextMethod(.Object, ...)
    .Object@stacks <- stacks
    .Object@sequences <- sequences
    .Object@referenceSequence <- referenceSequence
    return(.Object)
})

## ReferenceAlignmentsTrack Class --------------------------------------------

#' @describeIn AlignmentsTrack-class The file-based version of the `AlignmentsTrack-class`.
#' @exportClass ReferenceAlignmentsTrack
setClass("ReferenceAlignmentsTrack", contains = c("AlignmentsTrack", "ReferenceTrack"))

## Initialize ----------------------------------------------------------------

## This just needs to set the appropriate slots that are being inherited
## from ReferenceTrack because the multiple inheritance has some strange
## features with regards to method selection

#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom Biostrings DNAStringSet
#' @describeIn AlignmentsTrack-class Initialize.
#' @export
setMethod("initialize", "ReferenceAlignmentsTrack", function(.Object, stream, reference, mapping = list(),
                                                             args = list(), defaults = list(), stacks = numeric(),
                                                             stackRanges = GRanges(), sequences = Biostrings::DNAStringSet(),
                                                             referenceSequence = NULL, ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object = .Object, reference = reference, stream = stream,
        mapping = mapping, args = args, defaults = defaults)
    .Object <- callNextMethod(.Object, ...)
    .Object@referenceSequence <- referenceSequence
    return(.Object)
})

## Constructor ---------------------------------------------------------------

## Constructor
#' @describeIn AlignmentsTrack-class Constructor for `AlignmentsTrack-class`.
#' @export
AlignmentsTrack <- function(range = NULL, start = NULL, end = NULL, width = NULL, strand, chromosome, genome,
                            stacking = "squish", id, cigar, mapq, flag = scanBamFlag(isUnmappedQuery = FALSE), isize, groupid, status, md, seqs,
                            name = "AlignmentsTrack", isPaired = TRUE, importFunction, referenceSequence, ...) {
    ## Some defaults
    if (missing(importFunction)) {
        importFunction <- .import.bam.alignments
    }
    covars <- .getCovars(range)
    isStream <- FALSE
    if (!is.character(range)) {
        n <- max(c(length(start), length(end), length(width)), nrow(covars))
        id <- .covDefault(id, covars[["id"]], paste("read", seq_len(n), sep = "_"))
        cigar <- .covDefault(cigar, covars[["cigar"]], paste(if (is(range, "GRangesOrIRanges")) width(range) else width, "M", sep = ""))
        mapq <- .covDefault(mapq, covars[["mapq"]], rep(as.integer(NA), n))
        flag <- .covDefault(flag, covars[["flag"]], rep(as.integer(NA), n))
        isize <- .covDefault(isize, covars[["isize"]], rep(as.integer(NA), n))
        groupid <- .covDefault(groupid, covars[["groupid"]], seq_len(n))
        md <- .covDefault(md, covars[["md"]], rep(as.character(NA), n))
        status <- .covDefault(status, covars[["status"]], ifelse(groupid %in% groupid[duplicated(groupid)], "mated", "unmated"))
    }
    ## Build a GRanges object from the inputs
    .missingToNull(c(
        "strand", "chromosome", "importFunction", "genome", "id", "cigar", "mapq", "flag", "isize", "groupid", "status",
        "md", "seqs", "referenceSequence"
    ))
    args <- list(
        id = id, cigar = cigar, mapq = mapq, flag = flag, isize = isize, groupid = groupid, status = status, strand = strand, md = md,
        chromosome = chromosome, genome = genome
    )
    defs <- list(
        strand = "*", chromosome = "chrNA", genome = NA, id = as.character(NA), cigar = as.character(NA), mapq = as.integer(NA),
        flag = as.integer(NA), isize = as.integer(NA), groupid = as.character(NA), status = as.character(NA), md = as.character(NA)
    )
    range <- .buildRange(
        range = range, start = start, end = end, width = width,
        args = args, defaults = defs, chromosome = chromosome, trackType = "AlignmentsTrack",
        importFun = importFunction, stream = TRUE, autodetect = TRUE
    )
    ## This is going to be a list if we have to stream data from a file, otherwise we can compute some additional values
    if (is.list(range)) {
        isStream <- TRUE
        slist <- range
        range <- GRanges()
        stackRanges <- GRanges()
        stacks <- NULL
        seqs <- DNAStringSet()
    } else {
        if (is.null(seqs)) {
            seqs <- DNAStringSet(vapply(width(range), function(x) paste(rep("N", x), collapse = ""), character(1)))
        }
        addArgs <- list(...)
        if ("showIndels" %in% names(addArgs)) {
            showIndels <- addArgs$showIndels
        } else {
            showIndels <- FALSE
        }
        tmp <- .computeAlignments(range, drop.D.ranges = showIndels)
        range <- tmp$range
        stackRanges <- tmp$stackRange
        stacks <- tmp$stacks
    }
    ## If no chromosome was explicitly asked for we just take the first one in the GRanges object
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- if (length(range) > 0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    }
    ## And finally the object instantiation
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    if (!isStream) {
        return(new("AlignmentsTrack",
            chromosome = chromosome[1], range = range, stacks = stacks,
            name = name, genome = genome, stacking = stacking, stackRanges = stackRanges, sequences = seqs,
            referenceSequence = referenceSequence, ...
        ))
    } else {
        ## A bit hackish but for some functions we may want to know which track type we need but at the
        ## same time we do not want to enforce this as an additional argument
        e <- new.env()
        e[["._trackType"]] <- "AlignmentsTrack"
        e[["._isPaired"]] <- isPaired
        e[["._flag"]] <- flag
        environment(slist[["stream"]]) <- e
        return(new("ReferenceAlignmentsTrack",
            chromosome = chromosome[1], range = range, stackRanges = stackRanges,
            name = name, genome = genome, stacking = stacking, stream = slist[["stream"]], reference = slist[["reference"]],
            mapping = slist[["mapping"]], args = args, defaults = defs, stacks = stacks, referenceSequence = referenceSequence, ...
        ))
    }
}

## General accessors ---------------------------------------------------------

#' @describeIn AlignmentsTrack-class Return all additional annotation
#' information except for the genomic coordinates for the track items as
#' a `data.frame`.
#' @export
setMethod("values", "AlignmentsTrack", function(x) .dpOrDefault(x, ".__coverage"))

#' @describeIn AlignmentsTrack-class replace the value of the track's chromosome.
#' This has to be a valid UCSC chromosome identifier or an integer or character
#' scalar that can be reasonably coerced into one.
#' @export
setReplaceMethod("chromosome", "AlignmentsTrack", function(GdObject, value) {
    GdObject <- callNextMethod()
    if (!is.null(GdObject@referenceSequence)) {
        chromosome(GdObject@referenceSequence) <- value[1]
    }
    return(GdObject)
})



## Annotation Accessors ------------------------------------------------------
## Stacking ------------------------------------------------------------------

#' @describeIn AlignmentsTrack-class return the stack indices for each track item.
#' @export
setMethod(
    "stacks", "AlignmentsTrack",
    function(GdObject) if (length(GdObject)) ranges(GdObject)$stack else 0
)

#' @describeIn AlignmentsTrack-class recompute the stacks based on the available
#' space and on the object's track items and stacking settings.
#' @export
setMethod("setStacks", "AlignmentsTrack", function(GdObject, ...) {
    if (length(GdObject)) {
        bins <- if (!.needsStacking(GdObject)) rep(1, length(GdObject)) else disjointBins(GdObject@stackRanges)
        GdObject@stacks <- bins
        ranges(GdObject)$stack <- bins[match(ranges(GdObject)$groupid, names(GdObject@stackRanges))]
    }
    return(GdObject)
})


## Consolidate ---------------------------------------------------------------
## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------
## AlignmentTracks can be subset by using the information in the stackRanges slot, but for the actual reads we need to make sure that
## we keep all the bits that belong to a given group. We still want to record the requested ranges in the internal '.__plottingRange'
## display parameter.

#' @describeIn AlignmentsTrack-class Subset a `AlignmentsTrack` by coordinates
#' and sort if necessary.
#' @export
setMethod("subset", signature(x = "AlignmentsTrack"), function(x, from = NULL, to = NULL, stacks = FALSE, use.defaults = TRUE, ...) {
    ## Subset to a single chromosome first
    lx <- length(x)
    csel <- seqnames(x) != chromosome(x)
    if (any(csel)) {
        x <- x[!csel]
    }
    if (length(x)) {
        ## Nothing to do if everything is within the range, otherwise we subset the stackRanges and keep all group items
        ranges <- if (use.defaults) {
            .defaultRange(x, from = from, to = to)
        } else {
            c(
                from = ifelse(is.null(from), min(start(x@stackRanges)) - 1, from),
                to = ifelse(is.null(to), max(end(x@stackRanges)) + 1, to)
            )
        }
        displayPars(x) <- list(".__plottingRange" = ranges)
        sr <- subsetByOverlaps(x@stackRanges, GRanges(seqnames = chromosome(x)[1], ranges = IRanges(start = ranges["from"], end = ranges["to"])))
        if (length(sr) < length(x@stackRanges)) {
            x@stackRanges <- sr
            x@range <- x@range[x@range$groupid %in% names(sr)]
            if (stacks) {
                x <- setStacks(x)
            }
        }
    }
    return(x)
})

## ReferenceAlignmentsTracks need to stream the data from file and then pass the results on to the next method

#' @describeIn AlignmentsTrack-class Subset a `ReferenceAlignmentsTrack` by coordinates
#' and sort if necessary.
#' @export
setMethod("subset", signature(x = "ReferenceAlignmentsTrack"), function(x, from, to, chromosome, ...) {
    ## We only need to reach out into the referenced file once if the range is already contained in the object
    if (missing(from) || is.null(from) || missing(to) || is.null(to)) {
        stop("Need both start and end location to subset a ReferenceAnnotationTrack")
    }
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- Gviz::chromosome(x)
    } else {
        chromosome <- .chrName(chromosome)
    }
    subRegion <- GRanges(seqnames = chromosome[1], ranges = IRanges(start = from, end = to))
    oldRange <- .dpOrDefault(x, ".__plottingRange")
    isIn <- length(ranges(x)) != 0 && !is.null(oldRange) && ranges(subRegion) %within% IRanges(oldRange["from"], oldRange["to"])
    if (!isIn) {
        cMap <- .resolveColMapping(x@stream(x@reference, subRegion), x@args, x@mapping)
        seqs <- cMap$data$seq
        cMap$data$seq <- NULL
        range <- .computeAlignments(.buildRange(cMap$data, args = cMap$args, defaults = x@defaults, trackType = "AnnotationTrack"), drop.D.ranges = .dpOrDefault(x, "showIndels", FALSE))
        ranges(x) <- range$range
        x@stackRanges <- range$stackRanges
        x@stacks <- range$stacks
        x@sequences <- seqs
    } else {
        x@sequences <- subseq(x@sequences, start = from - (oldRange["from"] - 1), width = to - from + 1)
    }
    chromosome(x) <- chromosome[1]
    return(callNextMethod(x = x, from = from, to = to, drop = FALSE, ...))
})


## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------

setMethod("drawGrid", signature(GdObject = "AlignmentsTrack"), function(GdObject, from, to) {
    if (.dpOrDefault(GdObject, "grid", FALSE)) {
        yvals <- values(GdObject)
        ylim <- .dpOrDefault(GdObject, "ylim", if (!is.null(yvals) && length(yvals)) {
            range(yvals, na.rm = TRUE, finite = TRUE)
        } else {
            c(-1, 1)
        })
        if (diff(ylim) == 0) {
            ylim <- ylim + c(-1, 1)
        }
        yscale <- c(if (is.null(.dpOrDefault(GdObject, "transformation"))) 0 else min(ylim), max(ylim) + diff(range(ylim)) * 0.05)
        covHeight <- .dpOrDefault(GdObject, ".__coverageHeight", 0)
        pushViewport(viewport(y = 1 - covHeight["npc"], height = covHeight["npc"], just = c(0.5, 0), yscale = yscale, clip = TRUE))
        panel.grid(
            v = 0, h = .dpOrDefault(GdObject, "h", -1),
            col = .dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty = .dpOrDefault(GdObject, "lty.grid", 1),
            lwd = .dpOrDefault(GdObject, "lwd.grid", 1)
        )
        popViewport(1)
    }
})

## DrawAxis ------------------------------------------------------------------

#' @describeIn AlignmentsTrack-class add a y-axis to the title panel of a track.
#' @export
setMethod("drawAxis", signature(GdObject = "AlignmentsTrack"), function(GdObject, ...) {
    type <- match.arg(.dpOrDefault(GdObject, "type", .ALIGNMENT_TYPES), .ALIGNMENT_TYPES, several.ok = TRUE)
    if ("coverage" %in% type) {
        yvals <- values(GdObject)
        ylim <- .dpOrDefault(GdObject, "ylim", if (!is.null(yvals) && length(yvals)) {
            range(yvals, na.rm = TRUE, finite = TRUE)
        } else {
            c(-1, 1)
        })
        if (diff(ylim) == 0) {
            ylim <- ylim + c(-1, 1)
        }
        hSpaceAvail <- vpLocation()$isize["width"] / 6
        yscale <- c(if (is.null(.dpOrDefault(GdObject, "transformation"))) 0 else min(ylim), max(ylim) + diff(range(ylim)) * 0.05)
        col <- .dpOrDefault(GdObject, "col.axis", "white")
        acex <- .dpOrDefault(GdObject, "cex.axis")
        acol <- .dpOrDefault(GdObject, "col.axis", "white")
        at <- pretty(yscale)
        at <- at[at >= sort(ylim)[1] & at <= sort(ylim)[2]]
        covHeight <- .dpOrDefault(GdObject, ".__coverageHeight", c(npc = 0, points = 0))
        covSpace <- .dpOrDefault(GdObject, ".__coverageSpace", 0)
        pushViewport(viewport(y = 1 - (covHeight["npc"] + covSpace), height = covHeight["npc"], just = c(0.5, 0)))
        if (is.null(acex)) {
            vSpaceNeeded <- max(as.numeric(convertWidth(stringHeight(at), "inches"))) * length(at) * 1.5
            hSpaceNeeded <- max(as.numeric(convertWidth(stringWidth(at), "inches")))
            vSpaceAvail <- abs(diff(range(at))) / abs(diff(yscale)) * vpLocation()$isize["height"]
            acex <- max(0.6, min(vSpaceAvail / vSpaceNeeded, hSpaceAvail / hSpaceNeeded))
        }
        vpTitleAxis <- viewport(x = 0.95, width = 0.2, yscale = yscale, just = 0)
        pushViewport(vpTitleAxis)
        suppressWarnings(grid.yaxis(gp = gpar(col = acol, cex = acex), at = at))
        grid.lines(x = c(0, 0), y = ylim, gp = gpar(col = acol), default.units = "native")
        popViewport(2)
    } else {
        covHeight <- c(npc = 0, points = 0)
        covSpace <- 0
    }
    if ("sashimi" %in% type) {
        sash <- .dpOrDefault(GdObject, ".__sashimi", list(x = numeric(), y = numeric(), id = integer(), score = numeric()))
        yscale <- if (length(sash$y)) c(-(max(sash$y) + diff(range(sash$y)) * 0.05), 0) else c(-1, 0)
        ylim <- if (length(sash$y)) c(-max(sash$y), yscale[1] + max(sash$y)) else c(-1, 0)
        hSpaceAvail <- vpLocation()$isize["width"] / 6
        col <- .dpOrDefault(GdObject, "col.axis", "white")
        acex <- .dpOrDefault(GdObject, "cex.axis")
        acol <- .dpOrDefault(GdObject, "col.axis", "white")
        labs <- if (length(sash$score)) pretty(c(1, sash$score)) else pretty(c(1, .dpOrDefault(GdObject, ".__sashimiScore", 10)))
        at <- seq(ylim[1], ylim[2], length.out = length(labs))
        sashHeight <- .dpOrDefault(GdObject, ".__sashimiHeight", c(npc = 0, points = 0))
        sashSpace <- .dpOrDefault(GdObject, ".__sashimiSpace", 0)
        pushViewport(viewport(
            y = 1 - (sashHeight["npc"] + sashSpace + covHeight["npc"] + covSpace),
            height = sashHeight["npc"], just = c(0.5, 0)
        ))
        if (is.null(acex)) {
            vSpaceNeeded <- max(as.numeric(convertWidth(stringHeight(labs), "inches"))) * length(at) * 1.5
            hSpaceNeeded <- max(as.numeric(convertWidth(stringWidth(labs), "inches")))
            vSpaceAvail <- abs(diff(range(at))) / abs(diff(yscale)) * vpLocation()$isize["height"]
            acex <- max(0.6, min(vSpaceAvail / vSpaceNeeded, hSpaceAvail / hSpaceNeeded))
        }
        vpTitleAxis <- viewport(x = 0.75, width = 0.2, yscale = yscale, just = 0)
        pushViewport(vpTitleAxis)
        suppressWarnings(grid.yaxis(gp = gpar(col = acol, cex = acex), at = at, label = labs))
        grid.polygon(x = c(0, 0, 1), y = c(ylim[1], ylim[2], ylim[2]), default.units = "native", gp = gpar(col = acol, fill = acol))
        popViewport(2)
    }
    if (.dpOrDefault(GdObject, ".__isCropped", FALSE)) {
        vspacing <- as.numeric(convertHeight(unit(2, "points"), "npc"))
        vsize <- as.numeric(convertHeight(unit(4, "points"), "npc"))
        pushViewport(viewport(height = vsize, y = vspacing, just = c(0.5, 0)))
        .moreInd(direction = "down", lwd = 2)
        popViewport(1)
    }
})

## DrawGD --------------------------------------------------------------------

#' @describeIn AlignmentsTrack-class plot the object to a graphics device.
#' The return value of this method is the input object, potentially updated
#' during the plotting operation. Internally, there are two modes in which the
#' method can be called. Either in 'prepare' mode, in which case no plotting is
#' done but the object is preprocessed based on the available space, or in
#' 'plotting' mode, in which case the actual graphical output is created.
#' Since subsetting of the object can be potentially costly, this can be
#' switched off in case subsetting has already been performed before or
#' is not necessary.
#'
#' @importFrom GenomicAlignments cigarRangesAlongReferenceSpace
#' @export
setMethod("drawGD", signature("AlignmentsTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, subset = TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    imageMap(GdObject) <- NULL
    if (!length(GdObject)) {
        return(invisible(GdObject))
    }
    rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
    revs <- !.dpOrDefault(GdObject, "reverseStacking", FALSE)
    vSpacing <- as.numeric(convertHeight(unit(3, "points"), "npc"))
    type <- match.arg(.dpOrDefault(GdObject, "type", .ALIGNMENT_TYPES), .ALIGNMENT_TYPES, several.ok = TRUE)
    ylim <- .dpOrDefault(GdObject, "ylim")
    ## In prepare mode we need to make sure that the stacking information is updated from the optional display parameter (by calling
    ## the StackedTrack drawGD method), and compute the coverage vector which will be needed for the axis
    if (prepare) {
        if ((is.logical(debug) && debug) || "prepare" %in% debug) {
            browser()
        }
        GdObject <- callNextMethod(GdObject, ...)
        ## The mismatched bases need to be extracted from the read sequences and the reference sequence
        if (.dpOrDefault(GdObject, "showMismatches", TRUE) && !is.null(GdObject@referenceSequence)) {
            mm <- .findMismatches(GdObject)
            if (nrow(mm)) {
                displayPars(GdObject) <- list(".__mismatches" = mm)
            }
        }
        ## The coverage calculation and the height of the coverage section
        if ("coverage" %in% type) {
            covSpace <- as.numeric(convertHeight(unit(5, "points"), "npc"))
            if ("pileup" %in% type) {
                covHeight <- .dpOrDefault(GdObject, "coverageHeight", 0.1)
                if (covHeight > 0 && covHeight < 1) {
                    covHeight <- as.numeric(convertHeight(unit(covHeight, "npc"), "points"))
                }
                covHeight <- max(.dpOrDefault(GdObject, "minCoverageHeight", 50), covHeight)
                covHeight <- c(points = covHeight, npc = as.numeric(convertHeight(unit(covHeight, "points"), "npc")))
            } else if ("sashimi" %in% type) {
                covHeight <- c(npc = 0.5 - (covSpace * 2), points = 1)
            } else {
                covHeight <- c(npc = 1 - (covSpace * 2), points = 1)
            }
            coverage <- coverage(range(GdObject), width = maxBase)
            ## apply data transformation if one is set up
            trans <- .dpOrDefault(GdObject, "transformation")
            if (is.list(trans)) {
                trans <- trans[[1]]
            }
            if (!is.null(trans)) {
                if (!is.function(trans) || length(formals(trans)) != 1L) {
                    stop("Display parameter 'transformation' must be a function with a single argument")
                }
                test <- trans(coverage)
                if (!is(test, "Rle") || length(test) != length(coverage)) {
                    stop(
                        "The function in display parameter 'transformation' results in invalid output.\n",
                        "It has to return a numeric matrix with the same dimensions as the input data."
                    )
                }
                coverage <- test
            }
            displayPars(GdObject) <- list(
                ".__coverage" = coverage,
                ".__coverageHeight" = covHeight,
                ".__coverageSpace" = covSpace
            )
        } else {
            covHeight <- c(npc = 0, points = 0)
            covSpace <- 0
        }
        ## The sashimi calculation and the height of the sashimi section
        if ("sashimi" %in% type) {
            sashSpace <- as.numeric(convertHeight(unit(5, "points"), "npc"))
            if ("pileup" %in% type) {
                sashHeight <- .dpOrDefault(GdObject, "sashimiHeight", 0.1)
                if (sashHeight > 0 && sashHeight < 1) {
                    sashHeight <- as.numeric(convertHeight(unit(sashHeight, "npc"), "points"))
                }
                sashHeight <- max(.dpOrDefault(GdObject, "minSashimiHeight", 50), sashHeight)
                sashHeight <- c(points = sashHeight, npc = as.numeric(convertHeight(unit(sashHeight, "points"), "npc")))
            } else if ("coverage" %in% type) {
                sashHeight <- c(npc = 0.5 - (sashSpace * 2), points = 1)
            } else {
                sashHeight <- c(npc = 1 - (sashSpace * 2), points = 1)
            }
            sashScore <- .dpOrDefault(GdObject, "sashimiScore", 1L)
            sashLwdMax <- .dpOrDefault(GdObject, "lwd.sashimiMax", 10)
            sashStrand <- .dpOrDefault(GdObject, "sashimiStrand", "*")
            sashFilter <- .dpOrDefault(GdObject, "sashimiFilter", NULL)
            sashFilterTolerance <- .dpOrDefault(GdObject, "sashimiFilterTolerance", 0L)
            sashNumbers <- .dpOrDefault(GdObject, "sashimiNumbers", FALSE)
            sash <- .dpOrDefault(GdObject, "sashimiJunctions", NULL)
            if (is.null(sash)) {
                sash <- .create.summarizedJunctions.for.sashimi.junctions(ranges(GdObject))
            } else {
                if (!is(sash, "GRanges")) {
                    stop("\"sashimiJunctions\" object must be of \"GRanges\" class!")
                }
                sashMcolName <- if (sashStrand == "+") "plus_score" else if (sashStrand == "-") "minus_score" else "score"
                if (sum(colnames(mcols(sash)) == sashMcolName) != 1) {
                    stop(sprintf("\"mcols\" of \"sashimiJunctions\" object must contain column named \"%s\",\n which matches the specified (%s) \"sashimiStrand\"!", sashMcolName, sashStrand))
                }
            }
            sashTransform <- .dpOrDefault(GdObject, c("sashimiTransformation", "transformation"))
            sash <- .convert.summarizedJunctions.to.sashimi.junctions(
                juns = sash,
                score = sashScore,
                lwd.max = sashLwdMax,
                strand = sashStrand,
                filter = sashFilter,
                filterTolerance = sashFilterTolerance,
                trans = sashTransform
            )
            displayPars(GdObject) <- list(
                ".__sashimi" = sash,
                ".__sashimiHeight" = sashHeight,
                ".__sashimiSpace" = sashSpace,
                ".__sashimiNumbers" = sashNumbers
            )
        } else {
            sashHeight <- c(npc = 0, points = 0)
            sashSpace <- 0
        }
        if ("pileup" %in% type) {
            ## If there are more bins than we can plot we reduce the number until they fit
            pushViewport(viewport(height = 1 - (covHeight["npc"] + covSpace + sashHeight["npc"] + sashSpace) - vSpacing * 2, y = vSpacing, just = c(0.5, 0)))
            bins <- ranges(GdObject)$stack
            curVp <- vpLocation()
            mih <- min(curVp$size["height"], .dpOrDefault(GdObject, c("min.height", "max.height"), 3))
            mah <- min(curVp$size["height"], max(mih, .dpOrDefault(GdObject, c("max.height", "min.height"), 8)))
            bins <- stacks(GdObject)
            if (curVp$size["height"] / max(bins) < mih) {
                maxStack <- curVp$size["height"] %/% mih
                sel <- if (revs) bins <= maxStack else bins > max(bins) - maxStack
                ranges(GdObject) <- ranges(GdObject)[sel]
                displayPars(GdObject) <- list(".__isCropped" = TRUE)
                bins <- stacks(GdObject)
            }
            yrange <- range(bins) + c(-0.5, 0.5)
            add <- max(0, (curVp$size["height"] %/% mah) - max(bins))
            yrange <- if (revs) c(yrange[1], yrange[2] + add) else c(yrange[1] - add, yrange[2])
            displayPars(GdObject) <- list(".__yrange" = yrange)
            popViewport(1)
        }
        return(invisible(GdObject))
    }
    if ((is.logical(debug) && debug) || "draw" %in% debug) {
        browser()
    }
    mm <- .dpOrDefault(GdObject, ".__mismatches")
    ## The coverage plot first
    xscale <- if (!rev) c(minBase, maxBase) else c(maxBase, minBase)
    readInfo <- ranges(GdObject)
    covHeight <- .dpOrDefault(GdObject, ".__coverageHeight", c(npc = 0, points = 0))
    covSpace <- .dpOrDefault(GdObject, ".__coverageSpace", 0)
    if ("coverage" %in% type) {
        cov <- .dpOrDefault(GdObject, ".__coverage", Rle(lengths = maxBase, values = as.integer(0)))
        yscale <- c(if (is.null(.dpOrDefault(GdObject, "transformation"))) 0 else min(cov), max(cov) + diff(range(cov)) * 0.05)
        if (!is.null(ylim)) {
            yscale <- ylim
        }
        vp <- viewport(
            height = covHeight["npc"], y = 1 - (covHeight["npc"] + covSpace), just = c(0.5, 0), xscale = xscale,
            yscale = yscale, clip = TRUE
        )
        pushViewport(vp)
        res <- .pxResolution(coord = "x")
        gp <- gpar(
            col = .dpOrDefault(GdObject, c("col.coverage", "col"), .DEFAULT_SHADED_COL),
            fill = .dpOrDefault(GdObject, c("fill.coverage", "fill"), "#BABABA"),
            lwd = .dpOrDefault(GdObject, c("lwd.coverage", "lwd"), 1),
            lty = .dpOrDefault(GdObject, c("lty.coverage", "lty"), 1),
            alpha = .alpha(GdObject)
        )
        ## We can compact this if the resolution is not sufficient to speed up the drawing.
        mminBase <- max(1, minBase)
        if (res > 2) {
            brks <- ceiling((maxBase - mminBase) / res)
            x <- seq(mminBase, maxBase, len = brks)
            y <- tapply(as.integer(cov[mminBase:maxBase]), cut(mminBase:maxBase, breaks = brks), mean)
        } else {
            x <- mminBase:maxBase
            y <- as.integer(cov[x])
        }
        grid.polygon(c(minBase - max(1, res), 0, x, maxBase + max(1, res)), c(0, 0, y, 0), default.units = "native", gp = gp)
        grid.lines(y = c(0, 0), gp = gpar(col = gp$col, alpha = gp$alpha))
        if (!is.null(mm)) {
            fcol <- .dpOrDefault(GdObject@referenceSequence, "fontcolor", getBioColor("DNA_BASES_N"))
            vpos <- tapply(as.character(mm$base[mm$base != "."]), mm$position[mm$base != "."], table, simplify = FALSE)
            x <- rep(as.integer(names(vpos)), listLen(vpos))
            y <- unlist(lapply(vpos, cumsum), use.names = FALSE)
            col <- fcol[unlist(lapply(vpos, names), use.names = FALSE)]
            grid.rect(
                x = x, y = y, height = unlist(vpos, use.names = FALSE), width = 1, default.units = "native", just = c(0, 1),
                gp = gpar(col = "transparent", fill = col)
            )
        }
        popViewport(1)
        twoPx <- 2 * as.numeric(convertHeight(unit(1, "points"), "npc"))
        vp <- viewport(height = twoPx, y = 1 - (covHeight["npc"] + covSpace + twoPx), just = c(0.5, 0))
        pushViewport(vp)
        grid.rect(gp = gpar(
            fill = .dpOrDefault(GdObject, "background.title"), col = "transparent",
            alpha = .dpOrDefault(GdObject, "alpha.title")
        ), width = 2)
        popViewport(1)
    }
    ## The sashimi plot as second
    xscale <- if (!rev) c(minBase, maxBase) else c(maxBase, minBase)
    sashHeight <- .dpOrDefault(GdObject, ".__sashimiHeight", c(npc = 0, points = 0))
    sashSpace <- .dpOrDefault(GdObject, ".__sashimiSpace", 0)
    if ("sashimi" %in% type) {
        sash <- .dpOrDefault(GdObject, ".__sashimi", list(x = numeric(), y = numeric(), id = integer(), score = numeric(), scaled = numeric()))
        sashNumbers <- .dpOrDefault(GdObject, ".__sashimiNumbers", FALSE)
        yscale <- if (length(sash$y)) c(-(max(sash$y) + diff(range(sash$y)) * 0.15), 0) else c(-1, 0) # changed from 0.05 to 0.15 to make sure that numbers fit in the viewport
        vp <- viewport(
            height = sashHeight["npc"], y = 1 - (covHeight["npc"] + covSpace + sashHeight["npc"] + sashSpace), just = c(0.5, 0),
            xscale = xscale, yscale = yscale, clip = TRUE
        )
        pushViewport(vp)
        gp <- gpar(
            col = .dpOrDefault(GdObject, c("col.sashimi", "col"), .DEFAULT_SHADED_COL),
            fill = .dpOrDefault(GdObject, c("fill.sashimi", "fill"), "#FFFFFF"),
            lwd = .dpOrDefault(GdObject, c("lwd.sashimi", "lwd"), 1),
            lty = .dpOrDefault(GdObject, c("lty.sashimi", "lty"), 1),
            alpha = .alpha(GdObject)
        )
        if (length(sash$x)) {
            grid.xspline(sash$x, -sash$y,
                id = sash$id, shape = -1, open = TRUE,
                default.units = "native", gp = gpar(col = gp$col, lwd = sash$scaled)
            )
            ## print the number of reads together with the connecting lines (currently no scaling/resolution)
            if (sashNumbers) {
                grid.rect(sash$x[c(FALSE, TRUE, FALSE)], -sash$y[c(FALSE, TRUE, FALSE)],
                    width = convertUnit(stringWidth(sash$score) * 1.5, "inches"),
                    height = convertUnit(stringHeight(sash$score) * 1.5, "inches"),
                    default.units = "native", gp = gpar(col = gp$col, fill = gp$fill)
                )
                grid.text(
                    label = sash$score, sash$x[c(FALSE, TRUE, FALSE)], -sash$y[c(FALSE, TRUE, FALSE)],
                    default.units = "native", gp = gpar(col = gp$col)
                )
            }
        }
        popViewport(1)
    }
    ## Now the pileups
    if ("pileup" %in% type) {
        cex <- max(0.3, .dpOrDefault(GdObject, c("cex.mismatch", "cex"), 0.7))
        pushViewport(viewport(
            height = 1 - (covHeight["npc"] + covSpace + sashHeight["npc"] + sashSpace) - vSpacing * 4, y = vSpacing, just = c(0.5, 0),
            gp = .fontGp(GdObject, "mismatch", cex = cex)
        ))
        bins <- stacks(GdObject)
        stacks <- max(bins)
        yscale <- .dpOrDefault(GdObject, ".__yrange")
        if (revs) {
            yscale <- rev(yscale)
        }
        pushViewport(dataViewport(xscale = xscale, extension = 0, yscale = yscale, clip = TRUE))
        ## Figuring out resolution and the necessary plotting details
        res <- .pxResolution(coord = "x")
        ylim <- c(0, 1)
        h <- diff(ylim)
        middle <- mean(ylim)
        sh <- max(0, min(h, .dpOrDefault(GdObject, "stackHeight", 0.75))) / 2
        boxOnly <- res > 10
        if (boxOnly) {
            x <- c(start(readInfo), rep(end(readInfo) + 1, 2), start(readInfo))
            y <- c(rep(readInfo$stack + sh, 2), rep(readInfo$stack - sh, 2))
            id <- rep(readInfo$uid, 4)
        } else {
            ## We first precompute the coordinates for all the arrow polygons
            uid2strand <- setNames(as.character(readInfo$readStrand), as.character(readInfo$uid))
            arrowMap <- unlist(setNames(lapply(split(as.character(readInfo$uid), readInfo$entityId), function(x) {
                str <- uid2strand[x][1]
                setNames(str, if (str == "+") tail(x, 1) else head(x, 1))
            }), NULL))
            readInfo$arrow <- as.character(NA)
            readInfo$arrow[match(names(arrowMap), readInfo$uid)] <- arrowMap
            ## The parts that don't need arrow heads
            sel <- is.na(readInfo$arrow) | readInfo$arrow == "*"
            x <- c(start(readInfo)[sel], rep(end(readInfo)[sel] + 1, 2), start(readInfo)[sel])
            y <- c(rep(readInfo$stack[sel] + sh, 2), rep(readInfo$stack[sel] - sh, 2))
            id <- rep(readInfo$uid[sel], 4)
            ## The arrow heads facing right
            w <- Gviz:::.pxResolution(coord = "x", 5)
            sel <- readInfo$arrow == "+"
            ah <- pmax(start(readInfo)[sel], end(readInfo)[sel] + 1 - w)
            x <- c(x, start(readInfo)[sel], ah, end(readInfo)[sel] + 1, ah, start(readInfo)[sel])
            y <- c(y, rep(readInfo$stack[sel] + sh, 2), readInfo$stack[sel], rep(readInfo$stack[sel] - sh, 2))
            id <- c(id, rep(readInfo$uid[sel], 5))
            ## The arrow heads facing left
            sel <- readInfo$arrow == "-"
            ah <- pmin(end(readInfo)[sel], start(readInfo)[sel] + w)
            x <- c(x, start(readInfo)[sel], ah, rep(end(readInfo)[sel] + 1, 2), ah)
            y <- c(y, readInfo$stack[sel], rep(readInfo$stack[sel] + sh, 2), rep(readInfo$stack[sel] - sh, 2))
            id <- c(id, rep(readInfo$uid[sel], 5))
        }
        nn <- length(unique(readInfo$uid))
        gps <- data.frame(
            col = rep(.dpOrDefault(GdObject, c("col.reads", "col"), .DEFAULT_SHADED_COL), nn),
            fill = rep(.dpOrDefault(GdObject, c("fill.reads", "fill"), .DEFAULT_BRIGHT_SHADED_COL), nn),
            lwd = rep(.dpOrDefault(GdObject, c("lwd.reads", "lwd"), 1), nn),
            lty = rep(.dpOrDefault(GdObject, c("lty.reads", "lty"), 1), nn),
            alpha = rep(.alpha(GdObject, "reads"), nn), stringsAsFactors = FALSE
        )
        ## Finally we draw reads (we need to draw in two steps because of the different alpha levels, reads vs. mismatches)
        grid.polygon(x = x, y = y, id = id, default.units = "native", gp = gpar(col = gps$col, fill = gps$fill, lwd = gps$lwd, lty = gps$lty, alpha = unique(gps$alpha))) # fix for Sys.setenv(`_R_CHECK_LENGTH_1_CONDITION_`="true"); Sys.setenv(`_R_CHECK_LENGTH_1_LOGIC2_`="true")
        ## Now the coordinates for the connecting lines
        lineCoords <- NULL
        if (anyDuplicated(readInfo$entityId) != 0) {
            stTmp <- split(readInfo, readInfo$entityId)
            mateRanges <- unlist(range(stTmp))
            mateGaps <- gaps(GRanges(ranges = ranges(readInfo), seqnames = readInfo$entityId))
            rmap <- mateRanges[as.character(seqnames(mateGaps))]
            mateGaps <- mateGaps[start(rmap) <= start(mateGaps) & end(rmap) >= end(mateGaps)]
            gy <- readInfo$stack[match(as.character(seqnames(mateGaps)), readInfo$entityId)]
            lineCoords <- data.frame(
                x1 = start(mateGaps), y1 = gy, x2 = end(mateGaps) + 1, y2 = gy,
                col = .dpOrDefault(GdObject, c("col.gap", "col"), .DEFAULT_SHADED_COL),
                lwd = .dpOrDefault(GdObject, c("lwd.gap", "lwd"), 1),
                lty = .dpOrDefault(GdObject, c("lty.gap", "lty"), 1),
                alpha = .alpha(GdObject, "gap"), stringsAsFactors = FALSE
            )
            lineCoords <- lineCoords[!duplicated(lineCoords), ]
        } else {
            mateRanges <- setNames(readInfo, readInfo$entityId)
        }
        if (any(readInfo$status != "unmated") && anyDuplicated(readInfo$groupid) != 0) {
            pairGaps <- gaps(GRanges(ranges = ranges(mateRanges), seqnames = readInfo$groupid[match(names(mateRanges), readInfo$entityId)]))
            rmap <- GdObject@stackRanges[as.character(seqnames(pairGaps))]
            pairGaps <- pairGaps[start(rmap) <= start(pairGaps) & end(rmap) >= end(pairGaps)]
            gy <- readInfo$stack[match(as.character(seqnames(pairGaps)), readInfo$groupid)]
            if (length(pairGaps)) {
                pairsCoords <- data.frame(
                    x1 = start(pairGaps) - 1, y1 = gy, x2 = end(pairGaps) + 1, y2 = gy,
                    col = .dpOrDefault(GdObject, c("col.mates", "col"), .DEFAULT_BRIGHT_SHADED_COL),
                    lwd = .dpOrDefault(GdObject, c("lwd.mates", "lwd"), 1),
                    lty = .dpOrDefault(GdObject, c("lty.mates", "lty"), 1),
                    alpha = .alpha(GdObject, "mates"), stringsAsFactors = FALSE
                )
            } else {
                pairsCoords <- NULL
            }
            lineCoords <- rbind(lineCoords, pairsCoords[!duplicated(pairsCoords), ])
        }
        ## Consider the indels if needed
        ## - deletion as lines
        ## - insertions as vertical bars
        showIndels <- .dpOrDefault(GdObject, "showIndels", FALSE)
        delCoords <- NULL
        insCoords <- NULL
        if (showIndels) {
            cigarTmp <- DataFrame(cigar = readInfo$cigar, start = start(readInfo), entityId = readInfo$entityId, groupId = readInfo$groupid)
            cigarTmp <- cigarTmp[order(cigarTmp$entityId, cigarTmp$start), ]
            cigarTmp <- cigarTmp[!duplicated(cigarTmp$entityId), ]
            delGaps <- unlist(cigarRangesAlongReferenceSpace(cigarTmp$cigar, pos = cigarTmp$start, ops = "D", f = as.factor(cigarTmp$entityId)))
            gy <- readInfo$stack[match(names(delGaps), readInfo$entityId)]
            if (length(delGaps)) {
                delCoords <- data.frame(
                    x1 = start(delGaps) + 1, y1 = gy, x2 = end(delGaps) + 1, y2 = gy,
                    col = .dpOrDefault(GdObject, c("col.deletion", "col"), .DEFAULT_BRIGHT_SHADED_COL),
                    lwd = .dpOrDefault(GdObject, c("lwd.deletion", "lwd"), 1),
                    lty = .dpOrDefault(GdObject, c("lty.deletion", "lty"), 1),
                    alpha = .alpha(GdObject, "deletions"), stringsAsFactors = FALSE
                )
                lineCoords <- rbind(delCoords, lineCoords)
                lineCoords <- lineCoords[!duplicated(lineCoords[, c("x1", "y1", "x2", "y2")]), ]
            }
            insGaps <- unlist(cigarRangesAlongReferenceSpace(cigarTmp$cigar, pos = cigarTmp$start, ops = "I", f = as.factor(cigarTmp$entityId)))
            gy <- readInfo$stack[match(names(insGaps), readInfo$entityId)]
            if (length(insGaps)) {
                ## should both x coordinates be equal to start
                insCoords <- data.frame(
                    x1 = start(insGaps), y1 = gy - sh, x2 = start(insGaps), y2 = gy + sh,
                    col = .dpOrDefault(GdObject, c("col.insertion", "col"), .DEFAULT_BRIGHT_SHADED_COL),
                    lwd = .dpOrDefault(GdObject, c("lwd.insertion", "lwd"), 1),
                    lty = .dpOrDefault(GdObject, c("lty.insertion", "lty"), 1),
                    alpha = .alpha(GdObject, "insertions"), stringsAsFactors = FALSE
                )
            }
        }
        ## The mismatch information on the reads if needed
        mmLetters <- NULL
        if (!is.null(mm)) {
            fcol <- .dpOrDefault(GdObject@referenceSequence, "fontcolor", getBioColor("DNA_BASES_N"))
            ccol <- ifelse(rgb2hsv(col2rgb(fcol))["s", ] < 0.5, "black", "white")
            vpl <- vpLocation()
            lwidth <- max(as.numeric(convertUnit(stringWidth(DNA_ALPHABET), "inches"))) * 0.9337632
            lheight <- max(as.numeric(convertUnit(stringHeight(DNA_ALPHABET), "inches")))
            perLetterW <- vpl$isize["width"] / (maxBase - minBase + 1)
            perLetterH <- vpl$isize["height"] / abs(diff(current.viewport()$yscale))
            res <- .pxResolution(coord = "x")
            mw <- res * .dpOrDefault(GdObject, "min.width", 1)
            mwy <- max(1, mw)
            if (nrow(mm)) {
                x <- c(mm$position + rep(c(0, mwy, mwy, 0), each = nrow(mm)))
                y <- c(rep(mm$stack - sh, 2), rep(mm$stack + sh, 2))
                id <- c(rep(seq(max(id, na.rm = TRUE) + 1, len = nrow(mm)), 4))
                gps <- data.frame(
                    col = rep(if (!(lwidth < perLetterW && lheight < perLetterH)) {
                        "transparent"
                    } else {
                        .dpOrDefault(GdObject, "col.mismatch", .DEFAULT_SHADED_COL)
                    }, nrow(mm)),
                    fill = rep(fcol[as.character(mm$base)]),
                    lwd = rep(.dpOrDefault(GdObject, "lwd.mismatch", 1), nrow(mm)),
                    lty = rep(.dpOrDefault(GdObject, "lty.mismatch", 1), nrow(mm)),
                    alpha = rep(.dpOrDefault(GdObject, "alpha.mismatch", 1), nrow(mm)),
                    stringsAsFactors = FALSE
                )
                ## Finally we draw mm (we need to draw in two steps because of the different alpha levels, reads vs. mismatches)
                grid.polygon(x = x, y = y, id = id, default.units = "native", gp = gpar(col = gps$col, fill = gps$fill, lwd = gps$lwd, lty = gps$lty, alpha = unique(gps$alpha))) # fix for Sys.setenv(`_R_CHECK_LENGTH_1_CONDITION_`="true"); Sys.setenv(`_R_CHECK_LENGTH_1_LOGIC2_`="true")
                if (!.dpOrDefault(GdObject, "noLetters", FALSE) && lwidth < perLetterW && lheight < perLetterH) {
                    mmLetters <- data.frame(x = mm$position + 0.5, y = mm$stack, label = mm$base, col = ccol[mm$base], stringsAsFactors = FALSE)
                }
            }
        }
        if (!is.null(lineCoords)) {
            grid.segments(lineCoords$x1, lineCoords$y1, lineCoords$x2, lineCoords$y2,
                gp = gpar(
                    col = lineCoords$col, alpha = unique(lineCoords$alpha),
                    lwd = lineCoords$lwd, lty = lineCoords$lty
                ),
                default.units = "native"
            )
        }
        if (!is.null(insCoords)) {
            grid.segments(insCoords$x1, insCoords$y1, insCoords$x2, insCoords$y2,
                gp = gpar(
                    col = insCoords$col, alpha = unique(insCoords$alpha),
                    lwd = insCoords$lwd, lty = insCoords$lty
                ),
                default.units = "native"
            )
        }
        if (!is.null(mmLetters)) {
            grid.text(
                x = mmLetters$x, y = mmLetters$y, label = mmLetters$label,
                gp = gpar(col = mmLetters$col), default.units = "native"
            )
        }
        popViewport(2)
    }
    ## Eventually we set up the image map
    ## imageMap(GdObject) <- im
    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------
## Show ----------------------------------------------------------------------

#' @describeIn AlignmentsTrack-class Show method.
#' @export
setMethod(
    "show", signature(object = "AlignmentsTrack"),
    function(object) {
        cat(sprintf(
            paste("AlignmentsTrack track '%s' \n",
                "| genome: %s\n",
                "| active chromosome: %s\n",
                "| containing %i read%s\n",
                sep = ""
            ),
            names(object), genome(object),
            gsub("^chr", "", chromosome(object)),
            length(object),
            ifelse(length(object) == 1, "", "s")
        ))
    }
)

#' @describeIn AlignmentsTrack-class Show method.
#' @export
setMethod("show", signature(object = "ReferenceAlignmentsTrack"), function(object) {
    .referenceTrackInfo(object, "ReferenceAlignmentsTrack")
})
