#' @include AnnotationTrack-class.R
NULL

## GeneRegionTrack Class -----------------------------------------------------


#' GeneRegionTrack class and methods
#'
#'
#' A class to hold gene model data for a genomic region.
#'
#'
#' A track containing all gene models in a particular region. The data are
#' usually fetched dynamially from an online data store, but it is also
#' possible to manully construct objects from local data. Connections to
#' particular online data sources should be implemented as sub-classes, and
#' \code{GeneRegionTrack} is just the commone denominator that is being used
#' for plotting later on. There are several levels of data associated to a
#' \code{GeneRegionTrack}:
#'
#' \describe{
#'
#' \item{exon level:}{identifiers are stored in the exon column of the
#' \code{\linkS4class{GRanges}} object in the \code{range} slot. Data may be
#' extracted using the \code{exon} method.}
#'
#' \item{transcript level:}{identifiers are stored in the transcript column of
#' the \code{\linkS4class{GRanges}} object. Data may be extracted using the
#' \code{transcript} method.}
#'
#' \item{gene level:}{identifiers are stored in the gene column of the
#' \code{\linkS4class{GRanges}} object, more human-readable versions in the
#' symbol column. Data may be extracted using the \code{gene} or the
#' \code{symbol} methods.}
#'
#' \item{transcript-type level:}{information is stored in the feature column of
#' the \code{\linkS4class{GRanges}} object. If a display parameter of the same
#' name is specified, the software will use its value for the coloring.}
#'
#' }
#'
#' \code{GeneRegionTrack} objects also know about coding regions and non-coding
#' regions (e.g., UTRs) in a transcript, and will indicate those by using
#' different shapes (wide boxes for all coding regions, thinner boxes for
#' non-coding regions). This is archived by setting the \code{feature} values
#' of the object for non-coding elements to one of the options that are
#' provided in the \code{thinBoxFeature} display parameters. All other elements
#' are considered to be coding elements.
#'
#' @name GeneRegionTrack-class
#'
#' @param range
#'
#' An optional meta argument to handle the different input types. If the
#' \code{range} argument is missing, all the relevant information to create the
#' object has to be provided as individual function arguments (see below).
#'
#' The different input options for \code{range} are:
#'
#' \describe{
#'
#' \item{A \code{TxDb} object:}{ all the necessary gene model information
#' including exon locations, transcript groupings and associated gene ids are
#' contained in \code{TxDb} objects, and the coercion between the two is almost
#' completely automated. If desired, the data to be fetched from the
#' \code{TxDb} object can be restricted using the constructor's
#' \code{chromosome}, \code{start} and \code{end} arguments. See below for
#' details. A direct coercion method \code{as(obj, "GeneRegionTrack")} is also
#' available. A nice added benefit of this input option is that the UTR and
#' coding region information that is part of the original \code{TxDb} object is
#' retained in the \code{GeneRegionTrack}.}
#'
#' \item{A \code{GRanges} object:}{ the genomic ranges for the
#' \code{GeneRegion} track as well as the optional additional metadata columns
#' \code{feature}, \code{transcript}, \code{gene}, \code{exon} and
#' \code{symbol} (see description of the individual function parameters below
#' for details). Calling the constructor on a \code{GRanges} object without
#' further arguments, e.g.  \code{GeneRegionTrack(range=obj)} is equivalent to
#' calling the coerce method \code{as(obj, "GeneRegionTrack")}.}
#'
#' \item{A \code{GRangesList} object:}{ this is very similar to the previous
#' case, except that the grouping information that is part of the list
#' structure is preserved in the \code{GeneRegionTrack}. I.e., all the elements
#' within one list item receive the same group id. For consistancy, there is
#' also a coercion method from \code{GRangesLists} \code{as(obj,
#' "GeneRegionTrack")}. Please note that unless the necessary information about
#' gene ids, symbols, etc. is present in the individual \code{GRanges} meta
#' data slots, the object will not be particularly useful, because all the
#' identifiers will be set to a common default value.}
#'
#' \item{An \code{\linkS4class{IRanges}} object:}{ almost identical to the
#' \code{GRanges} case, except that the chromosome and strand information as
#' well as all additional data has to be provided in the separate
#' \code{chromosome}, \code{strand}, \code{feature}, \code{transcript},
#' \code{symbol}, \code{exon} or \code{gene} arguments, because it can not be
#' directly encoded in an \code{IRanges} object. Note that only the former two
#' are mandatory (if not provided explicitely the more or less reasonable
#' default values \code{chromosome=NA} and \code{strand=*} are used, but not
#' providing information about the gene-to-transcript relationship or the
#' human-readble symbols renders a lot of the class' functionality useles.}
#'
#' \item{A \code{data.frame} object:}{ the \code{data.frame} needs to contain
#' at least the two mandatory columns \code{start} and \code{end} with the
#' range coordinates. It may also contain a \code{chromosome} and a
#' \code{strand} column with the chromosome and strand information for each
#' range. If missing, this information will be drawn from the constructor's
#' \code{chromosome} or \code{strand} arguments. In addition, the
#' \code{feature}, \code{exon}, \code{transcript}, \code{gene} and
#' \code{symbol} data can be provided as columns in the \code{data.frame}. The
#' above comments about potential default values also apply here.}
#'
#' \item{A \code{character} scalar:}{ in this case the value of the
#' \code{range} argument is considered to be a file path to an annotation file
#' on disk. A range of file types are supported by the \code{Gviz} package as
#' identified by the file extension. See the \code{importFunction}
#' documentation below for further details.}
#'
#' }
#' @template GeneRegionTrack-class_param
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{GeneRegionTrack}.
#' @section Objects from the class:
#'
#' Objects can be created using the constructor function
#' \code{GeneRegionTrack}.
#'
#' @author Florian Hahne, Steve Lianoglou
#'
#' @inherit GdObject-class seealso
#'
#' @examples
#'
#'
#' ## The empty object
#' GeneRegionTrack()
#'
#' ## Load some sample data
#' data(cyp2b10)
#'
#' ## Construct the object
#' grTrack <- GeneRegionTrack(
#'     start = 26682683, end = 26711643,
#'     rstart = cyp2b10$start, rends = cyp2b10$end, chromosome = 7, genome = "mm9",
#'     transcript = cyp2b10$transcript, gene = cyp2b10$gene, symbol = cyp2b10$symbol,
#'     feature = cyp2b10$feature, exon = cyp2b10$exon,
#'     name = "Cyp2b10", strand = cyp2b10$strand
#' )
#'
#' ## Directly from the data.frame
#' grTrack <- GeneRegionTrack(cyp2b10)
#'
#' ## From a TxDb object
#' if (require(GenomicFeatures)) {
#'     samplefile <- system.file("extdata",
#'                               "hg19_knownGene_sample.sqlite",
#'                               package = "GenomicFeatures")
#'     txdb <- loadDb(samplefile)
#'     GeneRegionTrack(txdb)
#'     GeneRegionTrack(txdb, chromosome = "chr6", start = 35000000, end = 40000000)
#' }
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(grTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#' ## Plotting
#' plotTracks(grTrack)
#'
#' ## Track names
#' names(grTrack)
#' names(grTrack) <- "foo"
#' plotTracks(grTrack)
#'
#' ## Subsetting and splitting
#' subTrack <- subset(grTrack, from = 26700000, to = 26705000)
#' length(subTrack)
#' subTrack <- grTrack[transcript(grTrack) == "ENSMUST00000144140"]
#' split(grTrack, transcript(grTrack))
#'
#' ## Accessors
#' start(grTrack)
#' end(grTrack)
#' width(grTrack)
#' position(grTrack)
#' width(subTrack) <- width(subTrack) + 100
#'
#' strand(grTrack)
#' strand(subTrack) <- "-"
#'
#' chromosome(grTrack)
#' chromosome(subTrack) <- "chrX"
#'
#' genome(grTrack)
#' genome(subTrack) <- "hg19"
#'
#' range(grTrack)
#' ranges(grTrack)
#'
#' ## Annotation
#' identifier(grTrack)
#' identifier(grTrack, "lowest")
#' identifier(subTrack) <- "bar"
#'
#' feature(grTrack)
#' feature(subTrack) <- "foo"
#'
#' exon(grTrack)
#' exon(subTrack) <- letters[1:2]
#'
#' gene(grTrack)
#' gene(subTrack) <- "bar"
#'
#' symbol(grTrack)
#' symbol(subTrack) <- "foo"
#'
#' transcript(grTrack)
#' transcript(subTrack) <- c("foo", "bar")
#' chromosome(subTrack) <- "chr7"
#' plotTracks(subTrack)
#'
#' values(grTrack)
#'
#' ## Grouping
#' group(grTrack)
#' group(subTrack) <- "Group 1"
#' transcript(subTrack)
#' plotTracks(subTrack)
#'
#' ## Collapsing transcripts
#' plotTracks(grTrack,
#'     collapseTranscripts = TRUE, showId = TRUE,
#'     extend.left = 10000, shape = "arrow"
#' )
#'
#' ## Stacking
#' stacking(grTrack)
#' stacking(grTrack) <- "dense"
#' plotTracks(grTrack)
#'
#' ## coercion
#' as(grTrack, "data.frame")
#' as(grTrack, "UCSCData")
#'
#' ## HTML image map
#' coords(grTrack)
#' tags(grTrack)
#' grTrack <- plotTracks(grTrack)$foo
#' coords(grTrack)
#' tags(grTrack)
#' @exportClass GeneRegionTrack
setClass("GeneRegionTrack",
    contains = "AnnotationTrack",
    representation = representation(start = "NumericOrNULL", end = "NumericOrNULL"),
    prototype = prototype(
        columns = c("feature", "transcript", "symbol", "gene", "exon"),
        stacking = "squish",
        stacks = 0,
        start = 0,
        end = 0,
        name = "GeneRegionTrack",
        dp = DisplayPars(
            arrowHeadWidth = 10,
            arrowHeadMaxWidth = 20,
            col = NULL,
            collapseTranscripts = FALSE,
            exonAnnotation = NULL,
            fill = "orange",
            min.distance = 0,
            shape = c("smallArrow", "box"),
            showExonId = NULL,
            thinBoxFeature = .THIN_BOX_FEATURES,
            transcriptAnnotation = NULL
        )
    )
)

## Initialize ----------------------------------------------------------------

#' @describeIn GeneRegionTrack-class Initialize.
#' @export
setMethod("initialize", "GeneRegionTrack", function(.Object, start, end, ...) {
    if (is.null(list(...)$range) && is.null(list(...)$genome) && is.null(list(...)$chromosome)) {
        return(.Object)
    }
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "GeneRegionTrack")
    .Object@start <- ifelse(is.null(start), 0, start)
    .Object@end <- ifelse(is.null(end), 0, end)
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## ReferenceGeneRegionTrack Class --------------------------------------------

## The file-based version of the GeneRegionTrack class. This will mainly provide a means to dispatch to
## a special 'subset' method which should stream the necessary data from disk.

#' @describeIn GeneRegionTrack-class The file-based version of the `GeneRegionTrack-class`.
#' @exportClass ReferenceGeneRegionTrack
setClass("ReferenceGeneRegionTrack", contains = c("GeneRegionTrack", "ReferenceTrack"))

## Initialize ----------------------------------------------------------------


## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritence has some strange features with regards to method selection

#' @describeIn GeneRegionTrack-class Initialize.
#' @export
setMethod("initialize", "ReferenceGeneRegionTrack", function(.Object, stream, reference, mapping = list(),
                                                             args = list(), defaults = list(), ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object = .Object, reference = reference, stream = stream,
        mapping = mapping, args = args, defaults = defaults)
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

## Constructor. The following arguments are supported:
##    o range: one in a whole number of different potential inputs upon which the .buildRanges method will dispatch
##    o start, end: numeric vectors of the track start and end coordinates.
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o rstarts, rends, rwidths: integer vectors of exon start and end locations or widths, or a character vector of
##      comma-delimited exon locations, one vector element for each transcript
##    o strand, feature, exon, transcript, gene, symbol, chromosome: vectors of equal length containing
##       the exon strand, biotype, exon id, transcript id, gene id, human-readable gene symbol and chromosome information
##    o stacking: character controlling the stacking of overlapping items. One in 'hide',
##       'dense', 'squish', 'pack' or 'full'.
##    o name: the name of the track. This will be used for the title panel.
##    o exonList: boolean, causing the values in starts, rends or rwidths to be interpreted as delim-separated
##       lists that have to be exploded. All other annotation arguments will be repeated accordingly.
##    o delim: the delimiter if coordinates are in a list
## All additional items in ... are being treated as DisplayParameters

#' @describeIn GeneRegionTrack-class Constructor function for `GeneRegionTrack-class`.
#' @export
GeneRegionTrack <- function(range = NULL, rstarts = NULL, rends = NULL, rwidths = NULL, strand, feature, exon,
                            transcript, gene, symbol, chromosome, genome, stacking = "squish",
                            name = "GeneRegionTrack", start = NULL, end = NULL, importFunction, stream = FALSE, ...) {
    ## Some defaults
    covars <- if (is.data.frame(range)) range else if (is(range, "GRanges")) as.data.frame(mcols(range)) else data.frame()
    isStream <- FALSE
    if (!is.character(range)) {
        n <- if (is.null(range)) max(c(length(start), length(end), length(width))) else if (is(range, "data.frame")) nrow(range) else length(range)
        if (is.null(covars[["feature"]]) && missing(feature)) {
            feature <- paste("exon", seq_len(n), sep = "_")
        }
        if (is.null(covars[["exon"]]) && missing(exon)) {
            exon <- make.unique(rep(if (!missing(feature) && !is.null(feature)) as.character(feature) else covars[["feature"]], n)[seq_len(n)])
        }
        if (is.null(covars[["transcript"]]) && missing(transcript)) {
            transcript <- paste("transcript", seq_len(n), sep = "_")
        }
        if (is.null(covars[["gene"]]) && missing(gene)) {
            gene <- paste("gene", seq_len(n), sep = "_")
        }
    }
    ## Build a GRanges object from the inputs
    .missingToNull(c("feature", "exon", "transcript", "gene", "symbol", "strand", "chromosome", "importFunction", "genome"))
    args <- list(
        feature = feature, id = exon, exon = exon, transcript = transcript, gene = gene, symbol = symbol, strand = strand,
        chromosome = chromosome, genome = genome
    )
    defs <- list(
        feature = "unknown", id = "unknown", exon = "unknown", transcript = "unknown", genome = NA,
        gene = "unknown", symbol = "unknown", strand = "*", density = 1, chromosome = "chrNA"
    )
    range <- .buildRange(
        range = range, groupId = "transcript", start = rstarts, end = rends, width = rwidths, args = args, defaults = defs,
        chromosome = chromosome, tstart = start, tend = end, trackType = "GeneRegionTrack", importFun = importFunction,
        genome = genome
    )
    if (is.list(range)) {
        isStream <- TRUE
        slist <- range
        range <- GRanges()
    }
    if (is.null(start)) {
        start <- if (!length(range)) NULL else min(start(range))
    }
    if (is.null(end)) {
        end <- if (!length(range)) NULL else max(end(range))
    }
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- if (length(range) > 0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    }
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    if (!isStream) {
        return(new("GeneRegionTrack",
            start = start, end = end, chromosome = chromosome[1], range = range,
            name = name, genome = genome, stacking = stacking, ...
        ))
    } else {
        return(new("ReferenceGeneRegionTrack",
            start = start, end = end, chromosome = chromosome[1], range = range,
            name = name, genome = genome, stacking = stacking, stream = slist[["stream"]],
            reference = slist[["reference"]], mapping = slist[["mapping"]], args = args, defaults = defs, ...
        ))
    }
}

## .buildRange ---------------------------------------------------------------
##
## Helper methods to build a GRanges object from the input arguments.
##
## Coordinates for grouped elements may be passed in as comma-separated values
## (e.g. "1,5,9"), in which case we need to split and convert to numeric.
## This also implies that the additional arguments (like feature, group, etc.)
## have to be replicated accordingly. We handle this by passing along the repeat
## vector 'by' to the numeric method below.


## For TxDb objects we extract the grouping information and use the GRanges method

#' @importClassesFrom GenomicFeatures TxDb
#' @importMethodsFrom GenomicFeatures isActiveSeq "isActiveSeq<-" cdsBy exonsBy fiveUTRsByTranscript threeUTRsByTranscript transcriptsBy transcripts
setMethod(
    ".buildRange", signature("TxDb"),
    function(range, groupId = "transcript", tstart, tend, chromosome, args, ...) {
        ## If chromosome (and optional start and end) information is present we only extract parts of the annotation data
        noSubset <- is.null(tstart) && is.null(tend)
        if (!is.null(chromosome)) {
            chromosome <- .chrName(chromosome)
            ## Seems like TxDb objects use pass by reference for the active chromosomes, so we have to
            ## restore the old values after we are done
            oldAct <- seqlevels(range)
            oldRange <- range
            on.exit({
                restoreSeqlevels(oldRange)
                seqlevels(oldRange, pruning.mode = "coarse") <- oldAct
            })
            restoreSeqlevels(range)
            seqlevels(range, pruning.mode = "coarse") <- chromosome
            sl <- seqlengths(range)
            if (is.null(tstart)) {
                tstart <- rep(1, length(chromosome))
            }
            if (is.null(tend)) {
                tend <- sl[chromosome] + 1
                tend[is.na(tend)] <- tstart[is.na(tend)] + 1
            }
            sRange <- GRanges(seqnames = chromosome, ranges = IRanges(start = tstart, end = tend))
        }
        ## First the mapping of internal transcript ID to transcript name
        txs <- as.data.frame(values(transcripts(range, columns = c("tx_id", "tx_name"))))
        rownames(txs) <- txs[, "tx_id"]
        ## Now the CDS ranges
        t2c <- cdsBy(range, "tx")
        names(t2c) <- txs[names(t2c), 2]
        tids <- rep(names(t2c), elementNROWS(t2c))
        t2c <- unlist(t2c)
        if (length(t2c)) {
            t2c$tx_id <- tids
            t2c$feature_type <- "CDS"
        }
        ## And the 5'UTRS
        t2f <- fiveUTRsByTranscript(range)
        names(t2f) <- txs[names(t2f), 2]
        tids <- rep(names(t2f), elementNROWS(t2f))
        t2f <- unlist(t2f)
        if (length(t2f)) {
            t2f$tx_id <- tids
            t2f$feature_type <- "utr5"
        }
        ## And the 3'UTRS
        t2t <- threeUTRsByTranscript(range)
        names(t2t) <- txs[names(t2t), 2]
        tids <- rep(names(t2t), elementNROWS(t2t))
        t2t <- unlist(t2t)
        if (length(t2t)) {
            t2t$tx_id <- tids
            t2t$feature_type <- "utr3"
        }
        ## And finally all the non-coding transcripts
        nt2e <- exonsBy(range, "tx")
        names(nt2e) <- txs[names(nt2e), 2]
        nt2e <- nt2e[!names(nt2e) %in% c(values(t2c)$tx_id, values(t2f)$tx_id, values(t2t)$tx_id)]
        tids <- rep(names(nt2e), elementNROWS(nt2e))
        nt2e <- unlist(nt2e)
        if (length(nt2e)) {
            nt2e$tx_id <- tids
            nt2e$feature_type <- "ncRNA"
        }
        ## Now we can merge the three back together (we need to change the column names of t2c to make them all the same)
        colnames(values(t2c))[c(1, 2)] <- c("exon_id", "exon_name")
        ## t2e <- c(t2c, t2f, t2t, nt2e) ## This is super-slow, much more efficient if we build the GRanges object from the individual bits and pieces
        vals <- DataFrame(
            exon_id = c(values(t2c)$exon_id, values(t2f)$exon_id, values(t2t)$exon_id, values(nt2e)$exon_id),
            exon_name = c(values(t2c)$exon_name, values(t2f)$exon_name, values(t2t)$exon_name, values(nt2e)$exon_name),
            exon_rank = c(values(t2c)$exon_rank, values(t2f)$exon_rank, values(t2t)$exon_rank, values(nt2e)$exon_rank),
            tx_id = c(values(t2c)$tx_id, values(t2f)$tx_id, values(t2t)$tx_id, values(nt2e)$tx_id),
            feature_type = c(values(t2c)$feature_type, values(t2f)$feature_type, values(t2t)$feature_type, values(nt2e)$feature_type)
        )
        t2e <- GRanges(
            seqnames = c(seqnames(t2c), seqnames(t2f), seqnames(t2t), seqnames(nt2e)),
            ranges = IRanges(
                start = c(start(t2c), start(t2f), start(t2t), start(nt2e)),
                end = c(end(t2c), end(t2f), end(t2t), end(nt2e))
            ),
            strand = c(strand(t2c), strand(t2f), strand(t2t), strand(nt2e))
        )
        if (all(is.na(vals$exon_name))) {
            vals$exon_name <- paste(vals$tx_id, vals$exon_rank, sep = "_")
        }
        values(t2e) <- vals
        if (length(t2e) == 0) {
            return(GRanges())
        }
        ## Add the gene level annotation
        g2t <- transcriptsBy(range, "gene")
        gids <- rep(names(g2t), elementNROWS(g2t))
        g2t <- unlist(g2t)
        values(g2t)[["gene_id"]] <- gids
        values(t2e)$gene_id <- gids[match(values(t2e)$tx_id, as.character(txs[as.character(values(g2t)$tx_id), 2]))]
        vals <- values(t2e)[c("tx_id", "exon_name", "exon_rank", "feature_type", "tx_id", "gene_id")]
        colnames(vals) <- c("transcript", "exon", "rank", "feature", "symbol", "gene")
        ## Add the genome information
        genome(t2e) <- unique(genome(range))
        ## Finally we re-assign, subset if necessary, and sort
        range <- t2e
        values(range) <- vals
        if (!noSubset && !is.null(chromosome)) {
            ## We have to keep all exons for all the overlapping transcripts
            txSel <- unique(subsetByOverlaps(g2t, sRange)$tx_name)
            range <- range[range$transcript %in% txSel]
        }
        args <- list(genome = genome(range)[1])
        return(.buildRange(range = sort(range), chromosome = chromosome, args = args, ...))
    }
)

## For ensDb objects we extract the grouping information and use the GRanges method

#' @importClassesFrom ensembldb EnsDb
#' @importMethodsFrom ensembldb cdsBy exonsBy fiveUTRsByTranscript threeUTRsByTranscript transcriptsBy transcripts
setMethod(
    ".buildRange", signature("EnsDb"),
    function(range, groupId = "transcript", tstart, tend, chromosome, args, ...) {
        ## If chromosome (and optional start and end) information is present we only extract parts of the annotation data
        noSubset <- is.null(tstart) && is.null(tend)
        if (!is.null(chromosome)) {
            chromosome <- .chrName(chromosome)
            sl <- seqlengths(range)
            if (is.null(tstart)) {
                tstart <- rep(1, length(chromosome))
            }
            if (is.null(tend)) {
                tend <- sl[chromosome] + 1
                tend[is.na(tend)] <- tstart[is.na(tend)] + 1
            }
            sRange <- GRanges(seqnames = chromosome, ranges = IRanges(start = tstart, end = tend))
            ## sRange <- GRangesFilter(sRange, type = "any") can we filter directly?
        }
        ## First the mapping of internal transcript ID to transcript name
        txs <- as.data.frame(values(transcripts(range, columns = c("tx_id", "tx_name"))))
        rownames(txs) <- txs[, "tx_id"]
        ## Now the CDS ranges
        t2c <- cdsBy(range, "tx")
        names(t2c) <- txs[names(t2c), 2]
        tids <- rep(names(t2c), elementNROWS(t2c))
        t2c <- unlist(t2c)
        if (length(t2c)) {
            t2c$tx_id <- tids
            t2c$feature_type <- "CDS"
        }
        ## And the 5'UTRS
        t2f <- fiveUTRsByTranscript(range)
        names(t2f) <- txs[names(t2f), 2]
        tids <- rep(names(t2f), elementNROWS(t2f))
        t2f <- unlist(t2f)
        if (length(t2f)) {
            t2f$tx_id <- tids
            t2f$feature_type <- "utr5"
        }
        ## And the 3'UTRS
        t2t <- threeUTRsByTranscript(range)
        names(t2t) <- txs[names(t2t), 2]
        tids <- rep(names(t2t), elementNROWS(t2t))
        t2t <- unlist(t2t)
        if (length(t2t)) {
            t2t$tx_id <- tids
            t2t$feature_type <- "utr3"
        }
        ## And finally all the non-coding transcripts
        nt2e <- exonsBy(range, "tx")
        names(nt2e) <- txs[names(nt2e), 2]
        nt2e <- nt2e[!names(nt2e) %in% c(values(t2c)$tx_id, values(t2f)$tx_id, values(t2t)$tx_id)]
        tids <- rep(names(nt2e), elementNROWS(nt2e))
        nt2e <- unlist(nt2e)
        if (length(nt2e)) {
            nt2e$tx_id <- tids
            nt2e$feature_type <- "ncRNA"
        }
        ## Now we can merge the three back together (we need to change the column names of t2c to make them all the same)
        ## t2e <- c(t2c, t2f, t2t, nt2e) ## This is super-slow, much more efficient if we build the GRanges object from the individual bits and pieces
        vals <- DataFrame(
            exon_id = c(values(t2c)$exon_id, values(t2f)$exon_id, values(t2t)$exon_id, values(nt2e)$exon_id),
            exon_name = c(values(t2c)$exon_id, values(t2f)$exon_id, values(t2t)$exon_id, values(nt2e)$exon_id),
            exon_rank = c(values(t2c)$exon_rank, values(t2f)$exon_rank, values(t2t)$exon_rank, values(nt2e)$exon_rank),
            tx_id = c(values(t2c)$tx_id, values(t2f)$tx_id, values(t2t)$tx_id, values(nt2e)$tx_id),
            feature_type = c(values(t2c)$feature_type, values(t2f)$feature_type, values(t2t)$feature_type, values(nt2e)$feature_type)
        )
        t2e <- GRanges(
            seqnames = c(seqnames(t2c), seqnames(t2f), seqnames(t2t), seqnames(nt2e)),
            ranges = IRanges(
                start = c(start(t2c), start(t2f), start(t2t), start(nt2e)),
                end = c(end(t2c), end(t2f), end(t2t), end(nt2e))
            ),
            strand = c(strand(t2c), strand(t2f), strand(t2t), strand(nt2e))
        )
        if (all(is.na(vals$exon_name))) {
            vals$exon_name <- paste(vals$tx_id, vals$exon_rank, sep = "_")
        }
        values(t2e) <- vals
        if (length(t2e) == 0) {
            return(GRanges())
        }
        ## Add the gene level annotation
        g2t <- transcriptsBy(range, "gene")
        gids <- rep(names(g2t), elementNROWS(g2t))
        g2t <- unlist(g2t)
        values(g2t)[["gene_id"]] <- gids
        values(t2e)$gene_id <- gids[match(values(t2e)$tx_id, as.character(txs[as.character(values(g2t)$tx_id), 2]))]
        vals <- values(t2e)[c("tx_id", "exon_name", "exon_rank", "feature_type", "tx_id", "gene_id")]
        colnames(vals) <- c("transcript", "exon", "rank", "feature", "symbol", "gene")
        ## Add the genome information
        genome(t2e) <- unique(genome(range))
        ## Finally we re-assign, subset if necessary, and sort
        range <- t2e
        values(range) <- vals
        if (!noSubset && !is.null(chromosome)) {
            ## We have to keep all exons for all the overlapping transcripts
            txSel <- unique(subsetByOverlaps(g2t, sRange)$tx_name)
            range <- range[range$transcript %in% txSel]
        }
        args <- list(genome = genome(range)[1])
        return(.buildRange(range = sort(range), chromosome = chromosome, args = args, ...))
    }
)

## General accessors ---------------------------------------------------------
## Annotation Accessors ------------------------------------------------------

#' @describeIn GeneRegionTrack-class Extract the gene identifiers for all
#' gene models.
#' @export
setMethod("gene", signature(GdObject = "GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "gene"))

#' @describeIn GeneRegionTrack-class Replace the gene identifiers for all
#' gene models.
#' The replacement value must be a character of appropriate length or another
#' vector that can be coerced into such.
#' @export
setReplaceMethod("gene", signature("GeneRegionTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "gene"))

#' @describeIn GeneRegionTrack-class Extract the human-readble gene symbol
#' for all gene models.
#' @export
setMethod("symbol", signature(GdObject = "GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "symbol"))

#' @describeIn GeneRegionTrack-class Replace the human-readable gene symbol
#' for all gene models.
#' The replacement value must be a character of appropriate length or another
#' vector that can be coerced into such.
#' @export
setReplaceMethod("symbol", signature("GeneRegionTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "symbol"))

#' @describeIn GeneRegionTrack-class Extract the transcript identifiers for all
#' transcripts in the gene models.
#' @export
setMethod("transcript", signature(GdObject = "GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "transcript"))

#' @describeIn GeneRegionTrack-class Replace the transcript identifiers for all
#' transcripts in the gene model. The replacement value must be a character of
#' appropriate length or another vector that can be coerced into such.
#' @export
setReplaceMethod(
    "transcript", signature("GeneRegionTrack", "character"),
    function(GdObject, value) .setAnn(GdObject, value, "transcript")
)

#' @describeIn GeneRegionTrack-class Extract the exon identifiers for all exons
#' in the gene models.
#' @export
setMethod("exon", signature(GdObject = "GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "exon"))

#' @describeIn GeneRegionTrack-class replace the exon identifiers for all exons
#' in the gene model. The replacement value must be a character of appropriate
#' length or another vector that can be coerced into such.
#' @export
setReplaceMethod("exon", signature("GeneRegionTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "exon"))

#' @describeIn GeneRegionTrack-class extract the group membership for all track items.
#' @export
setMethod("group", "GeneRegionTrack", function(GdObject) transcript(GdObject))

#' @describeIn GeneRegionTrack-class replace the grouping information for track items.
#' The replacement value must be a factor of appropriate length or another
#' vector that can be coerced into such.
#' @export
setReplaceMethod(
    "group", signature("GeneRegionTrack", "character"),
    function(GdObject, value) .setAnn(GdObject, value, "transcript")
)

#' @describeIn GeneRegionTrack-class return track item identifiers.
#' Depending on the setting of the optional argument lowest, these are either
#' the group identifiers or the individual item identifiers.
#' export
setMethod("identifier", "GeneRegionTrack", function(GdObject, type = .dpOrDefault(GdObject, "transcriptAnnotation", "symbol")) {
    if (is.logical(type)) {
        type <- ifelse("symbol", "gene", type[1])
    }
    if (is.null(type)) {
        type <- "symbol"
    }
    id <- switch(as.character(type),
        "symbol" = symbol(GdObject),
        "gene" = gene(GdObject),
        "transcript" = transcript(GdObject),
        "feature" = feature(GdObject),
        "exon" = exon(GdObject),
        "lowest" = exon(GdObject),
        symbol(GdObject)
    )
    id[is.na(id)] <- "NA"
    return(id)
})

#' @describeIn GeneRegionTrack-class Set the track item identifiers.
#' The replacement value has to be a character vector of appropriate length.
#' This always replaces the group-level identifiers, so essentially it is
#' similar to `groups<-`.
#' @export
setReplaceMethod("identifier", c("GeneRegionTrack", "character"), function(GdObject, value) {
    type <- .dpOrDefault(GdObject, "transcriptAnnotation", "symbol")
    switch(as.character(type),
        "symbol" = symbol(GdObject) <- value,
        "gene" = gene(GdObject) <- value,
        "transcript" = transcript(GdObject) <- value,
        "feature" = feature(GdObject) <- value,
        "exon" = exon(GdObject) <- value,
        "lowest" = exon(GdObject) <- value,
        symbol(GdObject) <- value
    )
    return(GdObject)
})

## Stacking ------------------------------------------------------------------
## Consolidate ---------------------------------------------------------------
## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------

## FIXME: Still needs to be implemented

#' @describeIn GeneRegionTrack-class Subset a GeneRegionTrack by coordinates
#' and sort if necessary.
#' @export
setMethod("subset", signature(x = "ReferenceGeneRegionTrack"), function(x, ...) {
    warning("ReferenceGeneRegionTrack objects are not supported yet.")
    return(callNextMethod())
})

## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------
## DrawAxis ------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------

#' @describeIn GeneRegionTrack-class plot the object to a graphics device.
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
setMethod("drawGD", signature("GeneRegionTrack"), function(GdObject, ...) {
    displayPars(GdObject) <- list(showFeatureId = as.vector(.dpOrDefault(GdObject, "showExonId")))
    GdObject <- callNextMethod()
    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------

#' @importFrom rtracklayer GenomicData
setAs(
    "GeneRegionTrack", "UCSCData",
    function(from, to) {
        ranges <- cbind(as(as(ranges(from), "DataFrame"), "data.frame"),
            start = start(from),
            end = end(from), color = .getBiotypeColor(from), strand = strand(from)
        )
        ranges <- ranges[order(start(from)), ]
        ranges <- split(ranges, ranges[, "X.transcript"])
        start <- vapply(ranges, function(x) min(x$start), FUN.VALUE = numeric(1L))
        end <- vapply(ranges, function(x) max(x$end), FUN.VALUE = numeric(1L))
        name <- vapply(ranges, function(x) as.character(unique(x$X.symbol)), FUN.VALUE = character(1L))
        color <- vapply(ranges, function(x) as.character(unique(x$color)), FUN.VALUE = character(1L))
        strand <- vapply(ranges, function(x) as.character(unique(x$strand)), FUN.VALUE = character(1L))
        strand[strand == "*"] <- "+"
        id <- names(ranges)
        blocks <- vapply(ranges, nrow, FUN.VALUE = numeric(1L))
        bsizes <- vapply(ranges, function(x) paste(x$end - x$start + 1, collapse = ","), FUN.VALUE = character(1L))
        bstarts <- vapply(ranges, function(x) paste(x$start - min(x$start), collapse = ","), FUN.VALUE = character(1L))
        dcolor <- as.integer(col2rgb(.dpOrDefault(from, "col")))
        line <- new("BasicTrackLine",
            name = names(from),
            description = names(from),
            visibility = stacking(from), color = dcolor, itemRgb = TRUE
        )
        new("UCSCData", rtracklayer::GenomicData(IRanges(start, end),
            chrom = chromosome(from),
            id = id, name = name, itemRgb = color, blockCount = blocks,
            blockSizes = bsizes, blockStarts = bstarts,
            strand = strand
        ),
        trackLine = line
        )
    }
)

setAs("GRanges", "GeneRegionTrack", function(from, to) GeneRegionTrack(range = from))

setAs("GRangesList", "GeneRegionTrack", function(from, to) GeneRegionTrack(range = from))

setAs("TxDb", "GeneRegionTrack", function(from, to) GeneRegionTrack(range = from))

## Show ----------------------------------------------------------------------

#' @describeIn GeneRegionTrack-class Show method.
#' @export
setMethod("show", signature(object = "GeneRegionTrack"), function(object) {
    cat(sprintf("GeneRegionTrack '%s'\n%s\n", names(object), .annotationTrackInfo(object)))
})

#' @describeIn GeneRegionTrack-class Show method.
#' @export
setMethod("show", signature(object = "ReferenceGeneRegionTrack"), function(object) {
    .referenceTrackInfo(object, "ReferenceGeneRegionTrack")
})
