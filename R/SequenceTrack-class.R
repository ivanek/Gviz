#' @include GdObject-class.R
#' @include ReferenceTrack-class.R
NULL

## SequenceTrack Class -------------------------------------------------------

#' SequenceTrack class and methods
#'
#'
#' A track class to represent genomic sequences. The three child classes
#' \code{SequenceDNAStringSetTrack}, \code{SequenceRNAStringSetTrack} and
#' \code{SequenceBSgenomeTrack} do most of the work, however in practise they
#' are of no particular relevance to the user.
#'
#'
#' @name SequenceTrack-class
#' @param sequence
#'
#' A meta argument to handle the different input types, making the construction
#' of a \code{SequenceTrack} as flexible as possible.
#'
#' The different input options for \code{sequence} are:
#'
#' \describe{
#'
#' \item{An object of class \code{\linkS4class{DNAStringSet}}.}{ The individual
#' \code{\linkS4class{DNAString}}s are considered to be the different
#' chromosome sequences.}
#'
#' \item{An object of class \code{\linkS4class{BSgenome}}.}{ The \code{Gviz}
#' package tries to follow the \code{BSgenome} philosophy in that the
#' respective chromosome sequences are only realized once they are first
#' accessed.}
#'
#' \item{A \code{character} scalar:}{ in this case the value of the
#' \code{sequence} argument is considered to be a file path to an annotation
#' file on disk. A range of file types are supported by the \code{Gviz} package
#' as identified by the file extension. See the \code{importFunction}
#' documentation below for further details.}
#' }
#'
#' @template SequenceTrack-class_param
#'
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{SequenceDNAStringSetTrack}, \code{SequenceBSgenomeTrack} ore
#' \code{ReferenceSequenceTrack}, depending on the constructor arguments.
#' Typically the user will not have to be troubled with this distinction and
#' can rely on the constructor to make the right choice.
#' @section Objects from the class:
#'
#' Objects can be created using the constructor function \code{SequenceTrack}.
#' @author Florian Hahne
#' @inherit GdObject-class seealso
#'
#' @examples
#' ## An empty object
#' SequenceTrack()
#'
#' ## Construct from DNAStringSet
#' library(Biostrings)
#' letters <- c("A", "C", "T", "G", "N")
#' set.seed(999)
#' seqs <- DNAStringSet(c(chr1 = paste(sample(letters, 100000, TRUE),
#'     collapse = ""
#' ), chr2 = paste(sample(letters, 200000, TRUE), collapse = "")))
#' sTrack <- SequenceTrack(seqs, genome = "hg19")
#' sTrack
#'
#' ## Construct from BSGenome object
#' if (require(BSgenome.Hsapiens.UCSC.hg19)) {
#'     sTrack <- SequenceTrack(Hsapiens)
#'     sTrack
#' }
#'
#'
#' ## Set active chromosome
#' chromosome(sTrack)
#' chromosome(sTrack) <- "chr2"
#' head(seqnames(sTrack))
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(sTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#' ## Plotting
#' ## Sequences
#' plotTracks(sTrack, from = 199970, to = 200000)
#' ## Boxes
#' plotTracks(sTrack, from = 199800, to = 200000)
#' ## Line
#' plotTracks(sTrack, from = 1, to = 200000)
#' ## Force boxes
#' plotTracks(sTrack, from = 199970, to = 200000, noLetters = TRUE)
#' ## Direction indicator
#' plotTracks(sTrack, from = 199970, to = 200000, add53 = TRUE)
#' ## Sequence complement
#' plotTracks(sTrack, from = 199970, to = 200000, add53 = TRUE, complement = TRUE)
#' ## Colors
#' plotTracks(sTrack, from = 199970, to = 200000, add53 = TRUE, fontcolor = c(
#'     A = 1,
#'     C = 1, G = 1, T = 1, N = 1
#' ))
#'
#' ## Track names
#' names(sTrack)
#' names(sTrack) <- "foo"
#'
#' ## Accessors
#' genome(sTrack)
#' genome(sTrack) <- "mm9"
#' length(sTrack)
#'
#' ## Sequence extraction
#' subseq(sTrack, start = 100000, width = 20)
#' ## beyond the stored sequence range
#' subseq(sTrack, start = length(sTrack), width = 20)
#' @importClassesFrom Biostrings DNAStringSet RNAStringSet BStringSet DNAString RNAString BString
#' @importClassesFrom BSgenome  BSgenome MaskedBSgenome
#' @importFrom Biostrings DNAStringSet RNAStringSet BStringSet DNAString RNAString BString reverseComplement readDNAStringSet DNA_ALPHABET stackStrings
#' @importFrom BSgenome bsapply
#' @importFrom biovizBase getBioColor
#'
#' @exportClass SequenceTrack
setClass("SequenceTrack",
    representation = representation("VIRTUAL",
        chromosome = "character",
        genome = "character"
    ),
    contains = "GdObject",
    prototype = prototype(
        name = "Sequence",
        dp = DisplayPars(
            add53 = FALSE,
            background.title = "transparent",
            col = "darkgray",
            complement = FALSE,
            fontface = 2,
            fontsize = 10,
            lwd = 2,
            min.width = 2,
            noLetters = FALSE,
            showTitle = FALSE,
            size = NULL
        ),
        genome = as.character(NA),
        chromosome = "chrNA"
    )
)

## Essentially we just update the display parameters here and set the chromosome and the genome

#' @describeIn  SequenceTrack-class Initialize.
#' @export
setMethod("initialize", "SequenceTrack", function(.Object, chromosome, genome, ...) {
    ## the display parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "SequenceTrack")
    if (!missing(chromosome) &&
        !is.na(chromosome) &&
        !is.null(chromosome)) {
        if (!is.null(names(sequence))) {
            .Object@chromosome <- .chrName(names(sequence)[1])[1]
        } else {
            .Object@chromosome <- .chrName(chromosome)[1]
        }
    }
    if (missing(genome) || is.na(genome) || is.null(genome)) {
        genome <- as.character(NA)
    }
    .Object@genome <- genome
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructors --------------------------------------------------------------

#' @importMethodsFrom Biostrings seqtype "seqtype<-"
.SequenceTrack <- function(SeqTrackType, seqtype, sequence, chromosome,
                           genome, name = "SequenceTrack", importFunction,
                           stream = FALSE, ...) {
    .missingToNull(c("chromosome", "genome", "sequence"))
    if (is.null(sequence)) {
        return(new(SeqTrackType, chromosome = chromosome, genome = genome, name = name, ...))
    }
    if (is(sequence, "BSgenome")) {
        if (is.null(genome)) {
            genome <- .providerVersion(sequence)
        }
        if (is.null(chromosome)) {
            chromosome <- seqnames(sequence)[1]
        }
        obj <- new("SequenceBSgenomeTrack", sequence = sequence, chromosome = chromosome, genome = genome, name = name, ...)
    } else if (is(sequence, "XStringSet")) {
        if (seqtype(sequence) != seqtype) {
            seqtype(sequence) <- seqtype
        }
        if (is.null(names(sequence))) {
            stop("The sequences in the ", seqtype, "StringSet must be named")
        }
        if (any(duplicated(names(sequence)))) {
            stop("The sequence names in the ", seqtype, "StringSet must be unique")
        }
        if (is.null(chromosome)) {
            chromosome <- names(sequence)[1]
        }
        obj <- new(SeqTrackType, sequence = sequence, chromosome = chromosome, genome = genome, name = name, ...)
    } else if (is.character(sequence)) {
        sequence <- sequence[1]
        if (!file.exists(sequence)) {
            stop(sprintf("'%s' is not a valid file.", sequence))
        }
        ext <- .fileExtension(sequence)
        obj <- if (missing(importFunction) && ext %in% c("fa", "fasta")) {
            if (!file.exists(paste(sequence, "fai", sep = "."))) {
                new("SequenceDNAStringSetTrack",
                    sequence = readDNAStringSet(sequence), chromosome = chromosome,
                    genome = genome, name = name, ...
                )
            } else {
                new("ReferenceSequenceTrack",
                    chromosome = chromosome, genome = genome, name = name,
                    stream = .import.fasta, reference = path.expand(sequence), ...
                )
            }
        } else if (missing(importFunction) && ext == "2bit") {
            new("ReferenceSequenceTrack",
                chromosome = chromosome, genome = genome, name = name,
                stream = .import.2bit, reference = path.expand(sequence), ...
            )
        } else {
            if (missing(importFunction)) {
                stop(sprintf(
                    "No predefined import function exists for files with extension '%s'. Please manually provide an import function.",
                    ext
                ))
            } else {
                if (!stream) {
                    seq <- importFunction(file = sequence)
                    if (!is(seq, "DNAStringSet")) {
                        stop(
                            "The import function did not provide a valid DNAStringSet object. Unable to build track from file '",
                            sequence, "'"
                        )
                    }
                    new("SequenceDNAStringSetTrack",
                        sequence = importFunction(file = sequence), chromosome = chromosome,
                        genome = genome, name = name, ...
                    )
                } else {
                    new("ReferenceSequenceTrack",
                        chromosome = chromosome, genome = genome, name = name,
                        stream = importFunction, reference = path.expand(sequence), ...
                    )
                }
            }
        }
    } else {
        stop("Argument sequence must be of class 'BSgenome', 'XStringSet' or 'character'")
    }
    return(obj)
}

#' @describeIn SequenceTrack-class Constructor
#' @export
SequenceTrack <- function(sequence, chromosome, genome, name = "SequenceTrack", importFunction, stream = FALSE, ...) {
    .SequenceTrack(
        "SequenceDNAStringSetTrack",
        "DNA",
        sequence,
        chromosome,
        genome,
        ...
    )
}

#' @describeIn SequenceTrack-class Constructor
#' @export
RNASequenceTrack <- function(sequence, chromosome, genome, name = "SequenceTrack", importFunction, stream = FALSE, ...) {
    .SequenceTrack(
        "SequenceRNAStringSetTrack",
        "RNA",
        sequence,
        chromosome,
        genome,
        ...
    )
}

## SequenceDNAStringSetTrack -------------------------------------------------

## We want the following behaviour in the constructor:
##   a) sequence is missing (NULL) => build SequenceDNAStringSetTrack with chromosome NA and genome as supplied or NA if missing
##   b) sequence is DNAStringSet => build SequenceDNAStringSetTrack where chromosome is names(sequence)[1] or the supplied
##      chromosome if available, and genome as supplied or NA if missing
##   c) sequence is BSgenome => build SequenceBSgenomeTrack where chromosome is seqnames(sequence)[1] or the supplied
##      chromosome if available, and genome is the supplied genome or the one extracted from the BSgenome object

##
## SequenceDNAStringSetTrack: A track to visualize nucleotide sequences that are stored in a DNAStringSet
##
## Slots:
##    o sequence: a DNAStringSet object that contains all the sequence data

#' @describeIn SequenceTrack-class The `DNAStringSet`-based version of the `SequenceTrack-class`.
#' @exportClass SequenceDNAStringSetTrack
setClass("SequenceDNAStringSetTrack",
    representation = representation(sequence = "DNAStringSet"),
    contains = "SequenceTrack",
    prototype = prototype(
        sequence = Biostrings::DNAStringSet(),
        dp = DisplayPars(fontcolor = biovizBase::getBioColor("DNA_BASES_N"))
    )
)

#' @describeIn  SequenceTrack-class Initialize.
#' @export
setMethod("initialize", "SequenceDNAStringSetTrack", function(.Object, sequence, ...) {
    if (missing(sequence) || is.null(sequence)) {
        sequence <- Biostrings::DNAStringSet()
    }
    .Object@sequence <- sequence
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})


## SequenceRNAStringSetTrack -------------------------------------------------
##
## SequenceRNAStringSetTrack: A track to visualize nucleotide sequences that are stored in a RNAStringSet
## Slots:
##    o sequence: a RNAStringSet object that contains all the sequence data

#' @describeIn SequenceTrack-class The `RNAStringSet`-based version of the `SequenceTrack-class`.
#' @exportClass SequenceRNAStringSetTrack
setClass("SequenceRNAStringSetTrack",
    representation = representation(sequence = "RNAStringSet"),
    contains = "SequenceTrack",
    prototype = prototype(
        sequence = Biostrings::RNAStringSet(),
        dp = DisplayPars(fontcolor = biovizBase::getBioColor("RNA_BASES_N"))
    )
)

#' @describeIn SequenceTrack-class Initialize `RNAStringSet`-based version of the `SequenceTrack-class`.
#' @export
setMethod("initialize", "SequenceRNAStringSetTrack", function(.Object, sequence, ...) {
    if (missing(sequence) || is.null(sequence)) {
        sequence <- Biostrings::RNAStringSet()
    }
    .Object@sequence <- sequence
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})


## SequenceBSgenomeTrack -----------------------------------------------------
##
## SequenceBSgenomeTrack: A track to visualize nucleotide sequences that are stored in a BSgenome package
##
## Slots:
##    o sequence: a DNAStringSet object that contains all the sequence data
##    o pointerCache: an environment to hold pointers to the BSgenome sequences to prevent garbage collection. This
##       will only be filled once the individual sequences have been accessed for the first time

#' @describeIn SequenceTrack-class The `BSgenome`-based version of the `SequenceTrack-class`.
#' @exportClass SequenceBSgenomeTrack
setClass("SequenceBSgenomeTrack",
    representation = representation(sequence = "BSgenomeOrNULL", pointerCache = "environment"),
    contains = "SequenceTrack",
    prototype = prototype(
        sequence = NULL,
        dp = DisplayPars(fontcolor = biovizBase::getBioColor("DNA_BASES_N"))
    )
)

#' @describeIn  SequenceTrack-class Initialize.
#' @export
setMethod("initialize", "SequenceBSgenomeTrack", function(.Object, sequence = NULL, ...) {
    .Object@sequence <- sequence
    .Object@pointerCache <- new.env()
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})



## ReferenceSequenceTrack ----------------------------------------------------
##
## ReferenceSequenceTrack: The file-based version of the ReferenceTrack class.
##                         This will mainly provide a means to dispatch to a special 'subseq' method
##                         which should stream the necessary data from disk.

#' @describeIn SequenceTrack-class The file-based version of the `SequenceTrack-class`.
#' @exportClass ReferenceSequenceTrack
setClass("ReferenceSequenceTrack", contains = c("SequenceDNAStringSetTrack", "ReferenceTrack"))

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritance has some strange features with regards to method selection

#' @describeIn  SequenceTrack-class Initialize.
#' @export
setMethod("initialize", "ReferenceSequenceTrack", function(.Object, stream, reference, ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object = .Object, reference = reference, stream = stream)
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Initialize ----------------------------------------------------------------
## Constructor ---------------------------------------------------------------
## General accessors ---------------------------------------------------------

#' @describeIn SequenceTrack-class return the names (i.e., the chromosome)
#' of the sequences contained in the object.
#' @export
setMethod("seqnames", "SequenceTrack", function(x) as.character(names(x@sequence)))

# setMethod("seqnames", "SequenceDNAStringSetTrack", function(x) as.character(names(x@sequence)))
# setMethod("seqnames", "SequenceRNAStringSetTrack", function(x) as.character(names(x@sequence)))

#' @describeIn SequenceTrack-class return the names (i.e., the chromosome)
#' of the sequences contained in the object.
#' @export
setMethod("seqnames", "SequenceBSgenomeTrack", function(x) as.character(seqnames(x@sequence)))

#' @describeIn SequenceTrack-class return the names (i.e., the chromosome)
#' of the sequences contained in the object. Only those with length > 0.
#' @export
setMethod("seqlevels", "SequenceTrack", function(x) seqnames(x)[width(x@sequence) > 0])

# setMethod("seqlevels", "SequenceDNAStringSetTrack", function(x) seqnames(x)[width(x@sequence) > 0])
# setMethod("seqlevels", "SequenceRNAStringSetTrack", function(x) seqnames(x)[width(x@sequence) > 0])

#' @describeIn SequenceTrack-class return the names (i.e., the chromosome)
#' of the sequences contained in the object. Only those with length > 0.
#' @export
setMethod("seqlevels", "SequenceBSgenomeTrack", function(x) seqnames(x)[BSgenome::bsapply(new("BSParams", X = x@sequence, FUN = length, simplify = TRUE)) > 0]) # maybe seqnames only, to speed-up

#' @describeIn SequenceTrack-class return the start coordinates of the track
#' items.
#' @export
setMethod("start", "SequenceTrack", function(x) NULL)

#' @describeIn SequenceTrack-class return the end coordinates of the track
#' items.
#' @export
setMethod("end", "SequenceTrack", function(x) NULL)

#' @describeIn SequenceTrack-class return the with of the track items in
#' genomic coordinates.
#' @export
setMethod("width", "SequenceTrack", function(x) NULL)

#' @describeIn SequenceTrack-class return the length of the sequence for
#' active chromosome.
#' @export
setMethod("length", "SequenceTrack", function(x) {
    if (chromosome(x) %in% seqnames(x)) length(x@sequence[[chromosome(x)]]) else 0
})

#' @importMethodsFrom Biostrings unmasked complement
setMethod("subseq", "SequenceTrack", function(x, start = NA, end = NA, width = NA) {
    padding <- "-"
    if (!is.na(start[1] + end[1] + width[1])) {
        warning("All 'start', 'stop' and 'width' are provided, ignoring 'width'")
        width <- NA
    }
    ## We want start and end to be set if width is provided
    if (!is.na(width[1])) {
        if (is.na(start) && is.na(end)) {
            stop("Two out of the three in 'start', 'end' and 'width' have to be provided")
        }
        if (is.na(start)) {
            start <- end - width[1] + 1
        }
        if (is.na(end)) {
            end <- start + width[1] - 1
        }
    }
    if (is.na(start)) {
        start <- 1
    }
    w <- length(x)
    if (w > 0) {
        if (is.na(end)) {
            end <- w
        }
        rstart <- max(1, start[1], na.rm = TRUE)
        rend <- max(rstart, min(end[1], w, na.rm = TRUE))
    } else {
        if (is.na(end)) {
            end <- start
        }
        rend <- end
        rstart <- start
    }
    if (rend < rstart || end < start) {
        stop("'end' has to be bigger than 'start'")
    }
    if ((rend - rstart + 1) > 10e6) {
        stop("Sequence is too big! Unable to extract")
    }
    seqtype <- try(seqtype(x@sequence), silent = TRUE)
    if (is(seqtype, "try-error")) {
        seqtype <- "DNA"
    }
    class <- paste0(seqtype, "String")
    finalSeq <- rep(do.call(class, list(padding)), end - start + 1)
    if (chromosome(x) %in% seqnames(x) && rend >= rstart) {
        chrSeq <- x@sequence[[chromosome(x)]]
        seq <- subseq(chrSeq, start = rstart, end = rend)
        if (is(x, "SequenceBSgenomeTrack")) seq <- unmasked(seq)
        subseq(finalSeq, ifelse(start < 1, abs(start) + 2, 1), width = rend - rstart + 1) <- seq
    }
    if (is(x, "SequenceBSgenomeTrack") && chromosome(x) %in% seqnames(x)) {
        x@pointerCache[[chromosome(x)]] <- x@sequence[[chromosome(x)]]
    }
    if (.dpOrDefault(x, "complement", FALSE)) {
        finalSeq <- complement(finalSeq)
    }
    return(finalSeq)
})


setMethod("subseq", "ReferenceSequenceTrack", function(x, start = NA, end = NA, width = NA) {
    if (sum(c(is.na(start[1]), is.na(end[1]), is.na(width[1]))) >= 2) {
        stop("Two out of the three in 'start', 'end' and 'width' have to be provided")
    }
    if (!is.na(start[1] + end[1] + width[1])) {
        warning("All 'start', 'stop' and 'width' are provided, ignoring 'width'")
        width <- NA
    }
    ## We want start and end to be set if width is provided
    if (!is.na(width[1])) {
        if (is.na(start) && is.na(end)) {
            stop("Two out of the three in 'start', 'end' and 'width' have to be provided")
        }
        if (is.na(start)) {
            start <- end - width[1] + 1
        }
        if (is.na(end)) {
            end <- start + width[1] - 1
        }
    }
    x@sequence <- x@stream(file = x@reference, selection = GRanges(chromosome(x), ranges = IRanges(1, end)))
    return(callNextMethod())
})

#' @describeIn SequenceTrack-class return the chromosome for which the track
#' is defined.
#' @export
setMethod("chromosome", "SequenceTrack", function(GdObject) GdObject@chromosome)

#' @describeIn SequenceTrack-class replace the value of the track's chromosome.
#' This has to be a valid UCSC chromosome identifier or an integer or character
#' scalar that can be reasonably coerced into one.
#' @export
setReplaceMethod("chromosome", "SequenceTrack", function(GdObject, value) {
    GdObject@chromosome <- .chrName(value[1])
    return(GdObject)
})

#' @describeIn SequenceTrack-class Set the track's genome.
#' Usually this has to be a valid UCSC identifier, however this is not
#' formally enforced here.
#' @export
setMethod("genome", "SequenceTrack", function(x) x@genome)


## Annotation Accessors ------------------------------------------------------
## Consolidate ---------------------------------------------------------------
## Before starting of the plotting operation there are a bunch of housekeeping task that should be performed on each
## track, and the mileage may vary between track types, hence we add a layer of abstraction here by using a method.
##
## Available arguments are:
##    o GdObject: the input track object
##    o chromosome: the currently active chromosome which may have to be set for a RangeTrack or a SequenceTrack object
##    o ...: additional arguments that are considered to be display parameters

## For all track types we want to update the display parameters

## For RangeTracks and SequenceTracks we want to set the chromosome

#' @describeIn SequenceTrack-class Consolidate/
#' Determine whether there is `chromosome` settings or not, and add this information.
#' @param GdObject the input track object
#' @param chromosome the currently active chromosome which may have to be set
#' for a `RangeTrack` or a `SequenceTrack` object
#' parameters
# #' @keywords internal
#' @export
setMethod("consolidateTrack", signature(GdObject = "SequenceTrack"), function(GdObject, chromosome, ...) {
    if (!is.null(chromosome)) {
        chromosome(GdObject) <- chromosome
    }
    GdObject <- callNextMethod(GdObject, ...)
    return(GdObject)
})
## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------

#' @describeIn SequenceTrack-class plot the object to a graphics device.
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
setMethod("drawGD", signature("SequenceTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if ((is.logical(debug) && debug) || debug == "prepare") {
        browser()
    }
    fcol <- .dpOrDefault(GdObject, "fontcolor", getBioColor("DNA_BASES_N"))
    cex <- max(0.3, .dpOrDefault(GdObject, "cex", 1))
    xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
    pushViewport(viewport(xscale = xscale, clip = TRUE, gp = .fontGp(GdObject, cex = cex)))
    if (prepare) {
        pres <- .pxResolution()
        nsp <- max(as.numeric(convertHeight(stringHeight(stringWidth(DNA_ALPHABET)), "native")))
        nsp <- nsp / pres["y"] * 2
        displayPars(GdObject) <- list("neededVerticalSpace" = nsp)
        popViewport(1)
        return(invisible(GdObject))
    }
    if ((is.logical(debug) && debug) || debug == "draw") {
        browser()
    }
    imageMap(GdObject) <- NULL
    delta <- maxBase - minBase
    if (delta == 0) {
        return(invisible(GdObject))
    }
    lwidth <- max(as.numeric(convertUnit(stringWidth(DNA_ALPHABET), "inches")))
    perLetter <- vpLocation()$isize["width"] / (maxBase - minBase + 1)
    diff <- .pxResolution(.dpOrDefault(GdObject, "min.width", 2), coord = "x")
    ## FIXME: Need to deal with sequences that are too long.
    if (diff > 1 || (maxBase - minBase + 1) >= 10e6) {
        grid.lines(
            x = unit(c(minBase, maxBase), "native"), y = 0.5,
            gp = gpar(
                col = .dpOrDefault(GdObject, "col", "darkgray"),
                lwd = .dpOrDefault(GdObject, "lwd", 2)
            )
        )
    } else {
        sequence <- as.character(as(subseq(GdObject, start = minBase, end = maxBase - 1), "Rle"))
        at <- seq((minBase + 0.5), maxBase - 1 + 0.5, by = 1) # to align sequence (letters) with ticks position
        sequence[sequence == "-"] <- ""
        if (perLetter < 0.5 && .dpOrDefault(GdObject, "add53", FALSE)) {
            sequence[c(1, length(sequence))] <- ""
        }
        col <- fcol[toupper(sequence)]
        if (lwidth < perLetter && !.dpOrDefault(GdObject, "noLetters", FALSE)) {
            grid.text(
                x = unit(at, "native"), y = 0.5, label = sequence, rot = .dpOrDefault(GdObject, "rotation", 0),
                gp = gpar(col = col)
            )
        } else {
            grid.rect(
                x = unit(at, "native"), y = 0.05, width = unit(1, "native"), height = 0.9,
                gp = gpar(fill = col, col = "white"), just = c(0.5, 0)
            )
        }
    }
    ## The direction indicators
    if (.dpOrDefault(GdObject, "add53", FALSE)) {
        if (.dpOrDefault(GdObject, "complement", FALSE)) {
            grid.text(label = expression("3'"), x = unit(minBase + 0.1, "native"), just = c(0, 0.5), gp = gpar(col = "#808080", cex = 0.8))
            grid.text(label = expression("5'"), x = unit(maxBase - 0.1, "native"), just = c(1, 0.5), gp = gpar(col = "#808080", cex = 0.8))
        } else {
            grid.text(label = expression("5'"), x = unit(minBase + 0.1, "native"), just = c(0, 0.5), gp = gpar(col = "#808080", cex = 0.8))
            grid.text(label = expression("3'"), x = unit(maxBase - 0.1, "native"), just = c(1, 0.5), gp = gpar(col = "#808080", cex = 0.8))
        }
    }
    popViewport(1)
    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------

setAs("DNAString", "Rle", function(from, to) Rle(strsplit(as.character(from), "")[[1]]))

setAs("RNAString", "Rle", function(from, to) Rle(strsplit(as.character(from), "")[[1]]))

## Show ----------------------------------------------------------------------

.sequenceTrackInfo <- function(object) {
    msg <- sprintf(
        paste("Sequence track '%s':\n",
            "| genome: %s\n",
            "| chromosomes: %s\n",
            "| active chromosome: %s (%s nucleotides)\n",
            sep = ""
        ),
        names(object),
        genome(object),
        length(seqnames(object)),
        chromosome(object),
        length(object)
    )
    if (length(seqnames(object)) > 1) {
        msg <- paste(msg, "Call seqnames() to list all available chromosomes\n",
            "Call chromosome()<- to change the active chromosome\n",
            sep = ""
        )
    }
    return(msg)
}

## We need to show the name, genome, information about the source BSgenome object as well as the currently active chromosome


#' @describeIn  SequenceTrack-class Show method.
#' @export
setMethod(
    "show", signature(object = "SequenceBSgenomeTrack"),
    function(object) {
        cat(.sequenceTrackInfo(object),
            sprintf(
                paste("Parent BSgenome object:\n",
                    "| organism: %s\n",
                    "| provider: %s\n",
                    "| provider version: %s\n",
                    "| release date: %s\n",
                    # "| release name: %s\n",
                    "| package name: %s\n",
                    sep = ""
                ),
                organism(object@sequence),
                provider(object@sequence),
                .providerVersion(object@sequence),
                releaseDate(object@sequence),
                # releaseName(object@sequence),
                object@sequence@pkgname
            ),
            sep = ""
        )
    }
)

## Here we only need the name, genome and currently active chromosome information
#' @describeIn  SequenceTrack-class Show method.
#' @export
setMethod("show", signature(object = "SequenceDNAStringSetTrack"), function(object) cat(.sequenceTrackInfo(object)))

#' @describeIn  SequenceTrack-class Show method.
#' @export
setMethod("show", signature(object = "SequenceRNAStringSetTrack"), function(object) cat(.sequenceTrackInfo(object)))

#' @describeIn  SequenceTrack-class Show method.
#' @export
setMethod("show", signature(object = "ReferenceSequenceTrack"), function(object) {
    .referenceTrackInfo(object, "ReferenceSequenceTrack")
})
