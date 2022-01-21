#' @include GdObject-class.R
NULL

## ReferenceTrack-class ------------------------------------------------------

#' ReferenceTrack class and methods
#'
#'
#' A class allow for on-demand streaming of data off the file system.
#'
#' The `availableDefaultMappings` function can be used to find out whether
#' the package defines a mapping scheme between one of the many supported input
#' file types and the metadata columns of the tracks' `GRanges` objects.
#'
#' @name ReferenceTrack-class
#'
#' @slot stream Object of class function. The import function to stream data
#' of the file system. Needs to be able to handle the two mandatory arguments
#' `file` (a `character` containing a valid file path) and `selection`
#' (a `GRanges` object with the genomic region to plot).
#' @slot reference Object of class "character", the path to the file
#' containing the data.
#' @slot mapping Object of class `list`, a default mapping between the
#' metadata columns of the returned `GRanges` object from the import function
#' and the `elemenMetadata` columns that make up the final track object.
#' @slot args Object of class `list`, the passed in constructor arguments
#' during object instantiation. Those will be needed when fetching the data
#' in order to fill all necessary slots.
#' @slot defaults Object of class `list`, the relevant default values to be
#' used when neither mapping nor args provides the necessary information.
#'
#' @return
#' Constructor functions of `AnnotationTrack`, `DataTrack`, `SequenceTrack`
#' and `AlignmentsTrack`` can create a special object of corresponding
#' `Reference*Track` subclass with pointer to the referenced
#' file.
#'
#' @return A virtual class: No objects may be created from it.
#'
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @examples
#' # This is a reference class, below example from AlignmentsTrack
#'
#' afrom <- 2960000
#' ato <- 3160000
#' alTrack <- AlignmentsTrack(system.file(
#'     package = "Gviz", "extdata",
#'     "gapped.bam"
#' ), isPaired = TRUE)
#' plotTracks(alTrack, from = afrom, to = ato, chromosome = "chr12")
#' @exportClass ReferenceTrack
setClass("ReferenceTrack",
    representation = representation("VIRTUAL",
        stream = "function",
        reference = "character",
        mapping = "list",
        args = "list",
        defaults = "list"
    ),
    prototype = prototype(
        stream = function(x, selection) {},
        reference = "~",
        mapping = list()
    ),
    validity = function(object) {
        msg <- NULL
        if (!all(c("file", "selection") %in% names(formals(object@stream)))) {
            msg <- "The streaming function in the 'stream' slot needs to define two arguments, 'file' and 'selection'"
        }
        if (!file.exists(object@reference)) {
            msg <- c(msg, sprintf("The referenced file '%s' does not exist", object@reference))
        }
        return(if (is.null(msg)) TRUE else msg)
    }
)

## Initialize ----------------------------------------------------------------

#' @describeIn ReferenceTrack-class Initialize.
#' @export

setMethod("initialize", "ReferenceTrack", function(.Object, stream, reference, mapping = list(),
                                                   args = list(), defaults = list()) {
    .Object@stream <- stream
    .Object@reference <- reference
    .Object@mapping <- mapping
    .Object@args <- args
    .Object@defaults <- defaults
    validObject(.Object)
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

## For character scalars the data need to be extracted from a file and we have to deal with parser functions
## and column assignments here.
setMethod(
    ".buildRange", signature("character"),
    function(range, importFun = NULL, trackType, stream = FALSE, args, defaults, autodetect = is.null(importFun), ...) {
        .checkClass(range, "character", 1)
        .checkClass(importFun, c("NULL", "function"), mandatory = FALSE)
        .checkClass(stream, "logical", 1)
        ## We first check for the default column mapping and whether this is a streaming file
        defMap <- .defaultVarMap(tolower(.fileExtension(range)), trackType, stream, !autodetect)
        isStream <- !is.null(defMap[[".stream"]]) && defMap[[".stream"]]
        defMap[[".stream"]] <- NULL
        if (!isStream) {
            data <- if (is.null(importFun)) {
                .registerImportFun(range)
            } else {
                if (!"file" %in% names(formals(importFun))) {
                    stop("The user-defined import function needs to define a 'file' argument")
                }
                importFun(range)
            }
            if (!is(data, "GRanges")) {
                stop(
                    "The import function did not provide a valid GRanges object. Unable to build track from file '",
                    range, "'"
                )
            }
            if (trackType == "DataTrack") {
                ## For data tracks we take all numeric data columns regardless of any mapping
                mc <- .prepareDtData(as.data.frame(mcols(data)), length(data))
                mcols(data) <- t(mc)
            } else {
                ## For the rest we use the mapping as provided by the constructor if available
                cmap <- .resolveColMapping(data, args, defMap)
                args <- cmap$args
                data <- cmap$data
            }
            args[["chromosome"]] <- as.character(seqnames(data))
            args[["strand"]] <- as.character(strand(data))
            return(.buildRange(range = data, args = args, defaults = defaults, trackType = trackType, ...))
        } else {
            if (trackType != "DataTrack") {
                for (i in names(defMap)) {
                    if (is.character(args[[i]]) && length(args[[i]]) == 1) {
                        defMap[[i]] <- args[[i]]
                    }
                }
            }
            return(list(
                reference = path.expand(range), mapping = defMap,
                stream = if (is.null(importFun)) .registerImportFun(range) else importFun
            ))
        }
    }
)

## Show ----------------------------------------------------------------------
## A helper function to plot general information about a ReferenceTrack
.referenceTrackInfo <- function(object, type) {
    message(sprintf(
        "%s '%s'\n| genome: %s\n| active chromosome: %s\n| referenced file: %s",
        type,
        names(object),
        genome(object),
        chromosome(object),
        object@reference
    ))
    if (length(object@mapping) && type != "ReferenceDataTrack") {
        message(sprintf("| mapping: %s\n", paste(names(object@mapping), as.character(object@mapping), sep = "=", collapse = ", ")))
    }
}
