#' @include ImageMap-class.R
NULL

### DisplayPars-class --------------------------------------------------------

## The initial idea was for the class to uses pass by reference semantic,
## allowing to update the object without reassigning to a symbol. However this
## turned out to be problematic whenever a `GdObject` was copied and to avoid
## confusion with the users about unwanted side-effects we decided to deprecate
## this feature.

#' DisplayPars: A class to control the plotting parameters for GdObjects
#'
#' @name DisplayPars-class
#'
#' @description All tracks within this package are highly customizable. The
#' `DisplayPars` class facilitates this and provides a unified API to the
#' customization parameters.
#'
#' The individual parameters in a `DisplayParameters` class are stored as
#' pointers in an environment. This has the upshot of not having to copy the
#' whole track object when changing parameters, and parameters can be updated
#' without the need to explicitly reassign the track to a symbol (i.e.,
#' updating of parameters happens in place). The downside is that upon copying
#' of track objects, the parameter environment needs to be re-instantiated.
#'
#' Objects can be created using the constructor function `DisplayPars`.
#'
#' The default display parameters for a track object class can be queried using
#' the `availableDisplayPars` function.
#'
#' @slot pars an environment or a list containing parameter key value pairs.
#'
#' @template DisplayParams_param
#'
#' @return The return value of the constructor function is a new object of class
#' `DisplayPars`.
#'
#' @author Florian Hahne
#' @examples
#' ## Construct object
#' dp <- DisplayPars(col = "red", lwd = 2, transformation = log2)
#' dp
#'
#' ## Query parameters
#' displayPars(dp)
#' displayPars(dp, "col")
#' getPar(dp, c("col", "transformation"))
#'
#' ## Modify parameters
#' displayPars(dp) <- list(lty = 1, fontsize = 3)
#' setPar(dp, "pch", 20)
#' dp
#' @export
.DisplayPars <- setClass("DisplayPars", representation(pars = "ListOrEnv"))

### DisplayPars Constructor --------------------------------------------------

## The initializer needs to create the environment, we can't take this
## directly from the class prototype since this would result in using
## the same environment all the time

# #' @describeIn DisplayPars-class Initialize.
# #' @export
# setMethod("initialize", "DisplayPars", function(.Object, ...) {
#     e <- new.env(hash = TRUE)
#     args <- list(...)
#     n <- names(args)
#     if ((is.null(n) & length(args) == 1) || any(n == "")) {
#         stop("All supplied arguments must be named.")
#     }
#     .Object@pars <- args
#     return(.Object)
# })

## Constructor, all supplied arguments are added to the environment.

#' @describeIn DisplayPars-class Constructor function.
#' @export
DisplayPars <- function(...) {
    e <- new.env(hash = TRUE)
    args <- list(...)
    n <- names(args)
    if ((is.null(n) & length(args) == 1) || any(n == "")) {
        stop("All supplied arguments must be named.")
    }
    .DisplayPars(pars = args)
}

### DisplayPars Methods ------------------------------------------------------

## Update function and deprecation message to show that an old environment-based
## DisplayParameter object has been updated to a list-based object.
.updateDp <- function(x, interactive = TRUE) {
    if (interactive) {
        message(
            "Note that the behaviour of the 'setPar' method has changed. ",
            "You need to reassign the result to an object for the side effects ",
            "to happen. Pass-by-reference semantic is no longer supported."
        )
    }
    if (is.environment(x@pars)) {
        x@pars <- as.list(x@pars)
        message("The DisplayPars object has been updated to a list-based representation.")
    }
    return(x)
}

## The accessor methods for the DisplayPars class. (they need to be in
## here because they are being called in the prototypes). For the
## setter method, the input can either be a named list, or a single
## keyword/value pair. The DisplayPars class was first implemented as
## an environment, so we essentially had pass by reference semantic here,
## and the object could be modified without the need to reassign to a
## symbol. However this somewhat broke the the R paradigm, and 'setPar'
## has been deprecated in favour of the more standard 'DisplayPars<-'
## replacement method. The getter methods either return the list
## of all parameters, or a subset of parameters if their names are
## provided as a character vector.  Please note that for convenience
## the result is unlisted if only a single parameter is queried.
## Since version 1.7.3 we have introduced an alias table for display
## parameters in order to harmonize parameter names without loosing
## backwards compatibility.

### DisplayPars Methods Getters ----------------------------------------------

#' @describeIn DisplayPars-class Generics for `getPar`.
setGeneric("getPar", def = function(x, name, ...) standardGeneric("getPar"))

#' @describeIn DisplayPars-class Alias for the `displayPars` method.
#' @export
setMethod(
    "getPar", c("DisplayPars", "character"),
    function(x, name, asIs = FALSE) {
        aliasRes <- .dpAliasReverseTable[name]
        new <- is.na(aliasRes)
        aliasRes[new] <- name[new]
        if (is.environment(x@pars)) {
            name <- intersect(aliasRes, base::ls(x@pars, all.names = TRUE))
            tmp <- mget(name, x@pars)
        } else {
            name <- intersect(aliasRes, names(x@pars))
            tmp <- x@pars[name]
        }
        if (!asIs) {
            tmp <- if (is.list(tmp) && length(tmp) == 1) tmp[[1]] else tmp
        }
        return(if (length(tmp)) tmp else NULL)
    }
)

#' @describeIn DisplayPars-class Alias for the `displayPars` method.
#' @export
setMethod("getPar", c("DisplayPars", "missing"), function(x, hideInternal = TRUE) {
    pars <- as.list(x@pars)
    if (hideInternal) {
        pars <- pars[!grepl("^\\.__", names(pars))]
    }
    return(pars)
})

#' @describeIn DisplayPars-class Generics for `displayPars`.
setGeneric("displayPars", function(x, name, ...) standardGeneric("displayPars"))

#' @describeIn DisplayPars-class Returns all available display parameters.
#' @export
setMethod("displayPars", c("DisplayPars", "missing"), function(x, hideInternal = TRUE) getPar(x, hideInternal = hideInternal))

#' @describeIn DisplayPars-class Returns the value of a subset of display parameters, as identified by `name`.
#' @export
setMethod("displayPars", c("DisplayPars", "character"), function(x, name) getPar(x, name))

#' @describeIn DisplayPars-class Converts `DisplayPars` to `list`.
#' @export
setMethod("as.list", "DisplayPars", function(x) as(x, "list"))

setAs("DisplayPars", "list", function(from, to) if (!is.null(from)) as.list(from@pars) else list())

### DisplayPars Methods Setters ----------------------------------------------

#' @describeIn DisplayPars-class Generics for `SetPar`.
#'
#' SetPar generic function
#'
#' @param x object to set the displayPar value on
#' @param value named value to be set
#'
setGeneric("setPar", function(x, value, ...) standardGeneric("setPar"))

#' @describeIn DisplayPars-class Sets display parameters by the values of the
#' named `list` in value. Note that display parameters in the `DisplayPars-class`
#' are pass-by-reference, so no re-assignment to the symbol `obj` is necessary.
setMethod(
    "setPar", signature("DisplayPars", "list"),
    function(x, value, interactive = TRUE) {
        x <- .updateDp(x, interactive)
        aliasRes <- .dpAliasReverseTable[names(value)]
        new <- is.na(aliasRes)
        aliasRes[new] <- names(value)[new]
        x@pars[aliasRes] <- value
        return(x)
    }
)

#' @describeIn DisplayPars-class set the single display parameter name to value.
#' Note that display parameters in the DisplayPars class are pass-by-reference,
#' so no re-assignment to the symbol `obj` is necessary.
#' @export
setMethod(
    "setPar", signature("DisplayPars", "character"),
    function(x, name, value, interactive = TRUE) {
        if (!(length(name) == 1 && is.null(value)) && length(name) != length(value)) {
            stop("'name' and 'value' must be of equal length")
        }
        x <- .updateDp(x, interactive)
        for (i in seq_along(name)) {
            aliasRes <- .dpAliasReverseTable[name[i]]
            if (is.na(aliasRes)) {
                aliasRes <- name[i]
            }
            x@pars[[aliasRes]] <- value[[i]]
        }
        return(x)
    }
)

#' @describeIn DisplayPars-class Generics for `displayPars<-`.
setGeneric("displayPars<-",
    signature = c("x", "value"),
    function(x, recursive = FALSE, value) standardGeneric("displayPars<-")
)

#' @describeIn DisplayPars-class  Replaces or adds display parameters as provided by the named `list` items.
#' @export
setReplaceMethod("displayPars", signature("DisplayPars", "list"), function(x, recursive = FALSE, value) {
    x <- setPar(x, value, interactive = FALSE)
    return(x)
})

### Display parameters lookup table ------------------------------------------

## An alias table to define display parameter synonyms. Essentially this is a
## simple named list, where the element names represent the preferred parameter
## name under which the actual value is stored, and the element content is a
## character vector of synonyms. Please note that this structure is also turned
## into a reverse lookup table for fast access.

.dpAliasTable <- list(
    "fontcolor.title" = "col.title",
    "rotation.title" = c("rot.title", "rotate.title"),
    "fontcolor.group" = "col.group",
    "fontcolor.item" = "col.item",
    "just.group" = c("labelJust", "labelJustification"),
    "transcriptAnnotation" = "geneSymbols",
    "fontcolor.item" = c("fontcolor.exon", "fontcolor.feature"),
    "fontsize.item" = c("fontsize.exon", "fontsize.feature"),
    "fontfamily.item" = c("fontfamily.exon", "fontfamily.feature"),
    "fontface.item" = c("fontface.exon", "fontface.feature"),
    "lineheight.item" = c("lineheight.exon", "lineheight.feature"),
    "alpha.item" = c("alpha.exon", "alpha.feature"),
    "cex.item" = c("cex.exon", "cex.feature"),
    "col.mate" = c("col.mate", "col.mates", "col.pair", "col.pairs"),
    "lty.mate" = c("lty.mate", "lty.mates", "lty.pair", "lty.pairs"),
    "lwd.mate" = c("lwd.mate", "lwd.mates", "lwd.pair", "lwd.pairs"),
    "alpha.mate" = c("alpha.mate", "alpha.mates", "alpha.pair", "alpha.pairs"),
    "col.gap" = "col.gaps",
    "lwd.gap" = "lwd.gaps",
    "lty.gap" = "lty.gaps",
    "alpha.gap" = "alpha.gaps"
)
.dpAliasReverseTable <- character()
.dpAliasReverseTable[names(.dpAliasTable)] <- names(.dpAliasTable)
#' @importFrom Biobase listLen
.dpAliasReverseTable[unlist(.dpAliasTable, use.names = FALSE)] <- rep(names(.dpAliasTable), Biobase::listLen(.dpAliasTable))
.dpAliasReverseTable <- .dpAliasReverseTable[order(names(.dpAliasReverseTable))]

## Show ----------------------------------------------------------------------

#' @describeIn DisplayPars-class Show method.
#' @export
setMethod("show", "DisplayPars", function(object) {
    cat("Display parameters:\n")
    for (i in base::ls(object@pars))
    {
        cat(i, " = ", sep = "")
        o <- try(as.character(object@pars[[i]]), silent = TRUE)
        if (is(o, "try-error")) print(object@pars[[i]]) else cat(o, "\n")
    }
})

### InferredDisplayPars-class ------------------------------------------------

#' InferredDisplayPars-class
#'
#' A class to allow for querying of available display parameters.
#' Essentially this is a normal list with a bit of a fancier show method.
#'
#' @slot name character. The name of the class.
#' @slot inheritance character. A vector indicating the inheritance structure.
#'
#' @return list
#' @noRd
#' @keywords internal
setClass("InferredDisplayPars", representation(name = "character", inheritance = "character"), contains = "list")

#' @describeIn DisplayPars-class  Return `InferredDisplayPars` as a `list`.
#' @export
setMethod("as.list", "InferredDisplayPars", function(x) as(x, "list"))

setAs(
    "InferredDisplayPars", "list",
    function(from, to) {
        ll <- from@.Data
        names(ll) <- names(from)
        ll
    }
)

## Show ----------------------------------------------------------------------

#' @describeIn DisplayPars-class  Show method.
#' @export
setMethod("show", "InferredDisplayPars", function(object) {
    cat("\nThe following display parameters are available for '", object@name, "' objects:\n",
        "(see ? ", object@name, " for details on their usage)\n\n",
        sep = ""
    )
    for (i in names(object)) {
        cat(i, ifelse(object@inheritance[i] == object@name, "",
            paste(" (inherited from class '", object@inheritance[i], "')", sep = "")
        ),
        ": ",
        sep = ""
        )
        if (is.null(object[[i]])) {
            cat("NULL\n")
        } else {
            o <- try(as.character(object[[i]]), silent = TRUE)
            if (is(o, "try-error")) print(object[[i]]) else cat(o, "\n")
        }
    }
})

### availableDisplayPars function --------------------------------------------

#' @param class Either character scalar or object. Supported classes are:
#' `GdObject`, `GenomeAxisTrack`, `RangeTrack`, `NumericTrack`, `DataTrack`,
#' `IdeogramTrack`, `StackedTrack`, `AnnotationTrack`, `DetailsAnnotationTrack`,
#' `GeneRegionTrack`, `BiomartGeneRegionTrack`, `AlignmentsTrack`,
#' `SequenceTrack`, `SequenceBSgenomeTrack`, `SequenceDNAStringSetTrack`,
#' `SequenceRNAStringSetTrack`
#'
#' @return `availableDisplayPars` returns a list of the default display
#' parameters.
#'
#' @examples
#' ## Default parameters
#' availableDisplayPars("GenomeAxisTrack")
#' @describeIn DisplayPars-class Get default display parameters.
#' @export
availableDisplayPars <- function(class) {
    if (!is.character(class)) {
        class <- class(class)
    }
    class <- match.arg(class, c(
        "GdObject", "GenomeAxisTrack", "RangeTrack", "NumericTrack", "DataTrack",
        "IdeogramTrack", "StackedTrack", "AnnotationTrack", "DetailsAnnotationTrack",
        "GeneRegionTrack", "BiomartGeneRegionTrack", "AlignmentsTrack",
        "SequenceTrack", "SequenceBSgenomeTrack", "SequenceDNAStringSetTrack",
        "SequenceRNAStringSetTrack"
    ))
    parents <- names(getClassDef(class)@contains)
    .makeParMapping()
    pars <- .parMappings[c(parents, class)]
    finalPars <- inherited <- list()
    for (p in names(pars)) {
        finalPars[names(pars[[p]])] <- pars[[p]]
        inherited[names(pars[[p]])] <- p
    }
    finalPars <- finalPars[order(names(finalPars))]
    inherited <- inherited[order(names(inherited))]
    return(new("InferredDisplayPars",
        name = class,
        inheritance = unlist(inherited), finalPars
    ))
}
