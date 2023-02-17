#' @param GdObject Object of `GdObject-class`.
#' @param name Name of the retrieved parameter.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param value Value to be set.
#' @param recursive `logical`
#' @param object object
#' @param .Object .Object
#'
#' @param trackList A list of Gviz track objects that all have to inherit from
#' class \code{GdObject}.
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
#' \item{A \code{GRanges} object:}{ the genomic ranges for the highlighting
#' regions.}
#'
#' \item{An \code{\linkS4class{IRanges}} object:}{ almost identical to the
#' \code{GRanges} case, except that the chromosome information has to be
#' provided in the separate \code{chromosome} argument, because it can not be
#' directly encoded in an \code{IRanges} object.}
#'
#' \item{A \code{data.frame} object:}{ the \code{data.frame} needs to contain
#' at least the two mandatory columns \code{start} and \code{end} with the
#' range coordinates. It may also contain a \code{chromosome} column with the
#' chromosome information for each range. If missing, this information will be
#' drawn from the constructor's \code{chromosome} argument.}
#'
#' }
#' @param start,end An integer scalar with the genomic start or end coordinates
#' for the highlighting range. Can also be supplied as part of the \code{range}
#' argument.
#' @param width An integer vector of widths for highlighting ranges. This can
#' be used instead of either \code{start} or \code{end} to specify the range
#' coordinates.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' \code{options(ucscChromosomeNames=TRUE)}. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form \code{chrx}, where
#' \code{x} may be any possible string. The user has to make sure that the
#' respective chromosome is indeed defined for the the track's genome. If not
#' provided here, the constructor will try to build the chromosome information
#' based on the available inputs, and as a last resort will fall back to the
#' value \code{chrNA}. Please note that by definition all objects in the
#' \code{Gviz} package can only have a single active chromosome at a time
#' (although internally the information for more than one chromosome may be
#' present), and the user has to call the \code{chromosome<-} replacement
#' method in order to change to a different active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If not provided here the constructor will try to
#' extract this information from the provided inputs, and eventually will fall
#' back to the default value of \code{NA}.
#' @param name Character scalar of the track's name. This is not really used
#' and only exists fro completeness.
#' @param \dots All additional parameters are ignored.
