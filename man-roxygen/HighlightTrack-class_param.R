#' @param GdObject Object of class [`GdObject`][GdObject-class].
#' @param name Name of the retrieved parameter.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param value Value to be set.
#' @param recursive `logical`. For composite tracks, also set the display
#' parameters on each of the contained sub-tracks.
#' @param object Object of class [`HighlightTrack`][HighlightTrack-class].
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#'
#' @param trackList A list of `Gviz` track objects that all have to inherit
#' from class [`GdObject`][GdObject-class].
#' @param range An optional meta argument to handle the different input types.
#' If the `range` argument is missing, all the relevant information to create
#' the object has to be provided as individual function arguments (see below).
#'
#' The different input options for `range` are:
#'
#' * A [`GRanges`][GenomicRanges::GRanges-class] object: the genomic ranges for
#' the highlighting regions.
#' * An [`IRanges`][IRanges::IRanges-class] object: almost identical to the
#' [`GRanges`][GenomicRanges::GRanges-class] case, except that the chromosome
#' information has to be provided in the separate `chromosome` argument,
#' because it cannot be directly encoded in an
#' [`IRanges`][IRanges::IRanges-class] object.
#' * A `data.frame` object: the `data.frame` needs to contain at least the two
#' mandatory columns `start` and `end` with the range coordinates. It may also
#' contain a `chromosome` column with the chromosome information for each range.
#' If missing, this information will be drawn from the constructor's
#' `chromosome` argument.
#'
#' @param start,end An integer scalar with the genomic start or end coordinate
#' for the highlighting range. Can also be supplied as part of the `range`
#' argument.
#' @param width An integer vector of widths for the highlighting ranges. This
#' can be used instead of either `start` or `end` to specify the range
#' coordinates.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' `options(ucscChromosomeNames=TRUE)`. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form `chrx`, where `x` may
#' be any possible string. The user has to make sure that the respective
#' chromosome is indeed defined for the track's genome. If not provided here,
#' the constructor will try to build the chromosome information based on the
#' available inputs, and as a last resort will fall back to the value `chrNA`.
#' Please note that by definition all objects in the `Gviz` package can only
#' have a single active chromosome at a time (although internally the
#' information for more than one chromosome may be present), and the user has to
#' call the `chromosome<-` replacement method in order to change to a different
#' active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If not provided here, the constructor will try to
#' extract this information from the provided inputs, and eventually will fall
#' back to the default value of `NA`.
#' @param name Character scalar of the track's name. This is not really used
#' and only exists for completeness.
#' @param ... All additional parameters are ignored.
