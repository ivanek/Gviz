#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param GdObject Object of class [`GdObject`][GdObject-class].
#' @param value Value to be set.
#' @param from,to Numeric scalar, giving the range of genomic coordinates to
#' limit the tracks in. Note that `from` cannot be larger than `to`.
#' @param i Numeric scalar, index to subset.
#' @param j Numeric scalar, index to subset. Ignored.
#' @param f `factor` in the sense that `as.factor(f)` defines the grouping used
#' to split the track into a list of tracks.
#' @param sort `logical`. Sort the track's ranges by their genomic coordinates
#' after subsetting.
#' @param drop `logical`, indicating if levels that do not occur should be
#' dropped (if `f` is a factor).
#' @param use.defaults `logical`. Derive the subsetting range from the track's
#' own defaults, honouring the `min.width` and similar display parameters,
#' rather than using `from` and `to` verbatim.
#' @param ... Additional arguments.
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#' @param range A [`GRanges`][GenomicRanges::GRanges-class] object with the
#' genomic coordinates of the track items, and any additional annotation stored
#' in its metadata columns.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' `options(ucscChromosomeNames=TRUE)`. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form `chrx`, where `x` may
#' be any possible string. The user has to make sure that the respective
#' chromosome is indeed defined for the track's genome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point.
