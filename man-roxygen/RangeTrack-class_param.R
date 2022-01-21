#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param GdObject Object of `GdObject-class`.
#' @param value Value to be set.
#' @param from,to Numeric scalar, giving the range of genomic coordinates to
#' limit the tracks in. Note that `to` cannot be larger than `from.`
#' @param i Numeric scalar, index to subset.
#' @param j Numeric scalar, index to subset. Ignored.
#' @param f `factor` in the sense that `as.factor(f)` defines the grouping,
#' @param sort `logical`.
#' @param drop `logical`, indicating if levels that do not occur should be dropped (if `f` is a factor).
#' @param use.defaults `logical`.
#' @param ... Additional arguments.
#' @param .Object .Object
#' @param range range
#' @param chromosome chromosome
#' @param genome genome
