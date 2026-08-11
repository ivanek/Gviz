#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param GdObject Object of `GdObject-class`.
#' @param value Value to be set.
#' @param from,to Numeric scalar, giving the range of genomic coordinates to
#' limit the tracks in. Note that `to` cannot be larger than `from.`
#' @param stacks `logical`. Set if stacking should  be preserved.
#' @param i Numeric scalar, index to subset.
#' @param j Numeric scalar, index to subset. Ignored.
#' @param sort `logical`. Sort the track's ranges by their genomic coordinates
#' after subsetting.
#' @param drop `logical`, indicating if levels that do not occur should be
#' dropped (if `f` is a factor).
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#' @param stacking The stacking type for overlapping items of the track. One in
#' `c(hide, dense, squish, pack, full)`. Currently, only squish (make best use
#' of the available space), dense (no stacking, collapse overlapping ranges),
#' and hide (do not show any track items at all) are implemented.
#' @param ... Additional arguments.
