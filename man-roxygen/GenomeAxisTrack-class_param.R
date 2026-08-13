#' @param range Optional object of class
#' [`GRanges`][GenomicRanges::GRanges-class] or
#' [`IRanges`][IRanges::IRanges-class] containing regions to be highlighted on
#' the axis by colored boxes.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param id A character vector of the same length as `range` containing
#' identifiers for the ranges. If missing, the constructor will try to extract
#' the ids from `names(range)`.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.
