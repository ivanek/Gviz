#' @slot coords Object of class `matrix`, the image map coordinates, in the
#' order bottom-left x, bottom-left y, top-right x, top-right y. Row names are
#' mandatory for the matrix and have to be unique.
#' @slot tags Object of class `list`, the individual HTML tags for the image
#' map. The value of each list item has to be a named character vector, where
#' the names must match back into the row names of the `coords` matrix.
