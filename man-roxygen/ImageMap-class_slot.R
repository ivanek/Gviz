#' @slot coords Object of class `matrix`, the image map coordinates. In the
#' order x bl, y bl, x tr, y tr. Row names are mandatory for the matrix
#' and have to be unique.
#' @slot tags Object of class `list`, the individual HTML tags for the image map.
#' The value of each list item has to be a named character vector,
#' where the names must match back into the row names of the `coords` matrix
