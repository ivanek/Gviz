#' @include utils.R
NULL

### ImageMap-class -----------------------------------------------------------

#' ImageMap: HTML image map information
#'
#'
#' HTML image map information for annotation tracks.
#'
#' Objects of the `ImageMap-class` are usually not created by the user,
#' hence the constructor function `ImageMap` is not exported in the name space.
#'
#' @template ImageMap-class_param
#'
#' @template ImageMap-class_slot
#'
#' @name ImageMap-class
#'
#' @examples
#' ## Not provided. This is an internal structure.
#'
#' @export
setClass("ImageMap", representation(coords = "matrix", tags = "list"),
    prototype = prototype(coords = matrix(1, ncol = 4, nrow = 0), tags = list())
)

## Allow for NULL value slots
setClassUnion("ImageMapOrNULL", c("ImageMap", "NULL"))


### ImageMap Constructor -----------------------------------------------------

#' @noRd
#' @keywords internal
ImageMap <- function(coords, tags) {
    if (!(is.matrix(coords) && is.numeric(coords) && ncol(coords) == 4)) {
        stop("'coords' must be a numeric matrix with 4 columns")
    }
    rn <- rownames(coords)
    if (is.null(rn)) {
        stop("Rownames must be set for the matrix in 'coords'")
    }
    if (!is.list(tags) || is.null(names(tags)) || any(names(tags) == "") ||
        !all(vapply(tags, is.character, logical(1)))) {
        stop("'tags' must be a named list with character vector items.")
    }
    n <- unique(unlist(lapply(tags, names)))
    if (is.null(n) || any(n == "")) {
        stop("All items in the 'tags' list must be named character vectors.")
    }
    m <- n %in% rn
    if (!all(m)) {
        stop(
            "The following values in the 'tags' list ",
            "could not be mapped to the 'coords' matrix:\n",
            paste(n[!m], sep = "", collapse = ", ")
        )
    }
    new("ImageMap", coords = coords, tags = tags)
}

####  ImageMap Methods coords  and tags --------------------------------------

#' @describeIn ImageMap-class Generics for `coords`.
#' @exportMethod coords
#' @keywords internal
setGeneric("coords", function(ImageMap, ...) standardGeneric("coords"))

#' @describeIn ImageMap-class Returns the coordinates from the image map.
#' @export
setMethod("coords", "NULL", function(ImageMap) NULL)

#' @describeIn ImageMap-class Returns the coordinates from the image map.
#' @return Returns the coordinates from the image map.
#' @export
setMethod("coords", "ImageMap", function(ImageMap) ImageMap@coords)


#' @describeIn ImageMap-class Generics for `tags`.
#' @exportMethod tags
#' @keywords internal
setGeneric("tags", function(ImageMap, ...) standardGeneric("tags"))

#' @describeIn ImageMap-class Returns the tags from the image map
#' @export
setMethod("tags", "NULL", function(ImageMap) NULL)

#' @describeIn ImageMap-class Returns the tags from the image map
#' @return Returns the tags from the image map.
#' @export
setMethod("tags", "ImageMap", function(ImageMap) ImageMap@tags)
