#' Dynamic content based on the available resolution
#'
#' When plotting features linearly along genomic coordinates one frequently
#' runs into the problem of too little resolution to adequately display all
#' details. Most genome browsers try to reasonably reduce the amount of detail
#' that is shown based on the current zoom level.
#'
#' Most track classes in this package define an internal `collapseTrack`
#' method which tries to adjust the plotted content to the available
#' resolution, aims at reducing over-plotting and prevents rendering issues,
#' e.g. when lines are too thin to be plotted. This feature can be toggled on
#' or off using the `collapse` display parameter (see
#' `settings` for details on setting these parameters).
#'
#' In the simplest case (for `AnnotationTrack` objects) this
#' involves expanding all shown features to a minimum pixel width and height
#' (using display parameters `min.width` and `min.heigh`) and
#' collapsing overlapping annotation items (as defined by the parameter
#' `min.distance` into one single item to prevent over-plotting.
#'
#' For objects of class `DataTrack`, the data values
#' underlying collapsed regions will be summarized based on the `summary`
#' display parameter. See the class' documentation for more details.
#'
#' @name collapsing
#'
#' @seealso
#' \code{\linkS4class{AnnotationTrack}}
#'
#' \code{\linkS4class{DataTrack}}
#'
#' \code{\link{settings}}
NULL
