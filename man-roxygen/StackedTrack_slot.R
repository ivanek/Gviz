#' @slot stacking Object of class `character`, the stacking type of overlapping
#' items on the final plot. One in `c(hide, dense, squish, pack,full)`.
#' Currently, only `hide` (do not show the track items at all), `squish` (make
#' best use of the available space) and `dense` (no stacking at all) are
#' implemented.
#' @slot stacks Object of class `numeric`, holding the stack indices for each
#'  track item. This slot is usually populated by calling the `setStacks` method
#'  upon plotting, since the correct stacking is a function of the available
#'  plotting space.
