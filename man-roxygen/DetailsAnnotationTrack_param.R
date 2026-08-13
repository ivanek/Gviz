#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#' @param GdObject Object of class [`GdObject`][GdObject-class].
#' @param minBase,maxBase Numeric scalar, the start and end coordinates of the
#' plotting range.
#' @param prepare `logical`. Run the drawing method in preparation rather than
#' in production mode, i.e., only compute the track's layout without rendering
#' anything to the device.

#' @param fun A function that is being called for each entry in the
#' [`AnnotationTrack`][AnnotationTrack-class] object. See sections 'Details'
#' and 'Examples' for further information. When called internally by the
#' plotting machinery, a number of arguments are automatically passed on to this
#' function, and the user needs to make sure that they can all be digested
#' (i.e., either have all of them as formal named function arguments, or gobble
#' up everything that is not needed in `...`). These arguments are:
#'
#' * `start`: the genomic start coordinate of the range item.
#' * `end`: the genomic end coordinate of the range item.
#' * `strand`: the strand information for the range item.
#' * `chromosome`: the chromosome of the range item.
#' * `identifier`: the identifier of the range item, i.e., the result
#' of calling `identifier(DetailsAnnotationTrack, lowest=TRUE)`. Typically
#' those identifiers are passed on to the object constructor during
#' instantiation as the `id` argument.
#' * `index`: a counter enumerating the ranges. The
#' [`AnnotationTrack`][AnnotationTrack-class] object is sorted internally for
#' visibility, and the `index` argument refers to the index of plotting.
#' * `GdObject`: a reference to the currently plotted
#' [`DetailsAnnotationTrack`][DetailsAnnotationTrack-class] object.
#' * `GdObject.original`: a reference to the
#' [`DetailsAnnotationTrack`][DetailsAnnotationTrack-class] before any
#' processing like item collapsing has taken place. Essentially, this is the
#' track object as it exists in your working environment.
#'
#' Additional arguments can be passed to the plotting function by means of the
#' `detailsFunArgs` argument (see below). Note that the plot must use grid
#' graphics (e.g. functions in the `lattice` package or low-level grid
#' functions). To access a data object such as a matrix or data frame within the
#' function you can either store it as a variable in the global environment or,
#' to avoid name space conflicts, you can make it part of the function
#' environment by means of a closure. Alternatively, you may want to explicitly
#' stick it into an environment or pass it along in the `detailsFunArgs` list.
#' To figure out in your custom plotting function which annotation element is
#' currently being plotted you can either use the identifier, which has to be
#' unique for each range element, or you may want to use the genomic position
#' (start/end/strand/chromosome), e.g. if the data is stored in a
#' [`GRanges`][GenomicRanges::GRanges-class] object.
#' @param selectFun A function that is being called for each entry in the
#' [`AnnotationTrack`][AnnotationTrack-class] object with exactly the same
#' arguments as in `fun`. The purpose of this function is to decide for each
#' track element whether details should be drawn, and consequently it has to
#' return a single logical scalar. If the return value is `TRUE`, details will
#' be drawn for the item; if it is `FALSE`, the details strip for the item is
#' omitted.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.
