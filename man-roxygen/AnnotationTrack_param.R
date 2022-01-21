#' @param range An optional meta argument to handle the different input types.
#'  If the `range` argument is missing, all the relevant information to
#'  create the object has to be provided as individual function arguments
#'  (see below).
#'
#' The different input options for `range` are:
#'
#'  * A `GRanges` object: the genomic ranges for the
#' `Annotation` track as well as the optional additional metadata columns
#' `feature`, `group` and `id` (see description of the
#' individual function parameters below for details). Calling the constructor
#' on a `GRanges` object without further arguments, e.g.
#' `AnnotationTrack(range=obj)` is equivalent to calling the coerce method
#' `as(obj, "AnnotationTrack")`.
#' * A `GRangesList` object: this is very similar to the previous
#' case, except that the grouping information that is part of the list
#' structure is preserved in the `AnnotationTrack`. I.e., all the elements
#' within one list item receive the same group id. For consistency, there is
#' also a coercion method from `GRangesLists` `as(obj,"AnnotationTrack")`.
#' * An \code{\linkS4class{IRanges}} object: almost identical to the
#' `GRanges` case, except that the chromosome and strand information as
#' well as all additional metadata has to be provided in the separate
#' `chromosome`, `strand`, `feature`, `group` or `id`
#' arguments, because it can not be directly encoded in an `IRange`
#' object. Note that none of those inputs are mandatory, and if not provided
#' explicitly the more or less reasonable default values `chromosome=NA`
#' and `strand="*"` are used.
#' * A `data.frame` object: the `data.frame` needs to contain
#' at least the two mandatory columns `start` and `end` with the
#' range coordinates. It may also contain a `chromosome` and a
#' `strand` column with the chromosome and strand information for each
#' range. If missing it will be drawn from the separate `chromosome` or
#' `strand` arguments. In addition, the `feature`, `group` and
#' `id` data can be provided as additional columns. The above comments
#' about potential default values also apply here.
#' * A `character` scalar: in this case the value of the
#' `range` argument is considered to be a file path to an annotation file
#' on disk. A range of file types are supported by the `Gviz` package as
#' identified by the file extension. See the `importFunction`
#' documentation below for further details.
#'
#' @param start,end,width Integer vectors, giving the start and the end end
#' coordinates for the individual track items, or their width. Two of the three
#' need to be specified, and have to be of equal length or of length one, in
#' which case this single value will be recycled. Otherwise, the usual R
#' recycling rules for vectors do not apply here.
#' @param feature Factor (or other vector that can be coerced into one), giving
#' the feature types for the individual track items. When plotting the track to
#' the device, if a display parameter with the same name as the value of
#' `feature` is set, this will be used as the track item's fill colour. See
#' `grouping` for details.  Needs to be of equal length as the
#' provided genomic coordinates, or of length 1.
#' @param group Factor (or other vector that can be coerced into one), giving
#' the group memberships for the individual track items. When plotting to the
#' device, all items in the same group will be connected.  See
#' `grouping` for details. Needs to be of equal length as the
#' provided genomic coordinates, or of length 1.
#' @param id Character vector of track item identifiers. When plotting to the
#' device, it's value will be used as the identifier tag if the display
#' parameter `showFeatureId=TRUE`. Needs to be of equal length as the
#' provided genomic ranges, or of length 1.
#' @param strand Character vector, the strand information for the individual
#' track items. It may be provided in the form `+` for the Watson strand,
#' `-` for the Crick strand or `*` for either one of the two. Needs
#' to be of equal length as the provided genomic coordinates, or of length 1.
#' Please note that grouped items need to be on the same strand, and erroneous
#' entries will result in casting of an error.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' \code{options(ucscChromosomeNames=TRUE)}. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form `chrx`, where
#' `x` may be any possible string. The user has to make sure that the
#' respective chromosome is indeed defined for the the track's genome. If not
#' provided here, the constructor will try to construct the chromosome
#' information based on the available inputs, and as a last resort will fall
#' back to the value `chrNA`. Please note that by definition all objects
#' in the `Gviz` package can only have a single active chromosome at a
#' time (although internally the information for more than one chromosome may
#' be present), and the user has to call the `chromosome<-` replacement
#' method in order to change to a different active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If not provided here the constructor will try to
#' extract this information from the provided input, and eventually will fall
#' back to the default value of `NA`.
#' @param stacking The stacking type for overlapping items of the track. One in
#' `c(hide, dense, squish, pack,full)`. Currently, only squish (make best
#' use of the available space), dense (no stacking, collapse overlapping
#' ranges), and hide (do not show any track items at all) are implemented.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param fun A function that is being called for each entry in the
#' `AnnotationTrack` object. See section 'Details' and 'Examples' for
#' further information. When called internally by the plotting machinery, a
#' number of arguments are automatically passed on to this function, and the
#' user needs to make sure that they can all be digested (i.e., either have all
#' of them as formal named function arguments, or gobble up everything that is
#' not needed in `...`). These arguments are:
#'
#' * `start`: the genomic start coordinate of the range item.
#' * `end`: the genomic end coordinates of the range item.
#' * `strand`: the strand information for the range item.
#' * `chromosome`: the chromosome of the range item.
#' * `identifier`: the identifier of the range item, i.e., the result
#' of calling `identifier(DetailsAnnotationTrack, lowest=TRUE)`. Typically
#' those identifiers are passed on to the object constructor during
#' instantiation as the `id` argument.
#' * `index`: a counter enumerating the ranges. The `AnnotationTrack` object
#' is sorted internally for visibility, and the `index` argument refers to
#' the index of plotting.
#' * `GdObject`: a reference to the currently plotted `DetailsAnnotationTrack`
#' object.
#' * `GdObject.original`: a reference to the `DetailsAnnotationTrack` before
#' any processing like item collapsing has taken place. Essentially, this is
#' the track object as it exists in your working environment.
#'
#' Additional arguments can be passed to the plotting function by means of the
#' `detailsFunArgs` argument (see below). Note that the plot must use grid
#' graphics (e.g. function in the 'lattice' package or low-level grid
#' functions). To access a data object such a matrix or data frame within the
#' function you can either store it as a variable in the global environment or,
#' to avoid name space conflicts, you can make it part of the function
#' environment by means of a closure. Alternatively, you may want to
#' explicitely stick it into an environment or pass it along in the
#' `detailsFunArgs` list. To figure out in your custom plotting function
#' which annotation element is currently being plotted you can either use the
#' identifier which has to be unique for each range element, or you may want to
#' use the genomic position (start/end/strand/chromosome) e.g. if the data is
#' stored in a `GRanges` object.
#'
#' @param selectFun A function that is being called for each entry in the
#' `AnnotationTrack` object with exactly the same arguments as in
#' `fun`. The purpose of this function is to decide for each track element
#' whether details should be drawn, and consequently it has to return a single
#' logical scalar. If the return value is `TRUE`, details will be drawn
#' for the item, if it is `FALSE`, the details strip for the item is
#' omitted.
#' @param importFunction A user-defined function to be used to import the data
#' from a file. This only applies when the `range` argument is a character
#' string with the path to the input data file. The function needs to accept an
#' argument `x` containing the file path and has to return a proper
#' `GRanges` object with all the necessary metadata columns set. A set of
#' default import functions is already implemented in the package for a number
#' of different file types, and one of these defaults will be picked
#' automatically based on the extension of the input file name. If the
#' extension can not be mapped to any of the existing import function, an error
#' is raised asking for a user-defined import function via this argument.
#' Currently the following file types can be imported with the default
#' functions: `gff`, `gff1`, `gff2`, `gff3`, `bed`, `bam`.
#' @param stream A logical flag indicating that the user-provided import
#' function can deal with indexed files and knows how to process the additional
#' `selection` argument when accessing the data on disk. This causes the
#' constructor to return a `ReferenceAnnotationTrack` object which will
#' grab the necessary data on the fly during each plotting operation.
#' @param \dots Additional items which will all be interpreted as further
#' display parameters. See [`settings`](settings) and the "Display Parameters"
#' section below for details.
