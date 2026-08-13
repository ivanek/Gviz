#' @param GdObject Object of class [`GdObject`][GdObject-class].
#' @param value Value to be set.
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.

#' @param range An optional meta argument to handle the different input types.
#' If the `range` argument is missing, all the relevant information to create
#' the object has to be provided as individual function arguments (see below).
#'
#' The different input options for `range` are:
#'
#' * A [`GRanges`][GenomicRanges::GRanges-class] object: essentially all the
#' necessary information to create a [`DataTrack`][DataTrack-class] can be
#' contained in a single [`GRanges`][GenomicRanges::GRanges-class] object. The
#' track's coordinates are taken from the `start`, `end` and `seqnames` slots,
#' the genome information from the `genome` slot, and the numeric data values
#' can be extracted from additional metadata columns (please note that
#' non-numeric columns are being ignored with a warning). As a matter of fact,
#' calling the constructor on a [`GRanges`][GenomicRanges::GRanges-class]
#' object without further arguments, e.g. `DataTrack(range=obj)`, is equivalent
#' to calling the coerce method `as(obj, "DataTrack")`. Alternatively, the
#' [`GRanges`][GenomicRanges::GRanges-class] object may only contain the
#' coordinate information, in which case the numeric data part is expected to be
#' present in the separate `data` argument, and the ranges have to match the
#' dimensions of the data matrix. If `data` is not `NULL`, this will always take
#' precedence over anything defined in the `range` argument. See below for
#' details.
#' * An [`IRanges`][IRanges::IRanges-class] object: this is very similar to the
#' above case, except that the numeric data part now always has to be provided
#' in the separate `data` argument. Also the chromosome information must be
#' provided in the `chromosome` argument, because neither of the two can be
#' directly encoded in an [`IRanges`][IRanges::IRanges-class] object.
#' * A `data.frame` object: the `data.frame` needs to contain at least the two
#' mandatory columns `start` and `end` with the range coordinates. It may also
#' contain a `chromosome` column with the chromosome information for each range.
#' If missing, it will be drawn from the separate `chromosome` argument. All
#' additional numeric columns will be interpreted as data columns, unless the
#' `data` argument is explicitly provided.
#' * A `character` scalar: in this case the value of the `range` argument is
#' considered to be a file path to an annotation file on disk. A range of file
#' types is supported by the `Gviz` package as identified by the file extension.
#' See the `importFunction` documentation below for further details.
#'
#' @param start,end,width Integer vectors, giving the start and the end
#' coordinates for the individual track items, or their width. Two of the three
#' need to be specified, and have to be of equal length or of length one, in
#' which case the single value will be recycled accordingly. Otherwise, the
#' usual R recycling rules for vectors do not apply and the function will throw
#' an error.
#' @param data A numeric matrix of data points with the number of columns equal
#' to the number of coordinates in `range`, or a numeric vector of appropriate
#' length that will be coerced into such a one-row matrix. Each individual row
#' is supposed to contain data for a given sample, where the coordinates for
#' each single observation are constant across samples. Depending on the
#' plotting type of the data (see the 'Details' and 'Display Parameters'
#' sections), sample grouping or data aggregation may be available.
#' Alternatively, this can be a character vector of column names that point
#' into the element metadata of the `range` object for subsetting. Naturally,
#' this is only supported when the `range` argument is of class
#' [`GRanges`][GenomicRanges::GRanges-class].
#' @param strand Character vector, the strand information for the individual
#' track items. Currently this has to be unique for the whole track and does not
#' really have any visible consequences, but we might decide to make
#' [`DataTrack`][DataTrack-class] objects strand-specific at a later stage.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' `options(ucscChromosomeNames=TRUE)`. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form `chrx`, where `x` may
#' be any possible string. The user has to make sure that the respective
#' chromosome is indeed defined for the track's genome. If not provided here,
#' the constructor will try to construct the chromosome information based on
#' the available inputs, and as a last resort will fall back to the value
#' `chrNA`. Please note that by definition all objects in the `Gviz` package can
#' only have a single active chromosome at a time (although internally the
#' information for more than one chromosome may be present), and the user has to
#' call the `chromosome<-` replacement method in order to change to a different
#' active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If not provided here, the constructor will try to
#' extract this information from the provided input, and eventually will fall
#' back to the default value of `NA`.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param importFunction A user-defined function to be used to import the data
#' from a file. This only applies when the `range` argument is a character
#' string with the path to the input data file. The function needs to accept an
#' argument `file` containing the file path and has to return a proper
#' [`GRanges`][GenomicRanges::GRanges-class] object with the data part attached
#' as numeric metadata columns. Essentially the process is equivalent to
#' constructing a [`DataTrack`][DataTrack-class] directly from a
#' [`GRanges`][GenomicRanges::GRanges-class] object in that non-numeric columns
#' will be dropped, and further subsetting can be achieved by means of the
#' `data` argument. A set of default import functions is already implemented in
#' the package for a number of different file types, and one of these defaults
#' will be picked automatically based on the extension of the input file name.
#' If the extension cannot be mapped to any of the existing import functions, an
#' error is raised asking for a user-defined import function. Currently the
#' following file types can be imported with the default functions: `wig`,
#' `bigWig/bw`, `bedGraph` and `bam`.
#'
#' Some file types support indexing by genomic coordinates (e.g., `bigWig` and
#' `bam`), and it makes sense to only load the part of the file that is needed
#' for plotting. To this end, the `Gviz` package defines the derived
#' [`ReferenceDataTrack`][ReferenceDataTrack-class] class, which supports
#' streaming data from the file system. The user typically does not have to deal
#' with this distinction but may rely on the constructor function to make the
#' right choice as long as the default import functions are used. However, once
#' a user-defined import function has been provided and if this function adds
#' support for indexed files, you will have to make the constructor aware of
#' this fact by setting the `stream` argument to `TRUE`. Please note that in
#' this case the import function needs to accept a second mandatory argument
#' `selection`, which is a [`GRanges`][GenomicRanges::GRanges-class] object
#' containing the dimensions of the plotted genomic range. As before, the
#' function has to return an appropriate
#' [`GRanges`][GenomicRanges::GRanges-class] object.
#' @param stream A logical flag indicating that the user-provided import
#' function can deal with indexed files and knows how to process the additional
#' `selection` argument when accessing the data on disk. This causes the
#' constructor to return a [`ReferenceDataTrack`][ReferenceDataTrack-class]
#' object which will grab the necessary data on the fly during each plotting
#' operation.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.
