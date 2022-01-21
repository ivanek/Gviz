#' @param GdObject Object of `GdObject-class`.
#' @param value Value to be set.
#' @param .Object .Object

#' @param range An optional meta argument to handle the different input types.
#' If the `range` argument is missing, all the relevant information to create
#' the object has to be provided as individual function arguments (see below).
#'
#' The different input options for \code{range} are:
#'
#' \describe{
#'
#' \item{}{A \code{GRanges} object: essentially all the necessary information
#' to create a \code{DataTrack} can be contained in a single \code{GRanges}
#' object. The track's coordinates are taken from the \code{start}, \code{end}
#' and \code{seqnames} slots, the genome information from the genome slot, and
#' the numeric data values can be extracted from additional metadata columns
#' columns (please note that non-numeric columns are being ignored with a
#' warning). As a matter of fact, calling the constructor on a \code{GRanges}
#' object without further arguments, e.g.  \code{DataTrack(range=obj)} is
#' equivalent to calling the coerce method \code{as(obj, "DataTrack")}.
#' Alternatively, the \code{GRanges} object may only contain the coordinate
#' information, in which case the numeric data part is expected to be present
#' in the separate \code{data} argument, and the ranges have to match the
#' dimensions of the data matrix. If \code{data} is not \code{NULL}, this will
#' always take precedence over anything defined in the \code{range} argument.
#' See below for details.}
#'
#' \item{}{An \code{\linkS4class{IRanges}} object: this is very similar to the
#' above case, except that the numeric data part now always has to be provided
#' in the separate \code{data} argument. Also the chromosome information must
#' be provided in the \code{chromosome} argument, because neither of the two
#' can be directly encoded in an \code{IRange} object.}
#'
#' \item{}{A \code{data.frame} object: the \code{data.frame} needs to contain
#' at least the two mandatory columns \code{start} and \code{end} with the
#' range coordinates. It may also contain a \code{chromosome} column with the
#' chromosome information for each range. If missing it will be drawn from the
#' separate \code{chromosome} argument. All additional numeric columns will be
#' interpreted as data columns, unless the \code{data} argument is explicitely
#' provided.}
#'
#' \item{}{A \code{character} scalar: in this case the value of the
#' \code{range} argument is considered to be a file path to an annotation file
#' on disk. A range of file types are supported by the \code{Gviz} package as
#' identified by the file extension. See the \code{importFunction}
#' documentation below for further details.}
#'
#' }
#' @param start,end,width Integer vectors, giving the start and the end end
#' coordinates for the individual track items, or their width. Two of the three
#' need to be specified, and have to be of equal length or of length one, in
#' which case the single value will be recycled accordingly. Otherwise, the
#' usual R recycling rules for vectors do not apply and the function will cast
#' an error.
#' @param data A numeric matrix of data points with the number of columns equal
#' to the number of coordinates in \code{range}, or a numeric vector of
#' appropriate length that will be coerced into such a one-row matrix. Each
#' individual row is supposed to contain data for a given sample, where the
#' coordinates for each single observation are constant across samples.
#' Depending on the plotting type of the data (see 'Details' and 'Display
#' Parameters' sections), sample grouping or data aggregation may be available.
#' Alternatively, this can be a character vector of column names that point
#' into the element metadata of the \code{range} object for subsetting.
#' Naturally, this is only supported when the \code{range} argument is of class
#' \code{GRanges}.
#' @param strand Character vector, the strand information for the individual
#' track items. Currently this has to be unique for the whole track and doesn't
#' really have any visible consequences, but we might decide to make
#' \code{DataTracks} strand-specific at a later stage.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' \code{options(ucscChromosomeNames=TRUE)}. Please note that in this case only
#' syntactic checking takes place, i.e., the argument value needs to be an
#' integer, numeric character or a character of the form \code{chrx}, where
#' \code{x} may be any possible string. The user has to make sure that the
#' respective chromosome is indeed defined for the the track's genome. If not
#' provided here, the constructor will try to construct the chromosome
#' information based on the available inputs, and as a last resort will fall
#' back to the value \code{chrNA}. Please note that by definition all objects
#' in the \code{Gviz} package can only have a single active chromosome at a
#' time (although internally the information for more than one chromosome may
#' be present), and the user has to call the \code{chromosome<-} replacement
#' method in order to change to a different active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. If not provided here the constructor will try to
#' extract this information from the provided input, and eventually will fall
#' back to the default value of \code{NA}.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param importFunction A user-defined function to be used to import the data
#' from a file. This only applies when the \code{range} argument is a character
#' string with the path to the input data file. The function needs to accept an
#' argument \code{file} containing the file path and has to return a proper
#' \code{GRanges} object with the data part attached as numeric metadata
#' columns. Essentially the process is equivalent to constructing a
#' \code{DataTrack} directly from a \code{GRanges} object in that non-numeric
#' columns will be dropped, and further subsetting can be archived by means of
#' the \code{data} argument. A set of default import functions is already
#' implemented in the package for a number of different file types, and one of
#' these defaults will be picked automatically based on the extension of the
#' input file name. If the extension can not be mapped to any of the existing
#' import function, an error is raised asking for a user-defined import
#' function. Currently the following file types can be imported with the
#' default functions: \code{wig}, \code{bigWig/bw}, \code{bedGraph} and
#' \code{bam}.
#'
#' Some file types support indexing by genomic coordinates (e.g., \code{bigWig}
#' and \code{bam}), and it makes sense to only load the part of the file that
#' is needed for plotting. To this end, the \code{Gviz} package defines the
#' derived \code{ReferenceDataTrack} class, which supports streaming data from
#' the file system. The user typically does not have to deal with this
#' distinction but may rely on the constructor function to make the right
#' choice as long as the default import functions are used. However, once a
#' user-defined import function has been provided and if this function adds
#' support for indexed files, you will have to make the constructor aware of
#' this fact by setting the \code{stream} argument to \code{TRUE}. Please note
#' that in this case the import function needs to accept a second mandatory
#' argument \code{selection} which is a \code{GRanges} object containing the
#' dimensions of the plotted genomic range. As before, the function has to
#' return an appropriate \code{GRanges} object.
#' @param stream A logical flag indicating that the user-provided import
#' function can deal with indexed files and knows how to process the additional
#' \code{selection} argument when accessing the data on disk. This causes the
#' constructor to return a \code{ReferenceDataTrack} object which will grab the
#' necessary data on the fly during each plotting operation.
#' @param \dots Additional items which will all be interpreted as further
#' display parameters.
