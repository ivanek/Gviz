#' @param GdObject Object of `GdObject-class`.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param minBase Start of the sequence.
#' @param maxBase End of the sequence.
#' @param reference Name of the file (for streatming).
#' @param value Value to be set.
#' @param prepare `logical`
#' @param object object
#' @param .Object .Object

#' @param chromosome The currently active chromosome of the track. A valid UCSC
#' chromosome identifier if \code{options(ucscChromosomeNames=TRUE)}. Please
#' note that in this case only syntactic checking takes place, i.e., the
#' argument value needs to be an integer, numeric character or a character of
#' the form \code{chrx}, where \code{x} may be any possible string. The user
#' has to make sure that sequences for the respective chromosomes are indeed
#' part of the object. If not provided here, the constructor will set it to the
#' first available sequence. Please note that by definition all objects in the
#' \code{Gviz} package can only have a single active chromosome at a time
#' (although internally the information for more than one chromosome may be
#' present), and the user has to call the \code{chromosome<-} replacement
#' method in order to change to a different active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. For a \code{SequenceBSgenomeTrack} object, the genome
#' information is extracted from the input \code{BSgenome} package. For a
#' \code{DNAStringSet} it has too be provided or the constructor will fall back
#' to the default value of \code{NA}.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param importFunction A user-defined function to be used to import the
#' sequence data from a file. This only applies when the \code{sequence}
#' argument is a character string with the path to the input data file. The
#' function needs to accept an argument \code{file} containing the file path
#' and has to return a proper \code{DNAStringSet} object with the sequence
#' information per chromosome. A set of default import functions is already
#' implemented in the package for a number of different file types, and one of
#' these defaults will be picked automatically based on the extension of the
#' input file name. If the extension can not be mapped to any of the existing
#' import function, an error is raised asking for a user-defined import
#' function. Currently the following file types can be imported with the
#' default functions: \code{fa/fasta} and \code{2bit}.
#'
#' Both file types support indexing by genomic coordinates, and it makes sense
#' to only load the part of the file that is needed for plotting. To this end,
#' the \code{Gviz} package defines the derived \code{ReferenceSequenceTrack}
#' class, which supports streaming data from the file system. The user
#' typically does not have to deal with this distinction but may rely on the
#' constructor function to make the right choice as long as the default import
#' functions are used. However, once a user-defined import function has been
#' provided and if this function adds support for indexed files, you will have
#' to make the constructor aware of this fact by setting the \code{stream}
#' argument to \code{TRUE}. Please note that in this case the import function
#' needs to accept a second mandatory argument \code{selection} which is a
#' \code{GRanges} object containing the dimensions of the plotted genomic
#' range. As before, the function has to return an appropriate
#' \code{DNAStringSet} object.
#' @param stream A logical flag indicating that the user-provided import
#' function can deal with indexed files and knows how to process the additional
#' \code{selection} argument when accessing the data on disk. This causes the
#' constructor to return a \code{ReferenceSequenceTrack} object which will grab
#' the necessary data on the fly during each plotting operation.
#' @param \dots Additional items which will all be interpreted as further
#' display parameters. See \code{\link{settings}} and the "Display Parameters"
#' section below for details.
