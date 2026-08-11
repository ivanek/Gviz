#' @param GdObject Object of `GdObject-class`.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param minBase Start of the sequence.
#' @param maxBase End of the sequence.
#' @param value Value to be set.
#' @param prepare `logical`. Run the drawing method in preparation rather than
#' in production mode, i.e., only compute the track's layout without rendering
#' anything to the device.
#' @param object Object inheriting from `SequenceTrack-class`.
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.

#' @param chromosome The currently active chromosome of the track. A valid UCSC
#' chromosome identifier if `options(ucscChromosomeNames=TRUE)`. Please
#' note that in this case only syntactic checking takes place, i.e., the
#' argument value needs to be an integer, numeric character or a character of
#' the form `chrx`, where `x` may be any possible string. The user
#' has to make sure that sequences for the respective chromosomes are indeed
#' part of the object. If not provided here, the constructor will set it to the
#' first available sequence. Please note that by definition all objects in the
#' `Gviz` package can only have a single active chromosome at a time
#' (although internally the information for more than one chromosome may be
#' present), and the user has to call the `chromosome<-` replacement
#' method in order to change to a different active chromosome.
#' @param genome The genome on which the track's ranges are defined. Usually
#' this is a valid UCSC genome identifier, however this is not being formally
#' checked at this point. For a `SequenceBSgenomeTrack` object, the genome
#' information is extracted from the input `BSgenome` package. For a
#' `DNAStringSet` it has too be provided or the constructor will fall back
#' to the default value of `NA`.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param importFunction A user-defined function to be used to import the
#' sequence data from a file. This only applies when the `sequence`
#' argument is a character string with the path to the input data file. The
#' function needs to accept an argument `file` containing the file path
#' and has to return a proper `DNAStringSet` object with the sequence
#' information per chromosome. A set of default import functions is already
#' implemented in the package for a number of different file types, and one of
#' these defaults will be picked automatically based on the extension of the
#' input file name. If the extension can not be mapped to any of the existing
#' import function, an error is raised asking for a user-defined import
#' function. Currently the following file types can be imported with the
#' default functions: `fa/fasta` and `2bit`.
#' Both file types support indexing by genomic coordinates, and it makes sense
#' to only load the part of the file that is needed for plotting. To this end,
#' the `Gviz` package defines the derived `ReferenceSequenceTrack`
#' class, which supports streaming data from the file system. The user
#' typically does not have to deal with this distinction but may rely on the
#' constructor function to make the right choice as long as the default import
#' functions are used. However, once a user-defined import function has been
#' provided and if this function adds support for indexed files, you will have
#' to make the constructor aware of this fact by setting the `stream`
#' argument to `TRUE`. Please note that in this case the import function
#' needs to accept a second mandatory argument `selection` which is a
#' `GRanges` object containing the dimensions of the plotted genomic
#' range. As before, the function has to return an appropriate
#' `DNAStringSet` object.
#' @param stream A logical flag indicating that the user-provided import
#' function can deal with indexed files and knows how to process the additional
#' `selection` argument when accessing the data on disk. This causes the
#' constructor to return a `ReferenceSequenceTrack` object which will grab
#' the necessary data on the fly during each plotting operation.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.
