#' @param GdObject Object of `GdObject-class`.
#' @param .Object The object skeleton passed on by `new()` during class
#' instantiation, to be filled in by the `initialize` method.
#' @param stackRanges A `GRanges` object with the merged extent of each read
#' group, i.e., the region spanned by a read and its mate. This is what the
#' reads are stacked on, and its names match the `groupid` metadata column.
#' @param stacks `logical`. Set if stacking should  be preserved.
#' @param sequences A `DNAStringSet` with the read sequences.
#' @param stream A function to stream the data from an indexed file. It has to
#' accept the two arguments `file` and `selection`, and to return a `GRanges`
#' object.
#' @param reference A `character` scalar with the path to the referenced file.
#' @param mapping A named `list` mapping the columns of the imported data to
#' the metadata columns of the track's `GRanges` object.
#' @param args A `list` of the arguments the object has been constructed with.
#' @param defaults A `list` of the constructor's default arguments, used to
#' fill in whatever has not been provided in `args`.
#' @param x A valid track object class name, or the object itself, in which
#' case the class is derived directly from it.
#' @param from,to Numeric scalar, giving the range of genomic coordinates to
#' limit the tracks in. Note that `to` cannot be larger than `from.`
#' @param use.defaults `logical`. Derive the subsetting range from the track's
#' own defaults rather than using `from` and `to` verbatim.
#' @param minBase,maxBase Numeric scalar, the start and end coordinates of the
#' plotting range.
#' @param prepare `logical`. Run the drawing method in preparation rather than
#' in production mode, i.e., only compute the track's layout without rendering
#' anything to the device.
#' @param subset `logical`. Subset the track to the current plotting range
#' before drawing.
#' @param object Object inheriting from `AlignmentsTrack-class`.
#' @param value Value to be set.


#' @param start,end,width Integer vectors, giving the start and the end
#' coordinates for the individual track items, or their width. Two of the three
#' need to be specified, and have to be of equal length or of length one, in
#' which case this single value will be recycled. Otherwise, the usual R
#' recycling rules for vectors do not apply here.
#' @param id Character vector of read identifiers. Those identifiers have to be
#' unique, i.e., each range representing a read needs to have a unique
#' `id`.
#' @param cigar A character vector of valid CIGAR strings describing details of
#' the alignment. Typically those include alignment gaps or insertions and
#' deletions, but also hard and soft clipped read regions. If missing, a fully
#' mapped read without gaps or indels is assumed. Needs to be of equal length
#' as the provided genomic coordinates, or of length 1.
#' @param mapq A numeric vector of read mapping qualities. Needs to be of equal
#' length as the provided genomic coordinates, or of length 1.
#' @param flag A named integer vector of length 2, as is produced by
#' Rsamtools::scanBamFlag(), used to filter out undesirable reads. If missing,
#' all mapped reads will be included.
#' @param isize A numeric vector of empirical insert sizes. This only applies
#' if the reads are paired. Needs to be of equal length as the provided genomic
#' coordinates, or of length 1. Currently not used.
#' @param groupid A factor (or vector than can be coerced into one) defining
#' the read pairs. Reads with the same `groupid` are considered to be
#' mates. Please note that each read group may only have one or two members.
#' Needs to be of equal length as the provided genomic coordinates, or of
#' length 1.
#' @param status A factor describing the mapping status of a read. Has to be
#' one in `mated`, `unmated` or `ambiguous`. Needs to be of
#' equal length as the provided genomic coordinates, or of length 1.
#' @param md A character vector describing the mapping details. This is
#' effectively and alternative to the CIGAR encoding and it removes the
#' dependency on a reference sequence to figure out read mismatches. Needs to
#' be of equal length as the provided genomic coordinates, or of length 1.
#' Currently not used.
#' @param seqs `DNAStringSet` of read sequences.
#' @param strand Character vector, the strand information for the reads. It may
#' be provided in the form `+` for the Watson strand, `-` for the
#' Crick strand or `*` for either one of the two. Needs to be of equal
#' length as the provided genomic coordinates, or of length 1. Please note that
#' paired reads need to be on opposite strands, and erroneous entries will
#' result in casting of an error.
#' @param chromosome The chromosome on which the track's genomic ranges are
#' defined. A valid UCSC chromosome identifier if
#' `options(ucscChromosomeNames=TRUE)`. Please note that in this case only
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
#' `c(hide, dense, squish, pack, full)`. Currently, only squish (make best
#' use of the available space), dense (no stacking, collapse overlapping
#' ranges), and hide (do not show any track items at all) are implemented.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param isPaired A logical scalar to determine whether the reads are paired
#' or not. While this may be used to render paired-end data as single-end, the
#' oppsite will typically not have any effect because the appropriate
#' `groupid` settings will not be present.  Thus setting `isPaired`
#' to `TRUE` can usually be used to autodetect the pairing state of the
#' input data.
#' @param importFunction A user-defined function to be used to import the data
#' from a file. This only applies when the `range` argument is a character
#' string with the path to the input data file. The function needs to accept an
#' argument `x` containing the file path and a second argument
#' `selection` with the desired plotting ranges. It has to return a proper
#' `GRanges` object with all the necessary metadata columns set. A single
#' default import function is already implemented in the package for `BAM`
#' files.
#' @param referenceSequence An optional [`SequenceTrack`][SequenceTrack-class]
#' object containing the reference sequence against which the reads have been
#' aligned. This is only needed when mismatch information has to be added to the
#' plot (i.e., the `showMismatchs` display parameter is `TRUE`) because this is
#' normally not encoded in the `BAM` file. If not provided through this
#' argument, the [`plotTracks`] function is smart enough to detect the presence
#' of a [`SequenceTrack`][SequenceTrack-class] object in the track list and will
#' use that as a reference sequence.
#' @param ... Additional items which will all be interpreted as further
#' display parameters. See [`settings`] and the "Display Parameters"
#' section below for details.
