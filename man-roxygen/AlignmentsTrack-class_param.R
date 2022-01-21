#' @param GdObject Object of `GdObject-class`.
#' @param name Name of the retrieved parameter.
#' @param .Object .Object
#' @param stackRanges stackRanges
#' @param stacks stacks
#' @param sequences sequences
#' @param stream stream
#' @param reference reference
#' @param mapping mapping
#' @param args args
#' @param defaults defaults
#' @param x x
#' @param from,to from,to
#' @param use.defaults `logical`
#' @param minBase,maxBase minBase,maxBase
#' @param prepare `logical`
#' @param subset `logical`
#' @param object object
#' @param value value


#' @param start,end,width Integer vectors, giving the start and the end
#' coordinates for the individual track items, or their width. Two of the three
#' need to be specified, and have to be of equal length or of length one, in
#' which case this single value will be recycled. Otherwise, the usual R
#' recycling rules for vectors do not apply here.
#' @param id Character vector of read identifiers. Those identifiers have to be
#' unique, i.e., each range representing a read needs to have a unique
#' \code{id}.
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
#' the read pairs. Reads with the same \code{groupid} are considered to be
#' mates. Please note that each read group may only have one or two members.
#' Needs to be of equal length as the provided genomic coordinates, or of
#' length 1.
#' @param status A factor describing the mapping status of a read. Has to be
#' one in \code{mated}, \code{unmated} or \code{ambiguous}. Needs to be of
#' equal length as the provided genomic coordinates, or of length 1.
#' @param md A character vector describing the mapping details. This is
#' effectively and alternative to the CIGAR encoding and it removes the
#' dependency on a reference sequence to figure out read mismatches. Needs to
#' be of equal length as the provided genomic coordinates, or of length 1.
#' Currently not used.
#' @param seqs \code{DNAStringSet} of read sequences.
#' @param strand Character vector, the strand information for the reads. It may
#' be provided in the form \code{+} for the Watson strand, \code{-} for the
#' Crick strand or \code{*} for either one of the two. Needs to be of equal
#' length as the provided genomic coordinates, or of length 1. Please note that
#' paired reads need to be on opposite strands, and erroneous entries will
#' result in casting of an error.
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
#' @param stacking The stacking type for overlapping items of the track. One in
#' \code{c(hide, dense, squish, pack, full)}. Currently, only squish (make best
#' use of the available space), dense (no stacking, collapse overlapping
#' ranges), and hide (do not show any track items at all) are implemented.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting.
#' @param isPaired A logical scalar to determine whether the reads are paired
#' or not. While this may be used to render paired-end data as single-end, the
#' oppsite will typically not have any effect because the appropriate
#' \code{groupid} settings will not be present.  Thus setting \code{isPaired}
#' to \code{TRUE} can usually be used to autodetect the pairing state of the
#' input data.
#' @param importFunction A user-defined function to be used to import the data
#' from a file. This only applies when the \code{range} argument is a character
#' string with the path to the input data file. The function needs to accept an
#' argument \code{x} containing the file path and a second argument
#' \code{selection} with the desired plotting ranges. It has to return a proper
#' \code{GRanges} object with all the necessary metadata columns set. A single
#' default import function is already implemented in the package for \code{BAM}
#' files.
#' @param referenceSequence An optional \code{\linkS4class{SequenceTrack}}
#' object containing the reference sequence against which the reads have been
#' aligned. This is only needed when mismatch information has to be added to
#' the plot (i.e., the \code{showMismatchs} display parameter is \code{TRUE})
#' because this is normally not encoded in the \code{BAM} file. If not provided
#' through this argument, the \code{\link{plotTracks}} function is smart enough
#' to detect the presence of a \code{\linkS4class{SequenceTrack}} object in the
#' track list and will use that as a reference sequence.
#' @param \dots Additional items which will all be interpreted as further
#' display parameters. See \code{\link{settings}} and the "Display Parameters"
#' section below for details.
