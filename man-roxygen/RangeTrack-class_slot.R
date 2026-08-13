#' @slot range Object of class [`GRanges`][GenomicRanges::GRanges-class], the
#' genomic ranges of the track items as well as additional annotation
#' information in its `elementMetadata` slot. Please note that the slot is
#' actually implemented as a class union between
#' [`GRanges`][GenomicRanges::GRanges-class] and
#' [`IRanges`][IRanges::IRanges-class] to increase efficiency, for instance for
#' [`DataTrack`][DataTrack-class] objects. This usually does not concern the
#' user.
#' @slot chromosome Object of class `character`, the chromosome on which the
#' track is defined. There can only be a single chromosome for one track. For
#' certain subclasses, the space of allowed chromosome names is limited (e.g.,
#' only those chromosomes that exist for a particular genome). Throughout the
#' package, chromosome names have to be entered either as a single integer
#' scalar or as a character scalar of the form `chrXYZ`, where `XYZ` may be an
#' arbitrary character string.
#' @slot genome Object of class `character`, the genome for which the track is
#' defined. For most sub-classes this has to be a valid UCSC genome identifier,
#' however this may not always be formally checked upon object instantiation.
