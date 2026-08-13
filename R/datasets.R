#' Data sets
#'
#' Some sample data sets used for the illustrative examples and the vignette.
#'
#' @name datasets
#'
#' @aliases bmTrack cyp2b10 idTrack biomTrack biomTrack2 cpgIslands axTrack
#' @aliases conservation ensGenes denseAnnTrack geneModels iTrack itrack
#' @aliases idxTrack ideoTrack twoGroups from gcContent knownGenes refGenes
#' @aliases snpLocations to ctrack geneDetails dtHoriz bmt
#'
#' @format A heterogeneous collection of pre-built objects, one per data
#' set, used to avoid recomputing or re-downloading data in examples and
#' the vignette. Depending on the data set, the object is either a
#' ready-to-plot `Gviz` track (e.g. [`AnnotationTrack`][AnnotationTrack-class],
#' [`GeneRegionTrack`][GeneRegionTrack-class],
#' [`GenomeAxisTrack`][GenomeAxisTrack-class],
#' [`IdeogramTrack`][IdeogramTrack-class], [`DataTrack`][DataTrack-class] or
#' [`BiomartGeneRegionTrack`][BiomartGeneRegionTrack-class]), a
#' [`GRanges`][GenomicRanges::GRanges-class] or `data.frame` of annotation data,
#' or a numeric scalar (`from`, `to`) marking a genomic start or end coordinate.
#'
#' @docType data
#' @keywords datasets
NULL
