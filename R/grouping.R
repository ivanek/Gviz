#' Grouping of annotation features
#'
#' Many annotation tracks are actually composed of a number of grouped
#' sub-features, for instance exons in a gene model. This man page highlights
#' the use of grouping information to build informative annotation plots.
#'
#'
#' All track objects that inherit from class
#' [`AnnotationTrack`][AnnotationTrack-class] support the grouping feature. The
#' information is usually passed on to the constructor function (for
#' [`AnnotationTrack`][AnnotationTrack-class] via the `group` argument and for
#' [`GeneRegionTrack`][GeneRegionTrack-class] objects via the `exon` argument)
#' or automatically downloaded from an online annotation repository
#' ([`BiomartGeneRegionTrack`][BiomartGeneRegionTrack-class]). Group membership
#' is specified by a factor vector with as many items as there are annotation
#' items in the track (i.e., the value of `length(track)`). Upon plotting, the
#' grouped annotation features are displayed together and will not be separated
#' in the stacking of track items.
#'
#' @name grouping
#'
#' @return No return value, called for documentation purposes only.
#'
#' @author Florian Hahne
#' @seealso
#'
#' [AnnotationTrack-class]
#'
#' [BiomartGeneRegionTrack-class]
#'
#' [GeneRegionTrack-class]
NULL
