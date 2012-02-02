\name{grouping}
\Rdversion{1.1}
\alias{grouping}

\title{Grouping of annotation features}

\description{

  Many annotation tracks are actually composed of a number of grouped
  sub-features, for instance exons in a gene model. This man page
  highlights the use of grouping information to build informative
  annotation plots.

}

\details{
  
  All track objects that inherit from class
  \code{\linkS4class{AnnotationTrack}} support the grouping feature. The
  information is usually passed on to the constructor function (for
  \code{AnnotationTrack} via the \code{groups} argument and for
  \code{\linkS4class{GeneRegionTrack}} objects via the \code{exon}
  argument) or automatically downloaded from an online annotation
  repository (\code{\linkS4class{BiomartGeneRegionTrack}}). Group
  membership is specified by a factor vector with as many items as there
  are annotation items in the track (i.e., the value of
  \code{length(track)}. Upon plotting, the grouped annotation features
  are displayed together and will not be separated in the stacking of
  track items.
  
}

\author{Florian Hahne}
\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{BiomartGeneRegionTrack}}

  \code{\linkS4class{GeneRegionTrack}}



}




