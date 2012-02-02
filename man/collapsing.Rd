\name{collapsing}
\Rdversion{1.1}
\alias{collapsing}

\title{Dynamic content based on the available resolution}

\description{

  When plotting features linearily along genomic coordinates
  one frequently runs into the problem of too little resolution to
  adequatelty display all details. Most genome browsers try to
  reasonably reduce the amount of detail that is shown based on the
  current zoomn level.
  
}

\details{
  
  Most track classes in this package define an internal
  \code{collapseTrack} method which tries to adjust the plotted content
  to the available resolution, aims at reducing overplotting and
  prevents rendering issues, e.g. when lines are too thin to be
  plotted. This feature can be toggled on or off using the
  \code{collapse} display parameter (see \code{\link{settings}} for
  details on setting these parameters).

  In the simplest case (for
  \code{\linkS4class{AnnotationTrack}} objects) this involves expanding
  all shown features to a minimum pixel width and height (using
  display parameters \code{min.width} and \code{min.height}) and
  collapsing overlapping annotation items (as defined by the parameter
  \code{min.distance} into one single item to
  prevent overplotting.

  For objects of class \code{\linkS4class{DataTrack}}, the data values
  underlying collapsed regions will be summarized based on the
  \code{summary} display parameter. See the class' documentation for
  more details. 

}


\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{DataTrack}}

  \code{\link{settings}}

}






