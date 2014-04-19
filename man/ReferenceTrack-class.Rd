\name{ReferenceTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{ReferenceTrack-class}
\alias{ReferenceTrack}
\alias{availableDefaultMapping}
\alias{initialize,ReferenceTrack-method}

\title{ReferenceTrack class and methods}

\description{
  
  A class allow for on-demand streaming of data off the file system.

}

\usage{

availableDefaultMapping(file, trackType)

}

\arguments{
  
  \item{file}{A character scalar with a file name or just a file
    extension.}

  \item{trackType}{A character scalar with one of the available track
  types in the package.}


}

\section{Objects from the class}{

  A virtual class: No objects may be created from it.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{stream}:}{Object of class \code{function}. The import
      function to stream data of the file system. Needs to be able to
      handle the two mandatory arguments \code{file} (a \code{character}
      containing a valid file path) and \code{selection} (a
      \code{GRanges} object with the genomic region to plot).}
   
    \item{\code{reference}:}{Object of class \code{"character"}, the
      path to the file containing the data.}
    
    \item{\code{mapping}:}{Object of class \code{"list"}, a default
      mapping between the metadata columns of the returned
      \code{GRanges} object from the import function and the
      \code{elemenMetadata} columns that make up the final track object.}
    
    \item{\code{args}:}{Object of class \code{"list"}, the passed in
      constructor arguments during object instantiation. Those will be
      needed when fetching the data in order to fill all necessary
      slots.}

    \item{\code{defaults}:}{Object of class \code{"list"}, the relevant
      default values to be used when neither \code{mapping} nor
      \code{args} provides the necessary information.}
    
  }
  
}


\section{Methods}{
  
  \bold{\emph{Internal methods:}}	

  \describe{
    
    \item{initialize}{\code{signature(.Object="ReferenceTrack")}:
      initialize the object.}
    
  }
}

\details{

  The \code{availableDefaultMappings} function can be used to find out
  whether the package defines a mapping scheme between one of the many
  supported input file types and the metadata columns of
  the tracks's \code{GRanges} objects.
  
}


\author{Florian Hahne}



\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\linkS4class{GeneRegionTrack}}

  \code{\linkS4class{GRanges}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{IRanges}}

  \code{\linkS4class{RangeTrack}}

  \code{\link{DataTrack}}

}

\keyword{classes}






























