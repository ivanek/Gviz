\name{ImageMap-class}
\Rdversion{1.1}
\docType{class}
\alias{ImageMap-class}
\alias{coords,ImageMap-method}
\alias{tags,ImageMap-method}

\title{ImageMap class and methods}

\description{
  
  HTML image map information for annotation tracks.
  
}

\section{Objects from the Class}{

  Objects of the class are usually not created by the user, hence the
  constructor function \code{ImageMap} is not exported in the name
  space.
  
}
\section{Slots}{
  \describe{
    \item{\code{coords}:}{Object of class \code{"matrix"}, the image map
    coordinats. }
    \item{\code{tags}:}{Object of class \code{"list"}, the individual
      HTML tags for the image map. }
  }
}

\section{Extends}{
Class \code{"ImageMapOrNULL"}, directly.
}

\section{Methods}{
  
  \describe{
    
    \item{coords}{\code{signature(ImageMap="ImageMap")}: return the
      coordinates from the image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="ImageMap")}: return the tags from the
       image map.
      
      \emph{Usage:}
      
      \code{tags(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{tags(obj)}}

      }
    }
  }
}

\author{
Florian Hahne
}

\keyword{classes}





