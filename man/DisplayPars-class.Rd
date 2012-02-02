\name{DisplayPars-class}
\Rdversion{1.1}
\docType{class}
\alias{DisplayPars-class}
\alias{DisplayPars}
\alias{displayPars}
\alias{displayPars<-}
\alias{displayPars,DisplayPars,character-method}
\alias{displayPars,DisplayPars,missing-method}
\alias{displayPars<-,DisplayPars,list-method}
\alias{getPar}
\alias{getPar,DisplayPars,character-method}
\alias{getPar,DisplayPars,missing-method}
\alias{initialize,DisplayPars-method}
\alias{setPar}
\alias{setPar,DisplayPars,character-method}
\alias{setPar,DisplayPars,list-method}
\alias{availableDisplayPars}
\alias{as.list,InferredDisplayPars-method}
\alias{as.list,DisplayPars-method}
\alias{show,DisplayPars-method}
\alias{coerce,DisplayPars,list-method}
\alias{coerce,InferredDisplayPars,list-method}
\alias{head,InferredDisplayPars-method}
\alias{tail,InferredDisplayPars-method}

\title{DisplayPars class and method}

\description{

  All tracks within this package are highly customizable. The
  \code{DisplayPars} class facilitates this and provides a unified API
  to the customization parameters.
  
}

\section{Objects from the Class}{

   Objects can be created using the constructor function
  \code{DisplayPars}.
  
}



\usage{
DisplayPars(...)

availableDisplayPars(class)
}

\arguments{

  \item{\dots}{All named arguments are stored in the object's
    environment as individual parameters, regardless of their type.}

  \item{class}{A valid track object class name, or the object itself, in
    which case the class is derived directly from it.}
  
}

\details{

  The individual parameters in a \code{DisplayParameters} class are
  stored as pointers in an environment. This has the upshot of not
  having to copy the whole track object when changing parameters, and
  parameters can be updated without the need to explicietly reassign the
  track to a symbol (i.e., updating of parameters happens in place). The
  downside is that upon copying of track objects, the parameter
  emvironment needs to be reinstantiated.

  The default display parameters for a track object class can be queried
  using the \code{availableDisplayPars} function.
  
}


\value{
  
  The return value of the constructor function is a new object of class
  \code{DisplayPars}.
  
  \code{availableDisplayPars} returns a list of the default display
  parameters.
  
}

\section{Slots}{
  \describe{
    \item{\code{pars}:}{Object of class \code{"environment"}, the
      container for all customization parameters. }
  }
}

\section{Methods}{
  
  In the following code chunks, \code{obj} is considered to be an object
  of class \code{DisplayPars}.
  
  \bold{\emph{Exported in the name space:}}
  
  \describe{
    
    \item{displayPars}{\code{signature(x="DisplayPars",
	name="character")}: return the value of a subset of display
      parameters, as identified by \code{name}.
      
      \emph{Usage:}
      
      \code{displayPars(x, name)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj, c("foo", "bar"))}}
      
	\item{}{\code{displayPars(obj, "foobar")}}
    
      }
    }
    
    \item{displayPars}{\code{signature(x="DisplayPars",
	name="missing")}: return all available display parameters.
  
      \emph{Usage:}
      
      \code{displayPars(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}
      
      }
    }

    \item{getPar}{\code{signature(x="DisplayPars", name="character")}:
      alias for the \code{displayPars} method.

      \emph{Usage:}
      
      \code{getPar(x, name)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj, "col")}}

      }
    }
    
    \item{getPar}{\code{signature(x="DisplayPars", name="missing")}: alias
      for the \code{displayPars} method.

      \emph{Usage:}
      
      \code{getPar(x)}
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="DisplayPars",
	value="list")}: replace or add display parameters as provided by
      the named list items.
    
      \emph{Usage:}
      
      \code{displayPars<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj) <- list(foo="a", bar=2)}}
      
      }
    }

    \item{setPar}{\code{signature(x="DisplayPars", value="character")}: set
      the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{DisplayPars} class are
      pass-by-reference, so no re-assignmnet to the symbol \code{obj} is
      necessary.
      
      \emph{Usage:}
      
      \code{setPar(x, name, value)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{name}: the name of the display parameter to set.}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, "col", "red")}}

      }
    }
    
    \item{setPar}{\code{signature(x="DisplayPars", value="list")}: set
      display parameters by the values of the named list in
      \code{value}. Note that display parameters in the \code{DisplayPars}
      class are pass-by-reference, so no re-assignmnet to the symbol
      \code{obj} is necessary.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
  }

  \emph{Internal methods:}
  
  \describe{
    
    \item{initialize}{\code{signature(.Object = "DisplayPars")}:
      initialize the object. }
    
    \item{show}{\code{signature(object = "DisplayPars")}:  show a
      human-readable summary of the object. }
    
  }
}


\author{
Florian Hahne
}



\examples{

## Construct object
dp <- DisplayPars(col="red", lwd=2, transformation=log2)
dp

## Query parameters
displayPars(dp)
displayPars(dp, "col")
getPar(dp, c("col", "transformation"))

## Modify parameters
displayPars(dp) <- list(lty=1, fontsize=3)
setPar(dp, "pch", 20)
dp

## Default parameters
availableDisplayPars("GenomeAxisTrack")

}

\keyword{classes}





