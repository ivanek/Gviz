\name{GdObject-class}
\Rdversion{1.1}
\docType{class}
\alias{GdObject-class}
\alias{GdObject}
\alias{displayPars,GdObject,character-method}
\alias{displayPars,GdObject,missing-method}
\alias{displayPars<-,GdObject,list-method}
\alias{drawAxis,GdObject-method}
\alias{drawGrid,GdObject-method}
\alias{subset,GdObject-method}
\alias{getPar,GdObject,character-method}
\alias{getPar,GdObject,missing-method}
\alias{group,GdObject-method}
\alias{initialize,GdObject-method}
\alias{names,GdObject-method}
\alias{drawGD}
\alias{drawGD,DetailsAnnotationTrack-method}
\alias{drawGD,GeneRegionTrack-method}
\alias{drawGD,SequenceTrack-method}
\alias{names<-,GdObject,character-method}
\alias{genome<-,GdObject-method}
\alias{chromosome<-,GdObject-method}
\alias{isActiveSeq<-,GdObject-method}
\alias{setPar,GdObject,character-method}
\alias{setPar,GdObject,list-method}
\alias{coords,GdObject-method}
\alias{coords,NULL-method}
\alias{tags,GdObject-method}
\alias{chromosome,GdObject-method}
\alias{tags,NULL-method}
\alias{coords}
\alias{tags}
\alias{imageMap,GdObject-method}
\alias{imageMap}

\title{GdObject class and methods}

\description{
  
  The virtual parent class for all track items in the Gviz
  package.  This class definition contains all the common entities that
  are needed for a track to be plotted. During object instantiation for
  any of the sub-classes inheriting from \code{GdObject}, this class'
  global ininitializer has to be called in order to assure that all
  necessary settings are present.
  
}

\section{Objects from the class}{

  A virtual class: No objects may be created from it.
  
}
	
\section{Slots}{
  
  \describe{
  
    \item{\code{dp}:}{Object of class
      \code{\linkS4class{DisplayPars}}, the display settings
      controlling the look and feel of a track. See
      \code{\link{settings}} for details on setting graphical parameters
      for tracks.}
    
    \item{\code{name}:}{Object of class \code{"character"}, a
      human-readable name for the track that will be used in the track's
      annotation panel if necessary. }
    
    \item{\code{imageMap}:}{Object of class
      \code{\linkS4class{ImageMap}}, containing optional information
      for an HTML image map. This will be created by the \code{drawGD}
      methods when the track is plotted to a device and is usually not set
      by the user.}
    
  }
}


\section{Methods}{
  
  In the following code chunks, \code{obj} is considered to be an object
  of class \code{GdObject}.

  \bold{\emph{Exported in the name space:}}
    
  \describe{
    
    \item{displayPars}{\code{signature(x="GdObject", name="character")}:
      list the value of the display parameter \code{name}. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
      \emph{Usage:}
      
      \code{displayPars(x, name)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj, "col")}}

      }
    }
      
    
    \item{displayPars}{\code{signature(x="GdObject", name="missing")}:
      list the value of all available display parameters. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="GdObject", name="character")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Usage:}
      
      \code{getPar(x, name)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj, "col")}}

      }
    }
    
    \item{getPar}{\code{signature(x="GdObject", name="missing")}: alias
      for the \code{displayPars} method. See \code{\link{settings}} for
      details on display parameters and customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="GdObject", value="list")}:
      set display parameters using the values of the named list in
      \code{value}. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Usage:}
      
      \code{displayPars<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj) <- list(col="red", lwd=2)}}

      }
    }
    
    \item{setPar}{\code{signature(x="GdObject", value="character")}: set
      the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{GdObject} class are
      pass-by-reference, so no re-assignmnet to the symbol \code{obj} is
      necessary. See \code{\link{settings}} for details on display
      parameters and customization.
      
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
    
    \item{setPar}{\code{signature(x="GdObject", value="list")}: set
      display parameters by the values of the named list in
      \code{value}. Note that display parameters in the \code{GdObject}
      class are pass-by-reference, so no re-assignmnet to the symbol
      \code{obj} is necessary. See \code{\link{settings}} for details on
      display parameters and customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="GdObject")}: return
      grouping information for the individual items in the track. Unless
      overwritten in one of the sub-classes, this usualy returns
      \code{NULL}.
      
      \emph{Usage:}
      
      \code{group(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{group(obj)}}

      }
    }
    
    \item{names}{\code{signature(x="GdObject")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="GdObject", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="GdObject")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="GdObject")}: return the tags from the
      internal image map.
      
      \emph{Usage:}
      
      \code{tags(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{tags(obj)}}

      }
    }
    
    \item{subset}{\code{signature(x="GdObject")}: subset a
      \code{GdObject} by coordinates. Most of the respective sub-classes
      inheriting from \code{GdObject} overwrite this method, the default
      is to return the unaltered input object.
      
      \emph{Usage:}
      
      \code{subset(x, ...)}

      
      \emph{Additional Arguments:}
      
      \describe{

	\item{}{\code{\dots}: all further arguments are ignored.}

      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{subset(obj)}}

      }
    }
  }

  \bold{\emph{Internal methods:}}
    
  \describe{
    
    \item{drawAxis}{\code{signature(GdObject="GdObject")}: add a
      y-axis to the title panel of a track if necessary. Unless
      overwritten in one of the sub-classes this usualy does not plot
      anything and returns \code{NULL}.
      
      \emph{Usage:}
      
      \code{drawAxis(x, ...)}

      \emph{Additional Arguments:}
      
      \describe{

	\item{}{\code{\dots}: all further arguments are ignored.}

      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{Gviz:::drawAxis(obj)}}

      }
    }
    
    \item{drawGrid}{\code{signature(GdObject="GdObject")}: superpose a
      grid on top of a track if necessary. Unless overwritten in one of
      the sub-classes this usualy does not plot anything and returns
      \code{NULL}.
      
      \emph{Usage:}
      
      \code{drawGrid(GdObject, ...)}

      
      \emph{Additional Arguments:}
      
      \describe{

	\item{}{\code{\dots}: additional arguments are ignored.}

      }

      \emph{Examples:}
      
      \describe{

	\item{}{\code{Gviz:::drawGrid(obj)}}

      }
    }
    
    \item{initialize}{\code{signature(.Object="GdObject")}: initialize
      the object. This involves setting up a new environment for the
      display parameters and filling it up with the current
      settings. All arguments that have not been clobbered up by one of
      the sub-class initializers are considered to be additional display
      parameters and are also added to the environment.  See
      \code{\link{settings}} for details on setting graphical parameters
      for tracks.}
    
  }
}


\section{Display Parameters}{ 

  The following display parameters are set for objects of class
  \code{GdObject} upon instantiation, unless one or more of them have
  already been set by one of the optional sub-class initializers, which
  always get precedence over these global defaults. See \code{\link{settings}}
  for details on setting graphical parameters for tracks. \describe{ 

    \item{}{\code{alpha=1}: Numeric scalar. The transparency for all
      track items.} 

    \item{}{\code{background.panel="transparent"}: Integer or character
      scalar. The background color of the content panel.} 

    \item{}{\code{background.title="lightgray"}: Integer or character
      scalar. The background color for the title panels.} 

    \item{}{\code{cex=1}: Numeric scalar. The overall font expansion
      factor for all text.} 

    \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion factor
      for the axis annotation. Defaults to \code{NULL}, in which case
      it is computed based on the available space.} 

    \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion factor
      for the title panel. This effects the fontsize of both the title
      and the axis, if any. Defaults to \code{NULL}, which means that
      the text size is automatically adjusted to the available space.} 

    \item{}{\code{col="#0080FF"}: Integer or character scalar. Default
      line color setting for all plotting elements, unless there is a
      more specific control defined elsewhere.} 

    \item{}{\code{col.axis="white"}: Integer or character scalar. The
      font and line color for the y axis, if any.} 

    \item{}{\code{col.frame="lightgray"}: Integer or character scalar.
      The line color used for the panel frame, if \code{frame==TRUE}} 

    \item{}{\code{col.grid="#808080"}: Integer or character scalar.
      Default line color for grid lines, both when \code{type=="g"} in
      \code{\link{DataTrack}}s and when display parameter \code{grid==TRUE}.} 

    \item{}{\code{col.line=NULL}: Integer or character scalar. Default
      colors for plot lines. Usually the same as the global \code{col}
      parameter.} 

    \item{}{\code{col.symbol=NULL}: Integer or character scalar. Default
      colors for plot symbols. Usually the same as the global \code{col}
      parameter.} 

    \item{}{\code{col.title="white"}: Integer or character scalar. The
      font color for the title panels.} 

    \item{}{\code{collapse=TRUE}: Boolean controlling wether to collapse
      the content of the track to accomodate the minimum current device
      resolution. See \code{\link{collapsing}} for details.} 

    \item{}{\code{fill="lightgray"}: Integer or character scalar.
      Default fill color setting for all plotting elements, unless
      there is a more specific control defined elsewhere.} 

    \item{}{\code{fontcolor="black"}: Integer or character scalar. The
      font color for all text.} 

    \item{}{\code{fontface=1}: Integer or character scalar. The font
      face for all text.} 

    \item{}{\code{fontface.title=2}: Integer or character scalar. The
      font face for the title panels.} 

    \item{}{\code{fontfamily="sans"}: Integer or character scalar. The
      font family for all text.} 

    \item{}{\code{fontfamily.title="sans"}: Integer or character scalar.
      The font family for the title panels.} 

    \item{}{\code{fontsize=12}: Numeric scalar. The font size for all
      text.} 

    \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the track
      when plotting.} 

    \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
      of a grid.} 

    \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
      number of horizontal grid lines, see \code{\link{panel.grid}}
      for details.} 

    \item{}{\code{lineheight=1}: Numeric scalar. The font line height
      for all text.} 

    \item{}{\code{lty="solid"}: Numeric scalar. Default line type
      setting for all plotting elements, unless there is a more specific
      control defined elsewhere.} 

    \item{}{\code{lty.grid="solid"}: Integer or character scalar.
      Default line type for grid lines, both when \code{type=="g"} in
      \code{\link{DataTrack}}s and when display parameter \code{grid==TRUE}.} 

    \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
      for all plotting elements, unless there is a more specific control
      defined elsewhere.} 

    \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width for
      grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
      and when display parameter \code{grid==TRUE}.} 

    \item{}{\code{min.distance=1}: Numeric scalar. The minimum pixel
      distance before collapsing range items, only if \code{collapse==TRUE}.
      See \code{\link{collapsing}} for details.} 

    \item{}{\code{min.height=3}: Numeric scalar. The minimum range
      height in pixels to display. All ranges are expanded to this size
      in order to avoid rendering issues. See \code{\link{collapsing}}
      for details.} 

    \item{}{\code{min.width=1}: Numeric scalar. The minimum range width
      in pixels to display. All ranges are expanded to this size in
      order to avoid rendering issues. See \code{\link{collapsing}}
      for details.} 

    \item{}{\code{showAxis=TRUE}: Boolean controlling whether to plot
      a y axis (only applies to track types where axes are implemented).} 

    \item{}{\code{showTitle=TRUE}: Boolean controlling whether to plot
      a title panel. Although this can be set individually for each
      track, in multi-track plots as created by \code{\link{plotTracks}}
      there will still be an empty placeholder in case any of the other
      tracks include a title. The same holds true for axes. Note that
      the the title panel background color could be set to transparent
      in order to completely hide it.} 

    \item{}{\code{size=1}: Numeric scalar. The relative size of the
      track. Can be overridden in the \code{\link{plotTracks}} function.} 

    \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
      number of vertical grid lines, see \code{\link{panel.grid}} for
      details.} 

  } 

}


\author{Florian Hahne}


\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GeneRegionTrack}}

  \code{\linkS4class{ImageMap}}

  \code{\link{collapsing}}

  \code{\link{DataTrack}}

  \code{\link{panel.grid}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}

\keyword{classes}


































