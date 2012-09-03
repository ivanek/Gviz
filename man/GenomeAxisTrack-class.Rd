\name{GenomeAxisTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomeAxisTrack-class}
\alias{GenomeAxisTrack}
\alias{[,GenomeAxisTrack-method}
\alias{[,GenomeAxisTrack,ANY,ANY-method}
\alias{collapseTrack,GenomeAxisTrack-method}
\alias{drawGD,GenomeAxisTrack-method}
\alias{end,GenomeAxisTrack-method}
\alias{end<-,GenomeAxisTrack-method}
\alias{initialize,GenomeAxisTrack-method}
\alias{length,GenomeAxisTrack-method}
\alias{range,GenomeAxisTrack-method}
\alias{ranges,GenomeAxisTrack-method}
\alias{ranges<-,GenomeAxisTrack-method}
\alias{show,GenomeAxisTrack-method}
\alias{start,GenomeAxisTrack-method}
\alias{start<-,GenomeAxisTrack-method}
\alias{strand,GenomeAxisTrack-method}
\alias{subset,GenomeAxisTrack-method}
\alias{values,GenomeAxisTrack-method}
\alias{width,GenomeAxisTrack-method}

\title{GenomeAxisTrack class and methods}

\description{A class representing a customizable genomic axis.}

\section{Objects from the class}{

  Objects can be created using the constructor function
  \code{GenomeAxisTrack}.

}

\usage{
GenomeAxisTrack(range=NULL, name="Axis", id, ...)
}

\arguments{
  
  \item{range}{Optional \code{\linkS4class{GRanges}} or
    \code{\linkS4class{IRanges}} object to highlight certain regions on
    the axis.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

   \item{id}{A character vector of the same length as \code{range}
   containing identifiers for the ranges. If missing, the constructor
   will try to extract the ids from \code{names(range)}.}

  \item{\dots}{Additional items which will all be interpreted as further
    display parameters. See \code{\link{settings}} and the "Display
    Parameters" section below for details.}
  
}

\details{
  
  A \code{GenomeAxisTrack} can be customized using the familiar display
  parameters. By providing a \code{GRanges} or \code{IRanges} object to
  the constructor, ranges on the axis can be further highlighted.\\
  
  With the \code{scale} display parameter, a small scale indicator can
  be shown instead of the entire genomic axis. The scale can either be
  provided as a fraction of the plotting region (it will be rounded to
  the nearest human readable absolute value) or as an absolute value and
  is always displayed in bp, kb, mb or gb units. Note that most display
  parameters for the \code{GenomeAxisTrack} are ignored when a scale is
  used insterad of the full axis. In particular, only the parameters
  \code{exponent}, \code{alpha}, \code{lwd}, \code{col}, \code{cex},
  \code{distFromAxis} and \code{labelPos} are used.
  
}


\value{
  
  The return value of the constructor function is a new object of class
  \code{GenomeAxisTrack}.
  
}

\section{Slots}{
  \describe{
    
    \item{\code{range}:}{Object of class \code{\linkS4class{GRanges}},
      highlighted on the axis. }

    \item{\code{dp}:}{Object of class
      \code{\linkS4class{DisplayPars}}, inherited from class
      \code{\linkS4class{GdObject}}}
    
    \item{\code{name}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{GdObject}}}
    
     \item{\code{imageMap}:}{Object of class
      \code{\linkS4class{ImageMap}}, inherited from class
      \code{\linkS4class{GdObject}}}

  }
}

\section{Extends}{
Class \code{"\linkS4class{GdObject}"}, directly.
}

\section{Methods}{
  
  In the following code chunks, \code{obj} is considered to be an object of class \code{GenomeAxisTrack}.
	
  \bold{\emph{Exported in the name space:}}
  
  \describe{
    
    \item{[}{\code{signature(x="GenomeAxisTrack")}: subset the
      \code{GRanges} object in the \code{range} slot. For most
      applications, the \code{subset} method may be more appropriate.

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{i}: subsetting incides.}

      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{obj[1:5]}}

      }
    }

    \item{start, end, width}{\code{signature(x="GenomeAxisTrack")}: the
      start or end coordinates of the track items, or their width in
      genomic coordinates.
      
      \emph{Usage:}
      
      \code{start(x)}

      \code{end(x)}

      \code{width(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{start(obj)}}
      
	\item{}{\code{end(obj)}}
      
	\item{}{\code{width(obj)}}
	
      }
    }

    \item{range}{\code{signature(x="GenomeAxisTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="GenomeAxisTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }

    \item{strand}{\code{signature(x="GenomeAxisTrack")}: return a vector of
      strand specifiers for all track items, in the form '+' for the
      Watson strand, '-' for the Crick strand or '*' for either of the
      two.

      \emph{Usage:}
      
      \code{strand(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{strand(obj)}}

      }
    }

    \item{values}{\code{signature(x="GenomeAxisTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items.
      
      \emph{Usage:}
      
      \code{values(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }

    \item{subset}{\code{signature(x="GenomeAxisTrack")}: subset a
      \code{GenomeAxisTrack} by coordinates and sort if necessary.
      
      
      \emph{Usage:}
        
      \code{subset(x, from, to, sort=FALSE, ...)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: the coordinates range to subset
	  to.}
	
	\item{}{\code{sort}: sort the object after subsetting. Usually
	  not necessary.}
	
	\item{}{\code{\dots}: additional arguments are ignored.}
	
      }
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{subset(obj, from=10, to=20, sort=TRUE)}}

      }
    }
        
    \item{length}{\code{signature(x="GenomeAxisTrack")}: return the
      number of items stored in the \code{ranges} slot.

      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}
	
      }
    }
  }
  
  \emph{Internal methods:}
  
  \describe{

    \item{drawGD}{\code{signature(GdObject="GenomeAxisTrack")}: the
      workhorse function to plot the object.

      \emph{Usage:}

      \code{drawGD(GdObject, minBase, maxBase, prepare=FALSE,
	subset=TRUE, ...)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{minBase}, \code{maxBase}: the coordinate range to
	  plot.}
	
	\item{}{\code{prepare}: run method in preparation or in
	  production mode.}
	
	\item{}{\code{subset}: subset the object to the visible region
	  or skip the potentially expensive subsetting operation.}
	
	\item{}{\code{\dots}: all further arguments are ignored.}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{Gviz:::drawGD(obj)}}
	
	\item{}{\code{Gviz:::drawGD(obj, minBase=1, maxBase=100)}}
	
	\item{}{\code{Gviz:::drawGD(obj, prepare=TRUE,
	    subset=FALSE)}}
	
      }
    }
      
    \item{collapseTrack}{\code{signature(GdObject="GenomeAxisTrack")}:
      preprocess the track before plotting. This will collapse
      overlapping track items based on the available resolution and
      increase the width and height of all track objects to a minimum
      value to avoid rendering issues. See \code{\link{collapsing}} for
      details.
      
      \emph{Usage:}

      \code{collapseTrack(GdObject, diff=.pxResolution(coord="x"))}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{diff}: the minimum pixel width to display,
	  everything below that will be inflated to a width of
	  \code{diff}.}
	
      }

      \emph{Examples:}
      
      \describe{
  
	\item{}{\code{Gviz:::collapseTrack(obj)}}

      }
    }
    
    \item{initialize}{\code{signature(.Object="GenomeAxisTrack")}:
      initialize the object }
    
    \item{show}{\code{signature(object="GenomeAxisTrack")}: show a
      human-readable summary of the object }
    
  }

   \bold{\emph{Inherited:}}
    
  \describe{
    
    \item{displayPars}{\code{signature(x="GenomeAxisTrack", name="character")}:
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
      
    
    \item{displayPars}{\code{signature(x="GenomeAxisTrack", name="missing")}:
      list the value of all available display parameters. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="GenomeAxisTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="GenomeAxisTrack", name="missing")}: alias
      for the \code{displayPars} method. See \code{\link{settings}} for
      details on display parameters and customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="GenomeAxisTrack", value="list")}:
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
    
    \item{setPar}{\code{signature(x="GenomeAxisTrack", value="character")}: set
      the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{GenomeAxisTrack} class are
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
    
    \item{setPar}{\code{signature(x="GenomeAxisTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}. Note that display parameters in the \code{GenomeAxisTrack}
      class are pass-by-reference, so no re-assignmnet to the symbol
      \code{obj} is necessary. See \code{\link{settings}} for details on
      display parameters and customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="GenomeAxisTrack")}: return
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
    
    \item{names}{\code{signature(x="GenomeAxisTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="GenomeAxisTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="GenomeAxisTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="GenomeAxisTrack")}: return the tags from the
      internal image map.
      
      \emph{Usage:}
      
      \code{tags(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{tags(obj)}}

      }
    }
    
    \item{drawAxis}{\code{signature(GdObject="GenomeAxisTrack")}: add a
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
    
    \item{drawGrid}{\code{signature(GdObject="GenomeAxisTrack")}: superpose a
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
  }
}


\section{Display Parameters}{ 

  The following display parameters are set for objects of class
  \code{GenomeAxisTrack} upon instantiation, unless one or more of them
  have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters
  for tracks. \describe{ 

    \item{}{\code{add35=FALSE}: Logical scalar. Add 3' to 5' direction
      indicators.} 

    \item{}{\code{add53=FALSE}: Logical scalar. Add 5' to 3' direction
      indicators.} 

    \item{}{\code{background.title="transparent"}: Character scalar.
      The background color for the title panel. Defaults to omit the
      background.} 

    \item{}{\code{cex=0.8}: Numeric scalar. The overall font expansion
      factor for the axis annotation text.} 

    \item{}{\code{cex.id=0.7}: Numeric scalar. The text size for the
      optional range annotation.} 

    \item{}{\code{col="darkgray"}: Character scalar. The color for the
      axis lines and tickmarks.} 

    \item{}{\code{col.id="white"}: Character scalar. The text color
      for the optional range annotation.} 

    \item{}{\code{col.range="cornsilk4"}: Character scalar. The border
      color for highlighted regions on the axis.} 

    \item{}{\code{distFromAxis=1}: Numeric scalar. Control the distance
      of the axis annotation from the tick marks.} 

    \item{}{\code{exponent=NULL}: Numeric scalar. The exponent for the
      axis coordinates, e.g., 3 means mb, 6 means gb, etc. The default
      is to automatically determine the optimal exponent.} 

    \item{}{\code{fill.range="cornsilk3"}: Character scalar. The fill
      color for highlighted regions on the axis.} 

    \item{}{\code{fontcolor="#808080"}: Character scalar. The font
      color for the axis annotation text.} 

    \item{}{\code{fontsize=10}: Numeric scalar. Font size for the axis
      annotation text in points.} 

    \item{}{\code{labelPos="alternating"}: Character vector, one in
      "alternating", "revAlternating", "above" or "below". The vertical
      positioning of the axis labels. If \code{scale} is not \code{NULL},
      the possible values are "above", "below" and "beside". } 

    \item{}{\code{littleTicks=FALSE}: Logical scalar. Add more fine-grained
      tick marks.} 

    \item{}{\code{lwd=2}: Numeric scalar. The line width for the axis
      elementes.} 

    \item{}{\code{showId=FALSE}: Logical scalar. Show the optional
      range highlighting annotation.} 

    \item{}{\code{showTitle=FALSE}: Logical scalar. Plot a title panel.
      Defaults to omit the title panel.} 

    \item{}{\code{size=NULL}: Numeric scalar. The relative size of the
      track. Can be overridden in the \code{\link{plotTracks}} function.
      Defaults to the ideal size based on the other track settings.} 

    \item{}{\code{scale=NULL}: Numeric scalar. If not \code{NULL} a
     small scale is drawn instead of the full axis, if the value is
     between 0 and 1 it is interpreted as a fraction of the current
     plotting region, otherwise as an absolute length value in genomic
     coordinates.}

  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{GenomeAxisTrack} objects. 

    \describe{ 

    \item{}{\code{\linkS4class{GdObject}}: 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{fill="lightgray"}: Integer or character scalar.
          Default fill color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

      } 

    } 

  } 

}
  
\author{Florian Hahne}


\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\linkS4class{GRanges}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{IRanges}}

  \code{\linkS4class{RangeTrack}}

  \code{\linkS4class{StackedTrack}}

  \code{\link{collapsing}}

  \code{\link{DataTrack}}

  \code{\link{grouping}}

  \code{\link{panel.grid}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}
\examples{

## Construct object
axTrack <- GenomeAxisTrack(name="Axis",
range=IRanges::IRanges(start=c(100, 300, 800), end=c(150, 400, 1000)))

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(axTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}


## Plotting
plotTracks(axTrack, from=0, to=1100)

## Track names
names(axTrack)
names(axTrack) <- "foo"

## Subsetting and splitting
subTrack <- subset(axTrack, from=0, to=500)
length(subTrack)
subTrack[1]
split(axTrack, c(1,1,2))

## Accessors
start(axTrack)
end(axTrack)
width(axTrack)

strand(axTrack)

range(axTrack)
ranges(axTrack)

## Annotation
values(axTrack)

## Grouping
group(axTrack)

## HTML image map
coords(axTrack)
tags(axTrack)
axTrack <- plotTracks(axTrack)$foo
coords(axTrack)
tags(axTrack)

## adding an axis to another track
data(cyp2b10)
grTrack <- GeneRegionTrack(start=26682683, end=26711643,
		rstart=cyp2b10$start, rends=cyp2b10$end, chromosome=7, genome="mm9",
		transcript=cyp2b10$transcript, gene=cyp2b10$gene, symbol=cyp2b10$symbol,
		name="Cyp2b10", strand=cyp2b10$strand)
		
plotTracks(list(grTrack, GenomeAxisTrack()))
plotTracks(list(grTrack, GenomeAxisTrack(scale=0.1)))
plotTracks(list(grTrack, GenomeAxisTrack(scale=5000)))
plotTracks(list(grTrack, GenomeAxisTrack(scale=0.5, labelPos="below")))

}

\keyword{classes}































