\name{StackedTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{StackedTrack-class}
\alias{StackedTrack}
\alias{drawGD,StackedTrack-method}
\alias{initialize,StackedTrack-method}
\alias{setStacks,StackedTrack-method}
\alias{stacking,StackedTrack-method}
\alias{stacking<-,StackedTrack,character-method}
\alias{stacks,StackedTrack-method}
\alias{stacking}
\alias{stacking<-}
\alias{stacks}
\alias{subset,StackedTrack-method}
\alias{[,StackedTrack,ANY,ANY-method}

\title{StackedTrack class and methods}

\description{
  
  The virtual parent class for all track types in the Gviz
  package which contain potentially overlapping annotation items that
  have to be stacked when plotted.

}

\section{Objects from the Class}{

  A virtual Class: No objects may be created from it.
  
}

\section{Slots}{
  
  \describe{
  
    \item{\code{stacking}:}{Object of class \code{"character"}, the
      stacking type of overlapping items on the final plot.  One in
      \code{c(hide, dense, squish, pack,full)}. Currently, only
      \code{hide} (do not show the track items at all), \code{squish}
      (make best use of the available space) and \code{dense} (no
      stacking at all) are implemented. }
  
    \item{\code{stacks}:}{Object of class \code{"numeric"}, holding the
      stack indices for each track item. This slot is usually populated
      by calling the \code{setStacks} method upon plotting, since the
      correct stacking is a function of the available plotting space.}

    \item{\code{range}:}{Object of class
      \code{\linkS4class{GRanges}}, inherited from class
      \code{\linkS4class{RangeTrack}} }

    \item{\code{chromosome}:}{Object of class \code{"character"},
      inherited from class \code{\linkS4class{RangeTrack}} }

    \item{\code{genome}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{RangeTrack}} }

    \item{\code{dp}:}{Object of class
      \code{\linkS4class{DisplayPars}}, inherited from class
      \code{\linkS4class{GdObject}} }

    \item{\code{name}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{GdObject}}}

    \item{\code{imageMap}:}{Object of class \code{\linkS4class{ImageMap}},
      inherited from class \code{\linkS4class{GdObject}}}
  
  }
}

\section{Extends}{
  
  Class \code{"\linkS4class{RangeTrack}"}, directly.
  
  Class \code{"\linkS4class{GdObject}"}, by class "RangeTrack", distance 2.
  
}

\section{Methods}{

  In the following code chunks, \code{obj} is considered to be an object
  of class \code{StackedTrack}.
  
  \bold{\emph{Exported in the name space:}}

  \describe{
	  
    \item{stacking}{\code{signature(GdObject="StackedTrack")}: return
      the current stacking type.
    	
      \emph{Usage:}
      
      \code{stacking(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{stacking(obj)}}

      }
    }
    	
    \item{stacking<-}{\code{signature(GdObject="StackedTrack",
      value="character")}: set the object's stacking type to one in
      \code{c(hide, dense, squish, pack,full)}.
    	
      \emph{Usage:}
      
      \code{stacking<-(GdObject, value)}

      \emph{Additional Arguments:}
      
	\describe{
	
	  \item{}{\code{value}: replacement value.}
	  
	}

      \emph{Examples:}
      
      \describe{

	\item{}{\code{stacking(obj) <- "squish" }}

      }
    }
  }
    
  \bold{\emph{Internal methods:}}

  \describe{

    \item{drawGD}{\code{signature(GdObject="StackedTrack")}: plot the
      object to a graphics device. The return value of this method is
      the input object, potentially updated during the plotting
      operation. Internally, there are two modes in which the method can
      be called. Either in 'prepare' mode, in which case no plotting is
      done but the stacking information is updated based on the
      available space, or in 'plotting' mode, in which case the actual
      graphical output is created. Note that the method for this
      particular subclass is usually called through inheritance and not
      particularly useful on its own.
  
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
	
	\item{}{\code{Gviz:::drawGD(obj, prepare=FALSE)}}
	  
      }
    }

    \item{setStacks}{\code{signature(GdObject="StackedTrack")}:
      recompute the stacks based on the available space and on the
      object's track items and stacking settings.
      
      \emph{Usage:}
      
      \code{setStacks(GdObject, from, to)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: compute stacking within a
	  certain coordinates range. This needs to be supplied for the
	  plotting function to know the current genomic coordinates.}
	
      }

      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::setStacks(obj)}}

      }
    }
    
    \item{stacks}{\code{signature(GdObject="StackedTrack")}: return
      the stack indices for each track item.

      \emph{Usage:}
      
      \code{stacks(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::stacks(obj)}}

      }
    }
    
    \item{initialize}{\code{signature(.Object="StackedTrack")}:
      initialize the object. }
  
  }

  \bold{\emph{Inherited methods:}}
  
  \describe{
    
    \item{[}{\code{signature(x="StackedTrack", i="ANY", j="ANY",
	drop="ANY")}: subset the items in the \code{StackedTrack}
      object. This is essentially similar to subsetting of the
      \code{\linkS4class{GRanges}} object in the \code{range}
      slot. For most applications, the \code{subset} method may be
      more appropriate.
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{i}, \code{j}: subsetting indices, \code{j} is
	  ignored.}
	
	\item{}{\code{drop}: argument is ignored.}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{obj[1:5]}}
	
      }
    }
      
      
    \item{chromosome}{\code{signature(GdObject="StackedTrack")}:
      return the chromosome for which the track is defined.
      
      \emph{Usage:}
      
      \code{chromosome(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj)}}
	
      }
    }
    
    \item{chromosome<-}{\code{signature(GdObject="StackedTrack")}:
      replace the value of the track's chromosome. This has to be a
      valid UCSC chromosome identifier or an integer or character
      scalar that can be reasonably coerced into one.
      
      \emph{Usage:}
      
      \code{chromosome<-(GdObject, value)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{value}: replacement value.}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj) <- "chr12"}}
	
      }
    }
    
    \item{start, end, width}{\code{signature(x="StackedTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="StackedTrack")}:
      replace the start or end coordinates of the track items, or their
      width.
      
      \emph{Usage:}
      
      \code{start<-(x, value)}

      \code{end<-(x, value)}

      \code{width<-(x, value)}

      \emph{Additional Arguments:}
      
	\describe{
	
	  \item{}{\code{value}: replacement value.}
	  
	}

      \emph{Examples:}
      
      \describe{

	\item{}{\code{start(obj) <- 1:10}}
      
	\item{}{\code{end(obj) <- 20:30}}
      
	\item{}{\code{width(obj) <- 1}}

      }
    }
    
    \item{position}{\code{signature(GdObject="StackedTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{feature}{\code{signature(GdObject="StackedTrack")}: return the
      grouping information for track items. For certain sub-classes,
      groups may be indicated by different color schemes when
      plotting. See \code{\link{grouping}} or
      \code{\linkS4class{AnnotationTrack}} and
      \code{\linkS4class{GeneRegionTrack}} for details.
      
      \emph{Usage:}
      
      \code{feature(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{feature(obj)}}

      }
    }
    
    \item{feature<-}{\code{signature(gdObject="StackedTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{StackedTrack}. See \code{\link{grouping}} or
      \code{\linkS4class{AnnotationTrack}} and
      \code{\linkS4class{GeneRegionTrack}} for details.
      
      \emph{Usage:}
      
      \code{feature<-(GdObject, value)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{value}: replacement value.}
	  
      }
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{feature(obj) <- c("a", "a", "b", "c", "a")}}

      }
    }
    
    \item{genome}{\code{signature(x="StackedTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="StackedTrack")}: set the track's
      genome. Usually this has to be a valid UCSC identifier, however
      this is not formally enforced here.
	    
      \emph{Usage:}
      
      \code{genome<-(x, value)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{value}: replacement value.}
	  
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{genome(obj) <- "mm9"}}

      }
    }

    \item{length}{\code{signature(x="StackedTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="StackedTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="StackedTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{split}{\code{signature(x="StackedTrack")}: split a
      \code{StackedTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of objects of the same class as the input
      object, all inheriting from class \code{StackedTrack}.
      
      \emph{Usage:}
      
      \code{split(x, f, ...)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{f}: the splitting factor.}
	
	\item{}{\code{\dots}: all further arguments are ignored.}
	
      }

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{split(obj, c("a", "a", "b", "c", "a"))}}

      }
    }
    
    \item{strand}{\code{signature(x="StackedTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="StackedTrack")}: replace the
      strand information for the track items. The replacement value
      needs to be an appropriate scalar or vector of strand values.
      
      \emph{Usage:}
      
      \code{strand<-(x, value)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{value}: replacement value.}
	  
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{strand(obj) <- "+"}}

      }
    }
    
    \item{values}{\code{signature(x="StackedTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="StackedTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{subset}{\code{signature(x="StackedTrack")}: subset a
      \code{StackedTrack} by coordinates and sort if necessary.
      
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
    
    \item{displayPars}{\code{signature(x="StackedTrack",
      name="character")}: list the value of the display parameter
      \code{name}. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Usage:}
      
      \code{displayPars(x, name)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj, "col")}}
	
      }
    }
      
    
    \item{displayPars}{\code{signature(x="StackedTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="StackedTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="StackedTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="StackedTrack",
      value="list")}: set display parameters using the values of the
      named list in \code{value}. See \code{\link{settings}} for details
      on display parameters and customization.
      
      \emph{Usage:}
      
      \code{displayPars<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj) <- list(col="red", lwd=2)}}

      }
    }
    
    \item{setPar}{\code{signature(x="StackedTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{StackedTrack} class are
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
    
    \item{setPar}{\code{signature(x="StackedTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{StackedTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="StackedTrack")}: return
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
    
    \item{names}{\code{signature(x="StackedTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="StackedTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="StackedTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="StackedTrack")}: return the tags from the
      internal image map.
      
      \emph{Usage:}
      
      \code{tags(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{tags(obj)}}

      }
    }
  }
}


\section{Display Parameters}{

  The following display parameters are set for objects of class
  \code{AnnotationTrack} upon instantiation, unless one or more of them
  have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters
  for tracks.

  \describe{ 

    \item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse the
      y-ordering of stacked items. I.e., features that are plotted on
      the bottom-most stacks will be moved to the top-most stack and
      vice versa.}

  }

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{StackedTrack} objects. 

    \describe{ 

    \item{}{\code{\linkS4class{GdObject}}: 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panels.} 

        \item{}{\code{cex=1}: Numeric scalar. The overall font expansion
          factor for all text.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col="#0080FF"}: Integer or character scalar.
          Default line color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

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

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text.} 

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

        \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
          for all plotting elements, unless there is a more specific
          control defined elsewhere.} 

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

        \item{}{\code{showTitle=TRUE}: Boolean controlling whether to
          plot a title panel. Although this can be set individually
          for each track, in multi-track plots as created by
          \code{\link{plotTracks}} there will still be an empty
          placeholder in case any of the other tracks include a title.
          The same holds true for axes. Note that the the title panel
          background color could be set to transparent in order to
          completely hide it.} 

        \item{}{\code{size=1}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function.} 

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

  \code{\linkS4class{GeneRegionTrack}}

  \code{\linkS4class{GRanges}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{IRanges}}

  \code{\linkS4class{RangeTrack}}

  \code{\link{collapsing}}

  \code{\link{DataTrack}}

  \code{\link{grouping}}

  \code{\link{panel.grid}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}


\keyword{classes}






























