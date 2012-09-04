\name{IdeogramTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{IdeogramTrack-class}
\alias{IdeogramTrack}
\alias{drawGD,IdeogramTrack-method}
\alias{end,IdeogramTrack-method}
\alias{end<-,IdeogramTrack-method}
\alias{initialize,IdeogramTrack-method}
\alias{show,IdeogramTrack-method}
\alias{start,IdeogramTrack-method}
\alias{start<-,IdeogramTrack-method}
\alias{subset,IdeogramTrack-method}
\alias{width,IdeogramTrack-method}
\alias{width<-,IdeogramTrack-method}
\alias{length,IdeogramTrack-method}
\alias{[,IdeogramTrack-method}
\alias{[,IdeogramTrack,ANY,ANY-method}
\alias{chromosome<-,IdeogramTrack-method}
\alias{genome<-,IdeogramTrack-method}
\alias{position,IdeogramTrack-method}

\title{IdeogramTrack class and methods}

\description{
  
  A class to represent the schematic display of a chromosome, also known
  as an ideogram. The respective information is typically directly
  fetched from UCSC.

}

\section{Objects from the Class}{

  Objects can be created using the constructor function \code{IdeogramTrack}.

}

\usage{

IdeogramTrack(chromosome=NULL, genome, name=NULL, bands=NULL, ...)

}

\arguments{
  
  \item{chromosome}{The chromosome for which to create the ideogram. Has
    to be a valid UCSC chromosome identifier of the form \code{chrx}, or
    a single integer or numeric character unless
    \code{option(ucscChromosomeNames=FALSE)}. The user has to make sure
    that the respective chromosome is indeed defined for the the track's
    genome.}

  \item{genome}{The genome on which to create the ideogram. This has to
    be a valid UCSC genome identifier if the ideogram data is to be
    fetched from the UCSC repository.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting. Defaults to the selected chromosome.}

  \item{bands}{A \code{data.frame} with the cytoband information for all
    available chromosomes on the genome similar to the data that would
    be fetched from UCSC. The table needs to contain the mandatory
    columns \code{chrom}, \code{chromStart}, \code{chromEnd},
    \code{name} and \code{gieStain} with the chromosome name, cytoband
    start and end coordinates, cytoband name and coloring information,
    respectively. This can be used when no connection to the internet is
    available or when the cytoband information has been cached locally
    to avoid the somewhat slow connection to UCSC.}
    
  \item{\dots}{Additional items which will all be interpreted as further
    display parameters.}
  
}

\value{

  The return value of the constructor function is a new object of class
  \code{IdeogramTrack}.
  
}


\section{Slots}{
  
  \describe{
    
    \item{\code{range}:}{Object of class \code{\linkS4class{GRanges}},
      inherited from class \code{\linkS4class{StackedTrack}} containing
      the chromosome band information. This slot is filled automatically
      by the initializer method.}
    
    \item{\code{bandTable}:}{Object of
    class\code{\linkS4class{data.frame}} containing the chromosome band
    information in the format of UCSC. This slot is filled automatically
    by the initializer.}

    \item{\code{bandTable}:}{Object of class
      \code{\linkS4class{data.frame}} containing the chromosome band
      information for all chromosomes. This slot is filled automatically
      by the initializer method and only exists to prevent having to
      redo the rtracklayer query ervery time the chromosome is changed.}

    \item{\code{chromosome}:}{Object of class \code{"character"},
      inherited from class \code{\linkS4class{StackedTrack}} defining
      the ideogram's chromosome.}
    
    \item{\code{genome}:}{Object of class \code{"character"} inherited
      from class \code{\linkS4class{StackedTrack}} defining the
      ideogram's genome.}
    
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
  
  Class \code{"\linkS4class{RangeTrack}"}, directly.
  
  Class \code{"\linkS4class{GdObject}"}, by class "RangeTrack", distance
  2.

}

\section{Methods}{
  
   In the following code chunks, \code{obj} is considered to be an object of class \code{IdeogramTrack}.
  
  \bold{\emph{Exported in the name space:}}
  
   
  \describe{

    \item{start, end, width, position}{\code{signature(x/GdObject="IdeogramTrack")}:
      although \code{IdeogramTracks} inherit from \code{RangeTrack},
      the notion of coordinates is not particularly useful. Hence the
      coordinate methods all return \code{NULL}.
      
      \emph{Usage:}

      \code{start(x)}

      \code{end(x)}
      
      \code{width(x)}

      \code{position(GdObject)}

      \emph{Examples:}

      \describe{
       
	\item{}{\code{start(obj)}}
	
	\item{}{\code{end(obj)}}
      
	\item{}{\code{width(obj)}}

      }
    }

    \item{start<-, end<-, width<-}{\code{signature(x="RangeTrack")}:
      although \code{IdeogramTracks} inherit from \code{RangeTrack}, the
      notion of coordinates is not particularly useful. Hence the
      coordinate replacement methods all return the unaltered input
      object.
      
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

    \item{chromosome<-}{\code{signature(GdObject="IdeogramTrack")}:
      replace the value of the track's chromosome. This has to be a
      valid UCSC chromosome identifier or an integer or character
      scalar that can be reasonably coerced into one. The chromosome
      band information is updated automatically.
      
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

    \item{genome<-}{\code{signature(x="IdeogramTrack")}: set the track's
      genome. This has to be a valid UCSC identifier. The chromosome
      band information is updated automatically.
	    
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

    
    \item{subset}{\code{signature(x="IdeogramTrack")}: subsetting does
      not make much sense for these object, hence the unalered object is
      returned.
      
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

    \item{[}{\code{signature(x="IdeogramTrack", i="ANY", j="ANY",
	drop="ANY")}: subsetting of \code{IdeogramTrack} objects does
      not make much sense, hence the unaltered input argument is
      returned.
      
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
  }
  
  \bold{\emph{Internal methods:}}
  
  \describe{

    \item{drawGD}{\code{signature(gdObject="IdeogramTrack")}: plot the
      object to a graphics device. The return value of this method is
      the input object, potentially updated during the plotting
      operation. Internally, there are two modes in which the method can
      be called. Either in 'prepare' mode, in which case no plotting is
      done but the object is preprocessed based on the
      available space, or in 'plotting' mode, in which case the actual
      graphical output is created. Since subsetting of the object can be
      potentially costly, this can be switched off in case subsetting
      has already been performed before or is not necessary.
      
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
    
    \item{initialize}{\code{signature(.Object="IdeogramTrack")}:
      initialize the object. }
      
    \item{show}{\code{signature(object="IdeogramTrack")}: show a
	human-readable summary of the object. }
    }
    
    \bold{\emph{Inherited methods:}}
    
    \describe{
            
      \item{chromosome}{\code{signature(GdObject="IdeogramTrack")}:
	return the chromosome for which the track is defined.
      
	\emph{Usage:}
      
	\code{chromosome(GdObject)}

        \emph{Examples:}
      
	\describe{

	  \item{}{\code{chromosome(obj)}}

	}
      }
    
     
   
    \item{feature}{\code{signature(GdObject="IdeogramTrack")}: return the
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
    
    \item{feature<-}{\code{signature(gdObject="IdeogramTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{IdeogramTrack}. See \code{\link{grouping}} or
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
    
    \item{genome}{\code{signature(x="IdeogramTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{length}{\code{signature(x="IdeogramTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="IdeogramTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="IdeogramTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{split}{\code{signature(x="IdeogramTrack")}: splitting is not a
      useful operation for \code{IdeogramTrack} objects.  \emph{Usage:}
      
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
    
    \item{strand}{\code{signature(x="IdeogramTrack")}: strand
    information is not relevant for \code{IdeogramTrack} objects.
      
      \emph{Usage:}
      
      \code{strand(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{strand(obj)}}

      }
    }
    
    \item{strand<-}{\code{signature(x="IdeogramTrack")}: strand
      information is not relevant for \code{IdeogramTrack} objects.

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
    
    \item{values}{\code{signature(x="IdeogramTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="IdeogramTrack", to="data.frame")}:
	coerce the \code{\linkS4class{GRanges}} object in the
	\code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }
 
    \item{displayPars}{\code{signature(x="IdeogramTrack",
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
      
    
    \item{displayPars}{\code{signature(x="IdeogramTrack", name="missing")}:
      list the value of all available display parameters. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="IdeogramTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="IdeogramTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="IdeogramTrack", value="list")}:
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
    
    \item{setPar}{\code{signature(x="IdeogramTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{IdeogramTrack} class are
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
    
    \item{setPar}{\code{signature(x="IdeogramTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{IdeogramTrack} class are pass-by-reference, so no re-assignmnet
      to the symbol \code{obj} is necessary. See \code{\link{settings}}
      for details on display parameters and customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="IdeogramTrack")}: return
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
    
    \item{names}{\code{signature(x="IdeogramTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="IdeogramTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="IdeogramTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="IdeogramTrack")}: return the tags from the
      internal image map.
      
      \emph{Usage:}
      
      \code{tags(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{tags(obj)}}

      }
    }
    
    \item{drawAxis}{\code{signature(GdObject="IdeogramTrack")}: add a
      y-axis to the title panel of a track if necessary. For
      \code{IdeogramTrack} objects this does not plot anything and
      returns \code{NULL}.
      
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
    
    \item{drawGrid}{\code{signature(GdObject="IdeogramTrack")}:
      superpose a grid on top of a track if necessary. For
      \code{IdeogramTrack} objects this does not plot anything and
      returns \code{NULL}.
      
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
  
\details{

  Ideograms are schematic depictions of chromosomes, including
  chromosome band information and centromer location. The relevant data
  for various species is stored in the UCSC data base. The initializer
  method of the class will automatically fetch the respective data for a
  given genome and chromosome from UCSC and fill the appropriate object
  slots. When plotting \code{IdeogramTracks}, the current genomic
  location is indicated on the chromosome by a colored box.

}


\section{Display Parameters}{ 

  The following display parameters are set for objects of class
  \code{IdeogramTrack} upon instantiation, unless one or more of them
  have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters
  for tracks. \describe{ 

    \item{}{\code{background.title="transparent"}: Character scalar.
      The background color for the title panel. Defaults to omit the
      background.} 

    \item{}{\code{bevel=0.45}: Numeric scalar, between 0 and 1. The
      level of smoothness for the two ends of the ideogram.} 

    \item{}{\code{cex=0.8}: Numeric scalar. The overall font expansion
      factor for the chromosome name text.} 

    \item{}{\code{col="red"}: Character scalar. The border color used
      for the highlighting of the currently displayed genomic region.} 

    \item{}{\code{fill="#FFE3E6"}: Character scalar. The fill color
      used for the highlighting of the currently displayed genomic
      region.} 

    \item{}{\code{fontcolor="#808080"}: Character scalar. The font
      color for the chromosome name text.} 

    \item{}{\code{fontsize=10}: Numeric scalar. The font size for the
      chromosome name text.} 

    \item{}{\code{showId=TRUE}: Logical scalar. Indicate the chromosome
      name next to the ideogram.} 

    \item{}{\code{showTitle=FALSE}: Logical scalar. Plot a title panel.
      Defaults to omit the title panel.} 

    \item{}{\code{size=NULL}: Numeric scalar. The relative size of the
      track. Defaults to automatic size setting. Can be overridden in
      the \code{\link{plotTracks}} function.} 

  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{IdeogramTrack} objects. 

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

  \code{\linkS4class{StackedTrack}}

  \code{\link{clearSessionCache}}

  \code{\link{collapsing}}

  \code{\link{DataTrack}}

  \code{\link{grouping}}

  \code{\link{panel.grid}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}

\note{

  When fetching ideogram data from UCSC the results are cached for
  faster acces. See \code{\link{clearSessionCache}} on details to delete
  these cached items.
  
}

\examples{

\dontshow{
## Load some sample data
data(idTrack)
}

## Construct the object
\dontrun{
idTrack <- IdeogramTrack(chromosome=7, genome="mm9")
}

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(idTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
plotTracks(idTrack, from=5000000, to=9000000)

## Track names
names(idTrack)
names(idTrack) <- "foo"
plotTracks(idTrack, from=5000000, to=9000000)


## Accessors
chromosome(idTrack)
\dontrun{
chromosome(idTrack) <- "chrX"
}

genome(idTrack)
\dontrun{
genome(id) <- "hg19"
}

range(idTrack)
ranges(idTrack)

## Annotation
values(idTrack)

## coercion
as(idTrack, "data.frame")

}

\keyword{classes}






























