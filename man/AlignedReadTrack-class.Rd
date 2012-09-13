\name{AlignedReadTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{AlignedReadTrack-class}
\alias{AlignedReadTrack}
\alias{[,AlignedReadTrack-method}
\alias{[,AlignedReadTrack,ANY,ANY-method}
\alias{coverage,AlignedReadTrack-method}
\alias{drawAxis,AlignedReadTrack-method}
\alias{drawGD,AlignedReadTrack-method}
\alias{drawGrid,AlignedReadTrack-method}
\alias{initialize,AlignedReadTrack-method}
\alias{setCoverage,AlignedReadTrack-method}
\alias{show,AlignedReadTrack-method}
\alias{split,AlignedReadTrack-method}
\alias{subset,AlignedReadTrack-method}

\title{AlignedReadTrack class and methods (NOTE: THIS IS STILL IN
  DEVELOPMENT AND SUBJECT TO CHANGE)}

\description{
  
  A class to represent short sequences that have been aligned to a
  reference genome as they are typically generated in a next generation
  sequencing experiment.
  
}

\section{Objects from the Class}{

  Objects can be created using the constructor function \code{AlignedReadTrack}.

}


\usage{

AlignedReadTrack(range=NULL, start=NULL, end=NULL, width=NULL, chromosome, strand="+", genome,
stacking="squish", name="AlignedReadTrack", coverageOnly=FALSE, ...)

}

\arguments{

  \item{range}{An object of class \code{\linkS4class{GRanges}}, or a
    \code{data.frame} which will be coerced into one in which case it
    needs to contain at least the three columns:
		
    \describe{

      \item{}{\code{start}, \code{end}: the start and end coordinates
	for the track items.}
      
      \item{}{\code{strand}: the strand information for the track
	items. It may be provided in the form \code{+} for the Watson
	strand, \code{-} for the Crick strand or \code{*} for either one
	of the two.}
      
    }
    
    Alternatively, the \code{range} argument may be missing, in which
    case the relevant information has to be provided as individual
    function arguments (see below).
    
  }
	
  \item{start, end, width}{Integer vectors, giving the start and the end
    end coordinates for the individual track items, or their width. Two
    of the three need to be specified, and have to be of equal length or
    of length one, in which case this value will be recycled. Otherwise,
    the usual R recycling rules for vectors do not apply.}

  \item{strand}{Character vector, the strand information for the
    individual track items.  Needs to be of equal length as the
    \code{start, end} or \code{width} vectors, or of length 1. Please
    note that grouped items need to be on the same strand, and erroneous
    entries will result in casting of an error.}

  \item{chromosome}{The chromosome on which the track's genomic ranges
    are defined. A valid UCSC chromosome identifier. Please note that at
    this stage only syntactic checking takes place, i.e., the argument
    value needs to be a single integer, numeric character or a character
    of the form \code{chrx}, where \emph{x} may be any
    possible string. The user has to make sure that the respective
    chromosome is indeed defined for the the track's genome.}

  \item{genome}{The genome on which the track's ranges are
    defined. Usually this is a valid UCSC genome identifier, however
    this is not being formally checked at this point.}

  \item{stacking}{The stacking type for overlapping items of the
    track. One in \code{c(hide, dense, squish, pack,full)}. Currently,
    only hide (don't show the track items, squish (make best use of the
    available space) and dense (no stacking at all) are implemented.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

  \item{coverageOnly}{Instead of storing individual reads, just compute
    the coverage and store the resulting coverage vector.}

  \item{\dots}{Additional items which will all be interpreted as further
    display parameters.}
  
}
 
\value{

  The return value of the constructor function is a new object of class
  \code{AlignedReadTrack}.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{coverage}:}{Object of class \code{"list"}, a list of
      coverage vectors for the plus strand, the minus strand and for
      both strands combined.}

    \item{\code{coverageOnly}:}{Object of class \code{"logical"}, flag
      to determine whether the object stores read locations or the coverage
      vectors only.}

     \item{\code{stacking}:}{Object of class \code{"character"},
      inherited from class \code{\linkS4class{StackedTrack}}}

    \item{\code{stacks}:}{Object of class \code{"environment"},
      inherited from class \code{\linkS4class{StackedTrack}}}

    \item{\code{range}:}{Object of class \code{\linkS4class{GRanges}},
      inherited from class \code{\linkS4class{RangeTrack}}}

    \item{\code{chromosome}:}{Object of class \code{"character"},
      inherited from class \code{\linkS4class{RangeTrack}} }

    \item{\code{genome}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{RangeTrack}}}

    \item{\code{dp}:}{Object of class
      \code{\linkS4class{DisplayPars}}, inherited from class
      \code{\linkS4class{GdObject}} }

    \item{\code{name}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{GdObject}}}

    \item{\code{imageMap}:}{Object of class
      \code{\linkS4class{ImageMap}}, inherited from class
      \code{\linkS4class{GdObject}}}
    
  }  
}

\section{Extends}{
  
  Class \code{"\linkS4class{StackedTrack}"}, directly.
  
  Class \code{"\linkS4class{RangeTrack}"}, by class "StackedTrack",
  distance2.
  
  Class \code{"\linkS4class{GdObject}"}, by class "StackedTrack",
  distance3.
  
}

\section{Methods}{
  
  In the following code chunks, \code{obj} is considered to be an object
  of class \code{AlignedReadTrack}.
  
  \bold{\emph{Exported in the name space:}}
  
  \describe{
    
    \item{[}{\code{signature(x="AlignedReadTrack")}: subset the items in
      the \code{AlignedReadTrack}. This is essentially similar to
      subsetting of the \code{\linkS4class{GRanges}} object in the
      \code{range} slot. For most applications, the \code{subset} method
      may be more appropriate. The operation is only supported for
      objects that still contain all the read locations, i.e.,
      \code{coverageOnly=FALSE}.

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{i}: subsetting indices}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{obj[1:5]}}
	
      }
    }

    \item{subset}{\code{signature(x="AlignedReadTrack")}: subset a
      \code{AlignedReadTrack} by coordinates and sort if necessary.
      
      \emph{Usage:}
        
      \code{subset(x, from, to, sort=FALSE, stacks=FALSE)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: the coordinates range to subset
	  to.}
	
	\item{}{\code{sort}: sort the object after subsetting. Usually
	  not necessary.}
	
	\item{}{\code{stacks}: recompute the stacking after subsetting
	  which can be expensive and is not always necessary.}
	
      }
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{subset(obj, from=10, to=20)}}
	
	\item{}{\code{subset(obj, from=10, to=20, sort=TRUE, stacks=FALSE)}}

      }
    }
      
    \item{split}{\code{signature(x="AlignedReadTrack")}: split an
      \code{AlignedReadTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of \code{AlignedReadTrack} objects.
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{f}: the splitting factor.}
	
	\item{}{\code{\dots}: all further arguments are ignored.}
	
      }
      
      \emph{Usage:}
      
      \code{split(x, f, ...)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{split(obj, c("a", "a", "b", "c", "a"))}}
	
      }
    }
      
    \item{coverage}{\code{signature(x="AlignedReadTrack")}: return the
      coverage vector for one of the strands, or the combined vector.

      \emph{Usage:}
        
      \code{coverage(x, strand="*")}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{strand}: the selector for the strand, \code{+} for
	  the Watson strand, \code{-} for the Crick strand or \code{*}
	  for both strands.}
	
      }
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{coveraget(obj)}}
	    
	\item{}{\code{coverage(obj, strand="-")}}
	
      }
    }
  }

  \bold{\emph{Internal methods:}}
  
  \describe{
    
    
    \item{setCoverage}{\code{signature(GdObject="AlignedReadTrack")}:
      recompute the coverage on the plus and minus strand as well as for
      the combined strands and update the respective slot.
      
      \emph{Usage:}
      
      \code{setCoverage(GdObject)}

       \emph{Examples:}
      
      \describe{

	\item{}{\code{setCoverage(obj)}}

      }
    }

    \item{drawAxis}{\code{signature(GdObject="AlignedReadTrack")}: add a y-axis
      to the title panel of a track.

      \emph{Usage:}
      
      \code{drawAxis(GdObject, from, to, subset=FALSE, ...)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: compute axis range from the data
	  within a certain coordinates range only.}

	\item{}{\code{subset}: subset the object prior to calculating
	  the axis ranges. Can be expensive and is not always needed.}
	
	\item{}{\code{\dots}: all further arguments are ignored.}
	
      }

      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::drawAxis(obj)}}
	
      }
    }
   
    \item{drawGD}{\code{signature(gdObject="AlignedReadTrack")}: plot the
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

    \item{drawGrid}{\code{signature(GdObject="AlignedReadTrack")}: superpose a grid on top of a track.

      \emph{Usage:}
      
      \code{drawGrid(GdObject, from, to)}
       
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: draw grid within a certain
	  coordinates range. This needs to be supplied for the plotting
	  function to know the current genomic coordinates.}
	
      }

      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::drawGrid(obj, from=10, to=100)}}
	
      }
    }

     
    \item{initialize}{\code{signature(.Object="AligendReadTrack")}:
      initialize the object.}

    \item{show}{\code{signature(object="AlignedReadTrack")}: show a
      human-readable summary of the object.}
        
  }

  \bold{\emph{Inherited methods:}}

  \describe{
	  
    \item{stacking}{\code{signature(GdObject="AlignedReadTrack")}: return
      the current stacking type.
    	
      \emph{Usage:}
      
      \code{stacking(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{stacking(obj)}}

      }
    }
    	
    \item{stacking<-}{\code{signature(GdObject="AlignedReadTrack",
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
    
    \item{setStacks}{\code{signature(GdObject="AlignedReadTrack")}:
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
    
    \item{stacks}{\code{signature(GdObject="AlignedReadTrack")}: return
      the stack indices for each track item.

      \emph{Usage:}
      
      \code{stacks(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::stacks(obj)}}

      }
    }
    
    \item{chromosome}{\code{signature(GdObject="AlignedReadTrack")}:
      return the chromosome for which the track is defined.
      
      \emph{Usage:}
      
      \code{chromosome(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj)}}
	
      }
    }
    
    \item{chromosome<-}{\code{signature(GdObject="AlignedReadTrack")}:
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
    
    \item{start, end, width}{\code{signature(x="AlignedReadTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="AlignedReadTrack")}:
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
    
    \item{position}{\code{signature(GdObject="AlignedReadTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{feature}{\code{signature(GdObject="AlignedReadTrack")}: return the
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
    
    \item{feature<-}{\code{signature(gdObject="AlignedReadTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{AlignedReadTrack}. See \code{\link{grouping}} or
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
    
    \item{genome}{\code{signature(x="AlignedReadTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="AlignedReadTrack")}: set the track's
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

    \item{length}{\code{signature(x="AlignedReadTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="AlignedReadTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="AlignedReadTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{strand}{\code{signature(x="AlignedReadTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="AlignedReadTrack")}: replace the
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
    
    \item{values}{\code{signature(x="AlignedReadTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="AlignedReadTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{displayPars}{\code{signature(x="AlignedReadTrack",
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
      
    
    \item{displayPars}{\code{signature(x="AlignedReadTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="AlignedReadTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="AlignedReadTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="AlignedReadTrack",
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
    
    \item{setPar}{\code{signature(x="AlignedReadTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{AlignedReadTrack} class are
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
    
    \item{setPar}{\code{signature(x="AlignedReadTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{AlignedReadTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="AlignedReadTrack")}: return
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
    
    \item{names}{\code{signature(x="AlignedReadTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="AlignedReadTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="AlignedReadTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="AlignedReadTrack")}: return the tags from the
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
  \code{AlignedReadTrack} upon instantiation, unless one or more of
  them have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters
  for tracks. \describe{ 

    \item{}{\code{collapse=FALSE}: collapse overlapping ranges and
      aggregate the underlying data.} 

    \item{}{\code{detail="coverage"}: the amount of detail to plot the
      data. Either \code{coverage} to show the coverage only, or
      \code{reads} to show individual reads. For large data sets the
      latter can be very inefficient. Please note that \code{reads} is
      only available when the object has been created with option
      \code{coverageOnly=FALSE}.} 

    \item{}{\code{fill="#0080ff"}: the fill color for the coverage
      indicator.} 

    \item{}{\code{size=NULL}: the relative size of the track. Defaults
      to size selection based on the underlying data. Can be overridden
      in the \code{\link{plotTracks}} function.} 

    \item{}{\code{type="histogram"}: the plot type, one or several in
      \code{c("p","l", "b", "a", "s", "g", "r", "S", "smooth", "histogram",
      "mountain", "h", "boxplot", "gradient", "heatmap")}. See	the
      'Details' section in \code{\linkS4class{DataTrack}} for more
      information on the individual plotting types.} 

  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{AlignedReadTrack} objects. 

  \describe{

    \item{}{\code{\linkS4class{StackedTrack}}:
      
      \describe{ 

	\item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse the
	  y-ordering of stacked items. I.e., features that are plotted on
	  the bottom-most stacks will be moved to the top-most stack and
	  vice versa.}

      }
      
    }

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

  \code{\linkS4class{DataTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\linkS4class{GeneRegionTrack}}

  \code{\linkS4class{GRanges}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{IRanges}}

  \code{\linkS4class{RangeTrack}}

  \code{\linkS4class{StackedTrack}}

  \code{\link{collapsing}}

  \code{\link{grouping}}

  \code{\link{panel.grid}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}

\examples{

## Construct from individual arguments
arTrack <- AlignedReadTrack(start=runif(1000, 100, 200), width=24,
genome="mm9", chromosome=7, strand=sample(c("+", "-"), 1000, TRUE))

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(arTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
plotTracks(arTrack)

## Track names
names(arTrack)
names(arTrack) <- "foo"
plotTracks(arTrack)

## Subsetting and splitting
subTrack <- subset(arTrack, from=110, to=130)
length(subTrack)
subTrack[1:2]
split(arTrack, strand(arTrack))

## Accessors
start(arTrack)
end(arTrack)
width(arTrack)
position(arTrack)
width(subTrack) <- 30

strand(arTrack)
strand(subTrack) <- "-"

chromosome(arTrack)
chromosome(subTrack) <- "chrX"

genome(arTrack)
genome(subTrack) <- "mm9"

range(arTrack)
ranges(arTrack)

coverage(arTrack)

## Annotation
values(arTrack)

## Stacking
stacking(arTrack)
stacking(arTrack) <- "dense"

## coercion
as(arTrack, "data.frame")

}
\keyword{classes}































