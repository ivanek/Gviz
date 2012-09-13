\name{AnnotationTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{AnnotationTrack-class}
\alias{AnnotationTrack}
\alias{DetailsAnnotationTrack-class}
\alias{DetailsAnnotationTrack}
\alias{initialize,DetailsAnnotationTrack-method}
\alias{coerce,AnnotationTrack,UCSCData-method}
\alias{collapseTrack,AnnotationTrack-method}
\alias{drawGD,AnnotationTrack-method}
\alias{drawGrid,AnnotationTrack-method}
\alias{group}
\alias{group,AnnotationTrack-method}
\alias{group<-}
\alias{group<-,AnnotationTrack,character-method}
\alias{identifier,AnnotationTrack-method}
\alias{identifier<-,AnnotationTrack,character-method}
\alias{identifier}
\alias{identifier<-}
\alias{initialize,AnnotationTrack-method}
\alias{setStacks,AnnotationTrack-method}
\alias{show,AnnotationTrack-method}
\alias{subset,AnnotationTrack-method}
\alias{coerce,GRanges,AnnotationTrack-method}
\alias{coerce,GRangesList,AnnotationTrack-method}

\title{AnnotationTrack class and methods}

\description{
  
  A fairly generic track object for arbitrary genomic range annotations,
  with the option of grouped track items. The extended
  \code{DetailsAnnotationTrack} provides a flexible interface to add
  user-defined custom information for each range.

}

\section{Objects from the class}{
  
  Objects can be created using the constructor function
  \code{AnnotationTrack}.

}

\usage{

AnnotationTrack(range=NULL, start=NULL, end=NULL, width=NULL, feature, group,
 id, strand="*", chromosome, genome, stacking="squish", name="AnnotationTrack",
 fun, selectFun, ...)

}

\arguments{

  \code{AnnotationTrack} object can be created from a variety of
  different inputs in order to nicely embed the package into the
  existing Bioconductor landscape. Since the main components of this
  class are essentially genomic ranges, the obvious Bioconductor
  representation is most likely a \code{\linkS4class{GRanges}} object,
  or, for grouped elements, a \code{\linkS4class{GRangesList}}. However,
  in certain cases it may be desirable to build the object from
  individual function arguments.
  
  \item{range}{
    
    An optional meta argument to handle the different input types. If
    the \code{range} argument is missing, all the relevant information
    to create the object has to be provided as individual function
    arguments (see below).

    The different input options for \code{range} are:

    \describe{

      \item{}{A \code{GRanges} object: the genomic ranges for the
	\code{Annotation} track as well as the optional additional
	\code{elementMetadata} columns \code{feature}, \code{group} and
	\code{id} (see description of the individual function parameters
	below for details). Calling the constructor on a \code{GRanges}
	object without further arguments, e.g.
	\code{AnnotationTrack(range=obj)} is equivalent to calling the coerce
	method \code{as(obj, "AnnotationTrack")}.}

      \item{}{A \code{GRangesList} object: this is very similar to the
	previous case, except that the grouping information that is part
	of the list structure is preserved in the
	\code{AnnotationTrack}. I.e., all the elements within one list
	item receive the same group id. For consistancy, there is also a
	coercion method from \code{GRangesLists} \code{as(obj,
	"AnnotationTrack")}.}

      \item{}{An \code{\linkS4class{IRanges}} object: almost identical
	to the \code{GRanges} case, except that the chromosome and
	strand information as well as all additional data has to be
	provided in the separate \code{chromosome}, \code{strand},
	\code{feature}, \code{group} or \code{id} arguments, because it
	can not be directly encoded in an \code{IRange} object. Note
	that only the former two are mandatory (if not provided
	explicitely the more or less reasonable default values
	\code{chromosome=NA} and \code{strand=*} are used. }

      \item{}{A \code{data.frame} object: the \code{data.frame} needs to
	contain at least the two mandatory columns \code{start} and
	\code{end} with the range coordinates. It may also contain a
	\code{chromosome} and a \code{strand} column with the chromosome
	and strand information for each range. If missing it will be
	drawn from the separate \code{chromosome} or \code{strand}
	arguments. In addition, the \code{feature}, \code{group} and
	\code{id} data can be provided as additional columns. The above
	comments about potential default values also apply here.}

    }
  }
  
  \item{start, end, width}{Integer vectors, giving the start and the end
    end coordinates for the individual track items, or their width. Two
    of the three need to be specified, and have to be of equal length or
    of length one, in which case this single value will be
    recycled. Otherwise, the usual R recycling rules for vectors do not
    apply here.}

  \item{feature}{Factor (or other vector that can be coerced into one),
    giving the feature types for the individual track items. When
    plotting the track to the device, if a display parameter with the
    same name as the value of \code{feature} is set, this will be used
    as the track item's fill color. See \code{\link{grouping}} for
    details.  Needs to be of equal length as the provided genomic
    coordinates, or of length 1.}

  \item{group}{Factor (or other vector that can be coerced into one),
    giving the group memberships for the individual track items. When
    plotting to the device, all items in the same group will be
    connected.  See \code{\link{grouping}} for details. Needs to be of
    equal length as the provided genomic coordinates, or of length 1.}

  \item{id}{Character vector of track item identifiers. When plotting to
    the device, it's value will be used as the identifier tag if the
    display parameter \code{showFeatureId=TRUE}. Needs to be of equal
    length as the provided genomic ranges, or of length 1.}

  \item{strand}{Character vector, the strand information for the
    individual track items. It may be provided in the form \code{+} for
    the Watson strand, \code{-} for the Crick strand or \code{*} for
    either one of the two. Needs to be of equal length as the provided
    genomic coordinates, or of length 1. Please note that grouped items
    need to be on the same strand, and erroneous entries will result in
    casting of an error.}

  \item{chromosome}{The chromosome on which the track's genomic ranges
    are defined. A valid UCSC chromosome identifier if
    \code{options(ucscChromosomeNames=TRUE)}. Please note that in this
    case only syntactic checking takes place, i.e., the argument value
    needs to be an integer, numeric character or a character of the form
    \code{chrx}, where \code{x} may be any possible string. The user has
    to make sure that the respective chromosome is indeed defined for
    the the track's genome. If not provided here, the constructor will
    try to construct the chromosome information based on the available
    inputs, and as a last resort will fall back to the value
    \code{chrNA}. Please note that by definition all objects in the
    \code{Gviz} package can only have a single active chromosome at a
    time (although internally the information for more than one
    chromosome may be present), and the user has to call the
    \code{chromosome<-} replacement method in order to change to a
    different active chromosome.}

  \item{genome}{The genome on which the track's ranges are
    defined. Usually this is a valid UCSC genome identifier, however
    this is not being formally checked at this point. If not provided
    here the constructor will try to extract this information from the
    provided input, and eventually will fall back to the default value
    of \code{NA}.}

  \item{stacking}{The stacking type for overlapping items of the
    track. One in \code{c(hide, dense, squish, pack,full)}. Currently,
    only squish (make best use of the available space), dense (no
    stacking, collapse overlapping ranges), and hide (do not show any
    track items at all) are implemented.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

  \item{fun}{A function that is being called
    for each entry in the \code{AnnotationTrack} object. See section
    'Details' and 'Examples' for further information. When called
    internally by the plotting machinery, a number of arguments are
    automatically passed on to this function, and the user needs to make
    sure that they can all be digested (i.e., either have all of them as
    formal named function arguments, or gobble up everything that is not
    needed in \dots). These arguments are:
    
    \itemize{
      
      \item{\code{start}: the genomic start coordinate of the range
	item.}
      
      \item{\code{end}: the genomic end coordinates of the range
	item.}
      
      \item{\code{strand}: the strand information for the range item.}
      
      \item{\code{chromosome}: the chromosome of the range item.}
      
      \item{\code{identifier}: the identifier of the range item, i.e.,
	the result of calling \code{identifier(DetailsAnnotationTrack,
	  lowest=TRUE)}. Typically those identifiers are passed on to
	the object constructor during instantiation as the \code{id}
	argument.}

      \item{\code{index}: a counter enumerating the ranges. The
	\code{AnnotationTrack} object is sorted internally for
	visibility, and the \code{index} argument refers to the index
	of plotting.}

       \item{\code{GdObject}: a reference to the currently plotted
	 \code{DetailsAnnotationTrack} object.}

        \item{\code{GdObject.original}: a reference to the
          \code{DetailsAnnotationTrack} before any processing like item
          collapsing has taken place. Essentially, this is the track
          object as it exists in your working environment.}
      
    }
    
    Additional arguments can be passed to the plotting function by means
    of the \code{detailsFunArgs} argument (see below). Note that the
    plot must use grid graphics (e.g. function in the 'lattice' package
    or low-level grid functions). To access a data object such a matrix
    or data frame within the function you can either store it as a
    variable in the global environment or, to avoid name space
    conflicts, you can make it part of the function environment by means
    of a closure. Alternatively, you may want to explicitely stick it
    into an environment or pass it along in the \code{detailsFunArgs}
    list. To figure out in your custom plotting function which
    annotation element is currently being plotted you can either use the
    identifier which has to be unique for each range element, or you may
    want to use the genomic position (start/end/strand/chromosome)
    e.g. if the data is stored in a \code{GRanges} object.
  
  }

  \item{selectFun}{A function that is being called
    for each entry in the \code{AnnotationTrack} object with exactly
    the same arguments as in \code{fun}. The purpose of this function
    is to decide for each track element whether details should be
    drawn, and consequently it has to return a single logical
    scalar. If the return value is \code{TRUE}, details will be drawn
    for the item, if it is \code{FALSE}, the details strip for the
    item is omitted.
  
  }
  
  \item{\dots}{Additional items which will all be interpreted as further
    display parameters. See \code{\link{settings}} and the "Display
    Parameters" section below for details.}
}

\value{

  The return value of the constructor function is a new object of class
  \code{AnnotationTrack} or of class \code{DetailsAnnotationTrack},
  depending on the constructor arguments. Typically the user will not
  have to be troubled with this distinction and can rely on the
  constructor to make the right choice.
  
}


\section{details}{
	
  The \code{DetailsAnnotationTrack} class directly extends
  \code{AnnotationTrack}. The purpose of this track type is to add an
  arbitrarily detailed plot section (typically consisting of additional
  quantitative data) for each range element of an
  \code{AnnotationTrack}. This allows a locus wide view of annotation
  elements together with any kind of details per feature or element that
  may for instance provide insight on how some complex quantitative
  measurements change accoring to their position in a locus.  If the
  quantitative data is too complex for a \code{DataTrack} e.g. because
  it requires extra space or a trellis-like representation, a
  \code{DetailsAnnotationTrack} can be used instead. Example: An
  \code{AnnotationTrack} shows the positions of a number of probes from
  a microarray, and you want a histogram of the signal intensity
  distribution derived from all samples at each of these probe
  location. Another example usage would be to show for each element of
  an \code{AnnotationTrack} an xy-plot of the signal against some
  clinical measurement such as blood preassure. The limitation for
  applications of this type of track is basically only the available
  space of the device you are plotting to.

  This flexibility is possible by utilizing a simple function model to
  perform all the detailed plotting. The functionality of this plotting
  function \code{fun} is totally up to the user, and the function
  environment is prepared in a way that all necessary information about
  the plotted annotation feature is available. To restrict the details
  section to only selected number of annotation features one can supply
  another function \code{selectFun}, which decides for each feature
  separatly whether details are available or not. Finally, an arbitrary
  number of additional arguments can be passed on to these two function
  by means of the \code{detailsFunArgs} display parameter. This is
  expected to be a named list, and all list elements are passed along to
  the plotting function \code{fun} and to the selector function
  \code{selectFun} as additional named arguments. Please note that some
  argument names like \code{start}, \code{end} or \code{identifier} are
  reserved and can not be used in the \code{detailsFunArgs} list. For
  examples of plotting functions, see the 'Examples' section.
	
}


\section{Slots}{
  
  \describe{
  
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
  
    \item{\code{fun}:}{A function that is being called for each
      \code{AnnotationTrack} element to plot details.}

    \item{\code{selectFun}:}{A function that is being called for each
      \code{AnnotationTrack} element to decide whether details need to
      be plotted.}
  
  }
  
  Additional display parameters are being inherited from the
  \code{\linkS4class{StackedTrack}} parent class.
  
}

\section{Extends}{
  
  Class \code{"\linkS4class{StackedTrack}"}, directly.
  
  Class \code{"\linkS4class{RangeTrack}"}, by class "StackedTrack", distance 2.
  
  Class \code{"\linkS4class{GdObject}"}, by class "StackedTrack",
  distance3.

  \code{DetailsAnnotationTrack} directly extends \code{AnnotationTrack}.
  
}

\section{Methods}{

  In the following code chunks, \code{obj} is considered to be an object
  of class \code{AnnotationTrack} or \code{DetailsAnnotationTrack}.
  
  \bold{\emph{Exported in the name space:}}
  
  \describe{
    
    \item{group}{\code{signature(GdObject="AnnotationTrack")}: extract
      the group membership for all track items.
      
      \emph{Usage:}
        
      \code{group(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{group(obj)}}
	
      }
    }
		
    \item{group<-}{\code{signature(GdObject="AnnotationTrack",
      value="character")}: replace the grouping information for track
      items. The replacement value must be a factor of appropriate
      length or another vector that can be coerced into such.

      \emph{Usage:}
        
      \code{group<-(GdObject, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{group(obj) <- c("a", "a", "b", "c", "a")}}

      }
    }
	    
    \item{identifier}{\code{signature(GdObject="AnnotationTrack")}:
      return track item identifiers. Depending on the setting of the
      optional argument \code{lowest}, these are either the group
      identifiers or the individual item identifiers.

      \emph{Usage:}
        
      \code{identifier(GdObject, lowest=FALSE)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{lowest}: return the lowest-level identifier, i.e.,
	the item IDs, or the higher level group IDs which do not have to
	be unqiue.}
	
      }
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{identifier(obj)}}
	
	\item{}{\code{identifier(obj, lowest=TRUE)}}

      }
    }
	    
    \item{identifier<-}{\code{signature(GdObject="AnnotationTrack",
	value="character")}: Set the track item identifiers. The
      replacement value has to be a character vector of appropriate
      length. This always replaces the group-level identifiers, so
      essentially it is similar to \code{groups<-}.

      \emph{Usage:}
        
      \code{identifier<-(GdObject, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{identifier(obj) <- c("foo", "bar")}}

      }
    }
  }

  \bold{\emph{Internal methods:}}
  
  \describe{
    
    \item{coerce}{\code{signature(from="AnnotationTrack",
	to="UCSCData")}: coerce to a \code{UCSCData} object for export to
      the UCSC genome browser.

      \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "UCSCData")}}
	
      }
    }

    \item{collapseTrack}{\code{signature(GdObject="AnnotationTrack")}:
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

    \item{drawGD}{\code{signature(GdObject="AnnotationTrack")}: plot the
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

    \item{drawGrid}{\code{signature(GdObject="AnnotationTrack")}: superpose a grid on top of a track.

      \emph{Usage:}
      
      \code{drawGrid(GdObject, from, to)}
       
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: integer scalars, draw grid
	  within a certain coordinates range. This needs to be supplied
	  for the plotting function to know the current genomic
	  coordinates.}
	
      }

      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::drawGrid(obj, from=10, to=100)}}
	
      }
    }

    \item{setStacks}{\code{signature(GdObject="AnnotationTrack")}:
      recompute the stacks based on the available space and on the
      object's track items and stacking settings.

      \emph{Usage:}
      
      \code{setStacks(GdObject, from, to)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: integer scalars, compute
	  stacking within a certain coordinates range. This needs to be
	  supplied for the plotting function to know the current genomic
	  coordinates.}
	
      }

      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::setStacks(obj, from=1, to=100)}}

      }
    }
    
    \item{initialize}{\code{signature(.Object="AnnotationTrack")}:
      initialize the object }

    \item{show}{\code{signature(object="AnnotationTrack")}: show a
      human-readable summary of the object }
 
  }

  \bold{\emph{Inherited methods:}}

  \describe{
	  
    \item{stacking}{\code{signature(GdObject="AnnotationTrack")}: return
      the current stacking type.
    	
      \emph{Usage:}
      
      \code{stacking(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{stacking(obj)}}

      }
    }
    	
    \item{stacking<-}{\code{signature(GdObject="AnnotationTrack",
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

    \item{stacks}{\code{signature(GdObject="AnnotationTrack")}: return
      the stack indices for each track item.
      
      \emph{Usage:}
      
      \code{stacks(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::stacks(obj)}}
	
      }
    }
    
    \item{[}{\code{signature(x="AnnotationTrack", i="ANY", j="ANY",
	drop="ANY")}: subset the items in the \code{AnnotationTrack}
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
      
      
    \item{chromosome}{\code{signature(GdObject="AnnotationTrack")}:
       return the currently active chromosome for which the track is
       defined. For consistancy with other Bioconductor packages, the
       \code{isActiveSeq} alias is also provided.
      
      \emph{Usage:}
      
      \code{chromosome(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj)}}
	
      }
    }
    
    \item{chromosome<-}{\code{signature(GdObject="AnnotationTrack")}:
       replace the value of the track's active chromosome. This has to
       be a valid UCSC chromosome identifier or an integer or character
       scalar that can be reasonably coerced into one, unless
       \code{options(ucscChromosomeNames=FALSE)}. For consistancy with
       other Bioconductor packages, the \code{isActiveSeq<-} alias is
       also provided.
      
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
    
    \item{start, end, width}{\code{signature(x="AnnotationTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="AnnotationTrack")}:
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
    
    \item{position}{\code{signature(GdObject="AnnotationTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{feature}{\code{signature(GdObject="AnnotationTrack")}: return the
      grouping information for track items. For certain sub-classes,
      groups may be indicated by different color schemes when
      plotting. See \code{\link{grouping}} for details.
      
      \emph{Usage:}
      
      \code{feature(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{feature(obj)}}

      }
    }
    
    \item{feature<-}{\code{signature(gdObject="AnnotationTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{AnnotationTrack}. See \code{\link{grouping}}
      for details.
      
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
    
    \item{genome}{\code{signature(x="AnnotationTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="AnnotationTrack")}: set the track's
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

    \item{length}{\code{signature(x="AnnotationTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="AnnotationTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="AnnotationTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{split}{\code{signature(x="AnnotationTrack")}: split a
      \code{AnnotationTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of objects of the same class as the input
      object, all inheriting from class \code{AnnotationTrack}.
      
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
    
    \item{strand}{\code{signature(x="AnnotationTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="AnnotationTrack")}: replace the
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
    
    \item{values}{\code{signature(x="AnnotationTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="AnnotationTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{subset}{\code{signature(x="AnnotationTrack")}: subset a
      \code{AnnotationTrack} by coordinates and sort if necessary.
      
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
    
    \item{displayPars}{\code{signature(x="AnnotationTrack",
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
      
    
    \item{displayPars}{\code{signature(x="AnnotationTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="AnnotationTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="AnnotationTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="AnnotationTrack",
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
    
    \item{setPar}{\code{signature(x="AnnotationTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{AnnotationTrack} class are
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
    
    \item{setPar}{\code{signature(x="AnnotationTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{AnnotationTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{names}{\code{signature(x="AnnotationTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="AnnotationTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="AnnotationTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="AnnotationTrack")}: return the tags from the
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

    \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
      for item identifiers.} 

    \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
      factor for the group-level annotation.} 

    \item{}{\code{col="transparent"}: Character or integer scalar. The
      border color for all track items.}

    \item{}{\code{col.line="darkgray"}: Character scalar. The color used for
      connecting lines between grouped items. Defaults to a dark gray,
      but if set to \code{NULL} the same color as for the first item in
      the group is used.}

    \item{}{\code{fill="lightblue"}: Character or integer scalar. The
      fill color for untyped items. This is also used to connect grouped
      items. See \code{\link{grouping}} for details.} 

    \item{}{\code{fontcolor="white"}: Character or integer scalar. The
      font color for item identifiers.} 

    \item{}{\code{fontcolor.group="#808080"}: Character or integer
      scalar. The font color for the group-level annotation.} 

    \item{}{\code{fontface=1}: Integer scalar. The font face for item
      identifiers.} 

    \item{}{\code{fontface.group=2}: Numeric scalar. The font face for
      the group-level annotation.} 

    \item{}{\code{fontfamily="sans"}: Character scalar. The font family
      for item identifiers.} 

    \item{}{\code{fontsize=12}: Numeric scalar. The font size for item
      identifiers.} 

    \item{}{\code{lex=1}: Numeric scalar. The line expansion factor
      for all track items. This is also used to connect grouped items.
      See \code{\link{grouping}} for details.} 

    \item{}{\code{lineheight=1}: Numeric scalar. The font line height
      for item identifiers.} 

    \item{}{\code{lty="solid"}: Character or integer scalar. The line
      type for all track items. This is also used to connect grouped
      items. See \code{\link{grouping}} for details.} 

    \item{}{\code{lwd=1}: Integer scalar. The line width for all track
      items. This is also used to connect grouped items. See
      \code{\link{grouping}} for details.} 

    \item{}{\code{rotation=0}: Numeric scalar. The degree of text
      rotation for item identifiers.} 

    \item{}{\code{shape="arrow"}: Character scalar. The shape in which
      to display the track items. Currently only \code{box}, \code{arrow},
      \code{ellipse}, and \code{smallArrow} are implemented.} 

    \item{}{\code{showFeatureId=FALSE}: Logical scalar. Control whether
      to plot the individual track item identifiers.} 

    \item{}{\code{showId=FALSE}: Logical scalar. Control whether to
      annotate individual groups.} 

    \item{}{\code{showOverplotting=FALSE}: Logical scalar. Use a color
      gradient to show the amount of overplotting for collapsed items.
      This implies that \code{collapse==TRUE}} 

    \item{}{\code{size=1}: Numeric scalar. The relative size of the
      track. Can be overridden in the \code{\link{plotTracks}}
      function.} 

    \item{}{\code{mergeGroups=FALSE}: Logical scalar. Merge fully
      overlapping groups if \code{collapse==TRUE}.}

  } 

  \code{DetailsAnnotationTrack} adds the following additional display parameters:

  \describe{
    \item{}{\code{details.size=0.5}: Numeric scalar. The fraction of
      vertical space of the track used for the details section.}
    
    \item{}{\code{details.minWidth=100}: Numeric scalar. The minium
      width in pixels for a details panel, if less space is available no
      details are plotted.}
    
    \item{}{\code{detailsConnector.col="darkgray"}: Character or integer
      scalar. Color of the line connecting the \code{AnnotstionTrack}
      item with its details panel.}
  
    \item{}{\code{detailsConnector.lty="dashed"}: Character or integer
      scalar. Type of connecting line.}
    
    \item{}{\code{detailsConnector.lwd=1}: Integer scalar. Line width of
      the connector.}
    
    \item{}{\code{detailsConnector.pch=20}: Integer scalar. Type of the
      connector's ends.}
    
    \item{}{\code{detailsConnector.cex=1}: Numeric scalar. Relative size
      of the connector's end points.}
    
    \item{}{\code{detailsBorder.lty="solid"}: Character or integer
      scalar. Line type of the border around each details panel.}
    
    \item{}{\code{detailsBorder.lwd=1}: Integer scalar. Line width of
      the border.}
    
    \item{}{\code{detailsBorder.col="darkgray"}: Character or integer
      scalar. Line color of the border.}
    
    \item{}{\code{detailsBorder.fill="transparent"}: Character or
      integer scalar. Background color of the border.}
    
    \item{}{\code{details.ratio=Inf}: Numeric scalar. By default, the
      plotting method tries to fill all available space of the details
      panel tiles. Depending on the dimensions of your plot and the
      number of tiles this may lead to fairly stretched
      plots. Restricting the ration of width over height can help to
      fine tune for somewhat more sane graphics in these
      cases. Essentially this adds some white space in between
      individual tiles to force the desired ratio. Together with the
      \code{size} and \code{details.size} arguments, which control the
      vertical extension of the whole track and of the details section,
      this allows for some fairly generic resizing of the tiles.}
  
    \item{}{\code{detailsFunArgs=list()}: List. Additional arguments that
      get passed on the the details plotting function.}

    \item{}{\code{groupDetails=FALSE}: Logial scalar. Plot details for
      feature groups rather than for individual features.}
    
  }
  
  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{AnnotationTrack} \code{DetailsAnnotationTrack} objects. 

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

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

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
          size in order to avoid rendering issues. See
          \code{\link{collapsing}} for details. For feathered bars
          indicating the strandedness of grouped items this also
          controls the height of the arrow feathers.}

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

\author{Florian Hahne, Arne Mueller}

\seealso{

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
## An empty object
AnnotationTrack()

## Construct from individual arguments
st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
str <- c("-", "+", "-", "-")
gr <- c("Group1","Group2","Group1", "Group3")

annTrack <- AnnotationTrack(start=st, end=ed, strand=str, chromosome=7, genome="hg19", feature="test",
                            group=gr, id=paste("annTrack item", 1:4), name="generic annotation", stacking="squish")

## Or from a data.frame
df <- data.frame(start=st, end=ed, strand=str, id=paste("annTrack item", 1:4), feature="test",
                 group=gr)
annTrack <- AnnotationTrack(range=df, genome="hg19", chromosome=7, name="generic annotation", 
                             stacking="squish")

## Or from a GRanges object
gr <- GenomicRanges::GRanges(seqnames="chr7", range=IRanges::IRanges(start=df$start, end=df$end), strand=str)
genome(gr) <- "hg19"
IRanges::elementMetadata(gr) <- df[,-(1:3)]
annTrack <- AnnotationTrack(range=gr, name="generic annotation",
stacking="squish")

## Finally from a GRangesList
grl <- split(gr, values(gr)$group)
AnnotationTrack(grl)

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(annTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
plotTracks(annTrack)

## Track names
names(annTrack)
names(annTrack) <- "foo"
plotTracks(annTrack)

## Subsetting and splitting
subTrack <- subset(annTrack, to=2155000)
length(subTrack)
subTrack[1:2]
split(annTrack, c(1,2,1,2))

## Accessors
start(annTrack)
end(annTrack)
width(annTrack)
position(annTrack)
width(subTrack) <- width(subTrack)+1000

strand(annTrack)
strand(subTrack) <- "-"

chromosome(annTrack)
chromosome(subTrack) <- "chrX"

genome(annTrack)
genome(subTrack) <- "mm9"

range(annTrack)
ranges(annTrack)

## Annotation
identifier(annTrack)
identifier(annTrack, lowest=TRUE)
identifier(subTrack) <- "bar"

feature(annTrack)
feature(subTrack) <- "foo"

values(annTrack)

## Grouping
group(annTrack)
group(subTrack) <- "Group 1"
chromosome(subTrack) <- "chr7"
plotTracks(subTrack)

## Stacking
stacking(annTrack)
stacking(annTrack) <- "dense"
plotTracks(annTrack)

## coercion
as(annTrack, "data.frame")
as(annTrack, "UCSCData")

## HTML image map
coords(annTrack)
tags(annTrack)
annTrack <- plotTracks(annTrack)$foo
coords(annTrack)
tags(annTrack)

## DetailsAnnotationTrack
library(lattice) # need to use grid grapics

## generate two random distributions per row (probe/feature)
## the difference between the distributions increases from probe 1 to 4
m <- matrix(c(rgamma(400, 1)), ncol=100)
m[,51:100] <- m[,51:100] + 0:3
## rownames must be accessible by AnnotationTrack element identifier
rownames(m) <- identifier(annTrack, lowest=TRUE)

## create a lattice density plot for the values (signals) of the two groups
## as the chart must be placed into a pre-set grid view port we have to use
## print without calling plot.new! Note, use a common prefix for all lattice.
## Avoid wasting space by removing y-axis decorations.

## Note, in this example 'm' will be found in the environment the 'details'
## function is defined in. To avoid overwriting 'm' you should use a closure
## or environment to access 'm'. 
details <- function(identifier, ...) {
	d = data.frame(signal=m[identifier,], group=rep(c("grp1","grp2"), each=50))
    print(densityplot(~signal, group=group, data=d, main=identifier,
	scales=list(draw=FALSE, x=list(draw=TRUE)), ylab="", xlab="",
	), newpage=FALSE, prefix="plot")
}

deTrack <- AnnotationTrack(range=gr, genome="hg19", chromosome=7,
   name="generic annotation with details per entry", stacking="squish",
   fun=details, details.ratio=1)

plotTracks(deTrack)

set.seed(1234)
deTrack <- AnnotationTrack(range=gr, genome="hg19", chromosome=7,
   name="generic annotation with details per entry",
   stacking="squish",fun=details,
    details.ratio=1, selectFun=function(...){sample(c(FALSE, TRUE), 1)})

plotTracks(deTrack)

}


\keyword{classes}
































