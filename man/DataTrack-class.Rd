\name{DataTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{DataTrack-class}
\alias{DataTrack}
\alias{[,DataTrack-method}
\alias{[,DataTrack,ANY,ANY-method}
\alias{collapseTrack,DataTrack-method}
\alias{drawAxis,DataTrack-method}
\alias{drawGD,DataTrack-method}
\alias{initialize,DataTrack-method}
\alias{score,DataTrack-method}
\alias{show,DataTrack-method}
\alias{split,DataTrack-method}
\alias{values,DataTrack-method}
\alias{values<-,DataTrack-method}
\alias{feature,DataTrack-method}
\alias{feature<-,DataTrack,character-method}
\alias{range,DataTrack-method}
\alias{strand,DataTrack-method}
\alias{subset,DataTrack-method}
\alias{strand<-,DataTrack-method}
\alias{coerce,DataTrack,data.frame-method}
\alias{coerce,GRanges,DataTrack-method}

\title{DataTrack class and methods}

\description{

  A class to store numeric data values along genomic
  coordinates. Multiple samples as well as sample groupings are
  supported, with the restriction of equal genomic coordinates for a
  single observation across samples.
  
}

\section{Objects from the class}{
  
  Objects can be created using the constructor function
  \code{DataTrack}.

}

\usage{

DataTrack(range=NULL, start=NULL, end=NULL, width=NULL, data, chromosome, strand="*", genome,
name="DataTrack", ...)

}

\arguments{

  We tried to keep instantiation of \code{DataTrack} objects as flexible
  as possible to accomodate different use cases. For instance, one
  natural way to create a \code{DataTrack} is from an existing
  \code{\linkS4class{GRanges}} object. In other cases it might be more
  appropriate to build the object using individual function arguments.
  
  \item{range}{
    
    An optional meta argument to handle the different input types. If
    the \code{range} argument is missing, all the relevant information
    to create the object has to be provided as individual function
    arguments (see below).

    The different input options for \code{range} are:

    \describe{ 

      \item{}{A \code{GRanges} object: essentially all the necessary
	information to create a \code{DataTrack} can be contained in a
	single \code{GRanges} object. The track's coordinates are taken
	from the \code{start}, \code{end} and \code{seqnames} slots, the
	genome information from the genome slot, and the numeric data
	values can be extracted from additional \code{elementMetadata}
	columns (please note that non-numeric columns are being ignored
	with a warning). As a matter of fact, calling the constructor on
	a \code{GRanges} object without further arguments, e.g.
	\code{DataTrack(range=obj)} is equivalent to calling the coerce
	method \code{as(obj, "DataTrack")}. Alternatively, the
	\code{GRanges} object may only contain the coordinate
	information, in which case the numeric data part is expected to
	be present in the separate \code{data} argument, and the ranges
	have to match the dimensions of the data matrix. If \code{data}
	is not \code{NULL}, this will always take precedence over
	anything defined in the \code{range} argument. See below for
	details.}

      \item{}{An \code{\linkS4class{IRanges}} object: this is very
	similar to the above case, except that the numeric data part now
	always has to be provided in the separate \code{data}
	argument. Also the chromosome information must be provided in
	the \code{chromosome} argument, because neither of the two can
	be directly encoded in an \code{IRange} object.}

      \item{}{A \code{data.frame} object: the \code{data.frame} needs to
	contain at least the two mandatory columns \code{start} and
	\code{end} with the range coordinates. It may also contain a
	\code{chromosome} column with the chromosome information for
	each range. If missing it will be drawn from the separate
	\code{chromosome} argument. All additional numeric columns will
	be interpreted as data columns, unless the \code{data} argument
	is explicitely provided.}
	
    }
  }

  \item{start, end, width}{Integer vectors, giving the start and the end
    end coordinates for the individual track items, or their width. Two
    of the three need to be specified, and have to be of equal length or
    of length one, in which case the single value will be recycled
    accordingly. Otherwise, the usual R recycling rules for vectors do
    not apply and the function will cast an error.}
  
  \item{data}{A numeric matrix of data points with the number of columns
    equal to the number of coordinates in \code{range}, or a numeric
    vector of appropriate length that will be coerced into such a
    one-row matrix. Each individual row is supposed to contain data for
    a given sample, where the coordinates for each single observation
    are constant across samples. Depending on the plotting type of the
    data (see 'Details' and 'Display Parameters' sections), sample
    grouping or data aggregation may be available. Alternatively, this
    can be a character vector of column names that point into the
    element metadata of the \code{range} object for
    subsetting. Naturally, this is only supported when the \code{range}
    argument is of class \code{GRanges}.}
  
  \item{strand}{Character vector, the strand information for the
    individual track items. Currently this has to be unique for the
    whole track and doesn't really have any visible consequences, but we
    might decide to make \code{DataTracks} strand-specific at a later
    stage.}
  
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

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

  \item{\dots}{Additional items which will all be interpreted as further
    display parameters.}
  
}

\value{

  The return value of the constructor function is a new object of class
  \code{DataTrack}.
  
}
   
\section{Slots}{
  
  \describe{
  
    \item{\code{data}:}{Object of class \code{"matrix"}, containing the
      data values to be plotted. Individual rows of the matrix
      correspond to individual samples, and the number of columns has to
      be identical to the feature number of the \code{GRanges} object in
      the \code{range} slot.}

    \item{\code{strand}:}{Object of class \code{"character"}, the strand
      information for the track, in the form '+' for the Watson strand,
      '-' for the Crick strand or '*' for either of the two.}

    \item{\code{range}:}{Object of class \code{\linkS4class{IRanges}},
      inherited from class \code{\linkS4class{RangeTrack}}. The genomic
      coordinates for the data values. The length of the object needs to
      be identical to the number of columns of the data matrix in the
      \code{data} slot.}
    
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
  
  Class \code{"\linkS4class{NumericTrack}"}, directly.
  
  Class \code{"\linkS4class{RangeTrack}"}, by class "NumericTrack", distance 2.
  
  Class \code{"\linkS4class{GdObject}"}, by class "NumericTrack", distance 3.
  
}

\section{Methods}{
  
  In the following code chunks, \code{obj} is considered to be an
  object of class \code{DataTrack}.
  
  \bold{\emph{Exported in the name space:}}
  
  \describe{
    
    \item{[}{\code{signature(x="DataTrack")}: subsetting of the
      object, either to a subet of coordinates, or to a subset of
      samples.
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{i}, \code{j}: subsetting indices for coordinates
	  (\code{i}) or samples (\code{j}).}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{obj[1:3,]}}
	
	\item{}{\code{obj[,2:4]}}
	
      }
    }
    
    \item{values}{\code{signature(x="DataTrack")}: return the raw data
      values of the object, i.e., the data matrix in the \code{data}
      slot.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}
	
      }
    }

     \item{values<-}{\code{signature(x="DataTrack")}: replace the data
      matrix in the \code{data} slot.
      
      \emph{Usage:}
      
      \code{values<-(x, value)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{value}: replacement value.}
	  
      }

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj) <- matrix(1:10, ncol=2)}}
	
      }
    }
    
    \item{score}{\code{signature(x="DataTrack")}: return processed data
      values of the object exactly like they would be plotted to the
      device (modulo any potential aggregration or collapsing), i.e.,
      the raw data with optional transformations applied.
      
      \emph{Usage:}
      
      \code{score(x, from=NULL, to=NULL, sort=FALSE, transformation=TRUE)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: restrict to data
	  within a certain coordinates range.}
	
	\item{}{\code{sort}: sort the return values by
	  coordinates. This is usually not necessary since the data
	  should already be ordererd, however this is not formaly
	  checked anywhere and some operations strictly depend on
	  ordered data.}
	
	\item{}{\code{transformation}: apply a data transformation
	  in case one is defined as the \code{transformation} display
	  parameter.}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{score(obj)}}
	
	\item{}{\code{score(obj, from=100, to=10000)}}
	
	\item{}{\code{score(obj, sort=TRUE, transformation=FALSE)}}
	
      }
    }
      
    \item{split}{\code{signature(x="DataTrack")}: split a
      \code{DataTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of \code{DataTrack} objects.

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
    
    \item{range, ranges}{\code{signature(x="DataTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
       \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

	\item{}{\code{ranges(obj)}}

      }
    }

    \item{strand}{\code{signature(x="DataTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="DataTrack")}: replace the
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
    
    \item{feature}{\code{signature(GdObject="DataTrack")}: returns NULL
      since there is no grouping information for the ranges in a
      \code{DataTrack}.
      
      \emph{Usage:}
      
      \code{feature(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{feature(obj)}}

      }
    }
    
    \item{feature<-}{\code{signature(gdObject="DataTrack",
      value="character")}: this return the unaltered input object since
      there  is no grouping information for the ranges in a
      \code{DataTrack}.
      
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
  }
	
  \bold{\emph{Internal methods:}}

  \describe{
  
    \item{collapseTrack}{\code{signature(gdObject="DataTrack")}:
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
    
    \item{drawGD}{\code{signature(GdObject="DataTrack")}: plot the
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
    
    \item{drawAxis}{\code{signature(GdObject="DataTrack")}: add a y-axis
      to the title panel of a track.

      \emph{Usage:}
      
      \code{drawAxis(GdObject, from, to, ...)}
      
      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: compute axis range from the data
	  within a certain coordinates range only.}
	
	\item{}{\code{\dots}: all further arguments are ignored.}
	
      }

      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::drawAxis(obj)}}
	
      }
    }
    
    \item{initialize}{\code{signature(.Object="DataTrack")}:
      initialize the object }
    
    \item{show}{\code{signature(object="DataTrack")}: show a
      human-readable summary of the object }
    
  }
  
  \bold{\emph{Inherited methods:}}
  
  \describe{
    
    \item{drawGrid}{\code{signature(GdObject="DataTrack")}:
      superpose a grid on top of a track.
    
      \emph{Usage:}
      
      \code{drawGrid(GdObject, from, to, ...)}

      
      \emph{Additional Arguments:}
      
      \describe{

	\item{}{\code{from}, \code{to}: integer scalars, restrict to
	  coordinate range before computing the grid lines.}

      }

      \emph{Examples:}
      
      \describe{

	\item{}{\code{Gviz:::drawGrid(obj)}}

      }
    }

    \item{chromosome}{\code{signature(GdObject="DataTrack")}: return the
      currently active chromosome for which the track is defined. For
      consistancy with other Bioconductor packages, the
      \code{isActiveSeq} alias is also provided.
      
      \emph{Usage:}
      
      \code{chromosome(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj)}}
	
      }
    }
    
    \item{chromosome<-}{\code{signature(GdObject="DataTrack")}: replace
      the value of the track's active chromosome. This has to be a valid
      UCSC chromosome identifier or an integer or character scalar that
      can be reasonably coerced into one, unless
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
    
    \item{start, end, width}{\code{signature(x="DataTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="DataTrack")}:
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
    
    \item{position}{\code{signature(GdObject="DataTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{genome}{\code{signature(x="DataTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="DataTrack")}: set the track's
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

    \item{length}{\code{signature(x="DataTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="DataTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{subset}{\code{signature(x="DataTrack")}: subset a
      \code{NumericTrack} by coordinates and sort if necessary.
      
      \emph{Usage:}
        
      \code{subset(x, from, to, sort=FALSE, drop=TRUE, ...)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{from}, \code{to}: the coordinates range to subset
	  to.}
	
	\item{}{\code{sort}: sort the object after subsetting. Usually
	  not necessary.}

	\item{}{\code{drop}: drop unused regions on the other,
	  non-active chromosomes.w}
	
	\item{}{\code{\dots}: additional arguments are ignored.}
	
      }
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{subset(obj, from=10, to=20, sort=TRUE)}}

      }
    }
    
    \item{displayPars}{\code{signature(x="DataTrack",
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
      
    
    \item{displayPars}{\code{signature(x="DataTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="DataTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="DataTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="DataTrack",
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
    
    \item{setPar}{\code{signature(x="DataTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{DataTrack} class are
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
    
    \item{setPar}{\code{signature(x="DataTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{DataTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="DataTrack")}: return
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
    
    \item{names}{\code{signature(x="DataTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="DataTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="DataTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="DataTrack")}: return the tags from the
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
  \code{DataTrack} upon instantiation, unless one or more of them have
  already been set by one of the optional sub-class initializers, which
  always get precedence over these global defaults. See \code{\link{settings}}
  for details on setting graphical parameters for tracks. \describe{ 

    \item{}{\code{aggregation="mean"}: Function or character scalar.
      Used to aggregate values in windows or for collapsing overlapping
      items. The function has to accept a numeric vector as a single
      input parameter and has to return a numeric scalar with the
      aggregated value. Alternatively, one of the predefined options
      \code{mean}, \code{median} \code{sum}, \code{min}, \code{max} or
      \code{extreme} can be supplied as a character scalar. Defaults
      to \code{mean}.}
    
     \item{}{\code{aggregateGroups=FALSE}:Logical scalar. Aggregate the
       values within a sample group using the aggregation function
       specified in the \code{aggregate} parameter.}

    \item{}{\code{amount=NULL}: Numeric scalar. Amount of jittering in
      xy-type plots. See \code{\link{panel.xyplot}} for details.} 

    \item{}{\code{baseline=NULL}: Numeric scalar. Y-axis position of
      an optional baseline. This parameter has a special meaning for
      mountain-type plots, see the 'Details' section in
      \code{\linkS4class{DataTrack}} for more information.} 

    \item{}{\code{box.ratio=1}: Numeric scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{box.width=NULL}: Numeric scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{cex=0.7}: Numeric scalar. The default pixel size for
      plotting symbols.} 

    \item{}{\code{coef=1.5}: Numeric scalar. Parameter controlling the
      boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{col.baseline=NULL}: Character scalar. Color for the
      optional baseline, defaults to the setting of \code{col}.} 

    \item{}{\code{col.histogram="#808080"}: Character scalar. Line
      color in histogram-type plots.} 

    \item{}{\code{col.mountain=NULL}: Character scalar. Line color in
      mountain-type plots, defaults to the setting of \code{col}.} 

    \item{}{\code{collapse=FALSE}: Logical scalar. Collapse overlapping
      ranges and aggregate the underlying data.} 

    \item{}{\code{degree=1}: Numeric scalar. Parameter controlling the
      loess calculation for smooth and mountain-type plots. See
      \code{\link{panel.loess}} for details.} 

    \item{}{\code{do.out=TRUE}: Logical scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{evaluation=50}: Numeric scalar. Parameter controlling
      the loess calculation for smooth and mountain-type plots. See
      \code{\link{panel.loess}} for details.} 

    \item{}{\code{factor=0.5}: Numeric scalar. Factor to control amount
      of jittering in xy-type plots. See \code{\link{panel.xyplot}}
      for details.} 

    \item{}{\code{family="symmetric"}: Character scalar. Parameter
      controlling the loess calculation for smooth and mountain-type
      plots. See \code{\link{panel.loess}} for details.} 

    \item{}{\code{fill.histogram="lightgray"}: Character scalar. Fill
      color in histogram-type plots, defaults to the setting of
      \code{fill}.} 

    \item{}{\code{fill.mountain=c("#CCFFFF", "#FFCCFF")}: Character
      vector of length 2. Fill color in mountain-type plots.} 

    \item{}{\code{gradient=c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
      "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")}: Character
      vector. The base colors for the 'gradient' plotting type.} 

    \item{}{\code{groups=NULL}: Vector coercable to a factor. Optional
      sample grouping. See 'Details' section in \code{\linkS4class{DataTrack}}
      for further information.} 

    \item{}{\code{jitter.x=FALSE}: Logical scalar. Toggle on jittering
      on the x axis in xy-type plots. See \code{\link{panel.xyplot}}
      for details.} 

    \item{}{\code{jitter.y=FALSE}: Logical scalar. Toggle off jittering
      on the y axis in xy-type plots. See \code{\link{panel.xyplot}}
      for details.} 

    \item{}{\code{levels.fos=NULL}: Numeric scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{lty.baseline=NULL}: Character or numeric scalar. Line
      type of the optional baseline, defaults to the setting of
      \code{lty}.} 

    \item{}{\code{lty.mountain=NULL}: Character or numeric scalar. Line
      type in mountain-type plots, defaults to the setting of \code{lty}.} 

    \item{}{\code{lwd.baseline=NULL}: Numeric scalar. Line width of
      the optional baseline, defaults to the setting of \code{lwd}.} 

    \item{}{\code{lwd.mountain=NULL}: Numeric scalar. Line width in
      mountain-type plots, defaults to the setting of \code{lwd}.} 

    \item{}{\code{min.distance=0}: Numeric scalar. The mimimum distance
      in pixel below which to collapse ranges.} 

    \item{}{\code{na.rm=FALSE}: Boolean controlling whether to discard
      all NA values when plotting or to keep empty spaces for NAs} 

    \item{}{\code{ncolor=100}: Integer scalar. The number of colors
      for the 'gradient' plotting type} 

    \item{}{\code{notch=FALSE}: Logical scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{notch.frac=0.5}: Numeric scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{pch=20}: Integer scalar. The type of glyph used for
      plotting symbols.} 

    \item{}{\code{separator=0}: Numeric scalar. Number of pixels used
      to separate individual samples in heatmap-type plots.} 

    \item{}{\code{size=NULL}: Numeric scalar. The relative size of the
      track. Can be overridden in the \code{\link{plotTracks}} function.
      By default the size will be set automatically based on the selected
      plotting type.} 

    \item{}{\code{span=0.2}: Numeric scalar. Parameter controlling the
      loess calculation for smooth and mountain-type plots. See
      \code{\link{panel.loess}} for details.} 

    \item{}{\code{stackedBars=TRUE}: Logical scalar. When there are
      several data groups, draw the histogram-type plots as stacked
      barplots or grouped side by side.} 

    \item{}{\code{stats=X[[44]]}: Function. Parameter controlling the
      boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{transformation=NULL}: Function. Applied to the data
      matrix prior to plotting or when calling the \code{score} method.
      The function should accept exactly one input argument and its
      return value needs to be a numeric vector which can be coerced
      back into a data matrix of identical dimensionality as the input
      data.} 

    \item{}{\code{type="p"}: Character vector. The plot type, one or
      several in \code{c("p","l", "b", "a", "s", "g", "r", "S", "smooth",
      "histogram", "mountain", "h", "boxplot", "gradient", "heatmap")}.
      See 'Details' section in \code{\linkS4class{DataTrack}} for more
      information on the individual plotting types.} 

    \item{}{\code{varwidth=FALSE}: Logical scalar. Parameter controlling
      the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

    \item{}{\code{window=NULL}: Numeric or character scalar. Aggregate
      the rows values of the data matrix to \code{window} equally sized
      slices on the data range using the method defined in \code{aggregation}.
      If negative, apply a running window of size \code{windowSize}
      using the mean aggregation method. Alternatively, the special
      value \code{auto} causes the function to determine the optimal
      window size to avoid overplotting, and \code{fixed} uses
      fixed-size windows of size \code{windowSize}.} 

    \item{}{\code{windowSize=NULL}: Numeric scalar. The size of the
      running window when the value of \code{window} is negative.} 

    \item{}{\code{ylim=NULL}: Numeric vector of length 2. The range of
      the y-axis scale.} 

  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{DataTrack} objects. 

    \describe{ 

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


\details{

  Depending on the setting of the \code{type} display parameter, the
  data can be plotted in various different forms as well as combinations
  thereof. Supported plotting types are:

  \describe{
    
    \item{}{\code{p}: simple xy-plot.}
    
    \item{}{\code{l}: lines plot. In the case of multiple samples this
      plotting type is not overly usefull since the points in the data
      matrix are connected in column-wise order. Type \code{a} might be
      more appropriate in these situations.}

    \item{}{\code{b}: combination of xy-plot and lines plot.}

    \item{}{\code{a}: lines plot of the column-wise average values.}

    \item{}{\code{s}: sort and connect data points along the x-axis}

    \item{}{\code{S}: sort and connect data points along the y-axis}

    \item{}{\code{g}: add grid lines. To ensure a consitant look and
      feel across multiple tracks, grid lines should preferentially be
      added by using the \code{grid} display parameter.}

    \item{}{\code{r}: add a regression line to the plot.}
    
    \item{}{\code{h}: histogram-like vertical lines centered in the
      middle of the coordinate ranges.}
        
    \item{}{\code{smooth}: add a loess fit to the plot. The following
      display parameters can be used to control the loess calculation:
      \code{span, degree, family, evaluation}. See
      \code{\link{panel.loess}} for details.}

    \item{}{\code{histogram}: plot data as a histogram, where the width
      of the histogram bars reflects the width of the genomic ranges in
      the \code{range} slot.}

    \item{}{\code{mountain}: plot a smoothed version of the data
      relative to a baseline, as defined by the \code{baseline} display
      parameter. The following display parameters can be used to control
      the smoothing: \code{span, degree, family, evaluation}. See
      \code{\link{panel.loess}} for details. The layout of the plot can
      be further customized via the following display parameters:
      \code{col.mountain, lwd.mountain, lty.mountain, fill.mountain}.}


    \item{}{\code{boxplot}: plot the data as box-and-whisker
      plots. The layout of the plot can be further customized via the
      following display parameters: \code{box.ratio, box.width, varwidt,
      notch, notch.frac, levels.fos, stats, coef, do.out}. See
      \code{\link{panel.bwplot}} for details.}
    
    
    \item{}{\code{gradient}: collapse the data across samples and plot this
      average value as a color-coded gradient. Essenitally this is
      similar to the heatmap-type plot of a single sample. The layout of
      the plot can be further customized via the display parameters
      \code{ncolor} and \code{gradient} which control the number of
      gradient colors as well as the gradient base colors,
      respectively.}
    
    \item{}{\code{heatmap}: plot the color-coded values for all samples
      in the form of a heatmap. The data for individual samples can be
      visually separated by setting the \code{separator} display
      parameter. It's value is taken as the amount of spacing in pixels
      in between two heatmap rows. The layout of the plot can be further
      customized via the display parameters \code{ncolor} and
      \code{gradient} which control the number of gradient colors as
      well as the gradient base colors, respectively.}
    
  }

  For some of the above plotting-types the \code{groups} display
  parameter can be used to indicate sample sub-groupings. Its value is
  supposed to be a factor vector of similar length as the number of
  samples. In most cases, the groups are shown in different plotting
  colors and data aggregation operations are done in a stratified
  fashion.

  The \code{window} display parameter can be used to aggregate the data
  prior to plotting. Its value is taken as the number of equal-sized
  windows along the genomic coordinates of the track for which to
  compute average values. The special value \code{auto} can be used to
  automatically determine a reasonable number of windows which can be
  particularly useful when plotting very large genomic regions with many
  data points.

  The \code{aggregation} parameter can be set to define the aggregation
  function to be used when averaging in windows or across collapsed
  items. It takes the form of either a function which should condense a
  numeric vector into a single number, or one of the predefined options
  as character scalars \code{"mean"}, \code{"median"} or \code{"sum"}
  for mean, median or summation, respectively. Defaults to computing
  mean values for each sample. Note that the predefined options can be
  much faster because they are optimized to work on large numeric
  tables.
  
}

\author{Florian Hahne}



\seealso{

  \code{\linkS4class{DataTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\linkS4class{GRanges}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{IRanges}}

  \code{\linkS4class{NumericTrack}}

  \code{\linkS4class{RangeTrack}}

  \code{\link{collapsing}}

  \code{\link{grouping}}

  \code{\link{panel.bwplot}}

  \code{\link{panel.grid}}

  \code{\link{panel.loess}}

  \code{\link{panel.xyplot}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}

\examples{

## Object construction:

## An empty object
DataTrack()

## from individual arguments
dat <- matrix(runif(400), nrow=4)
dtTrack <- DataTrack(start=seq(1,1000, len=100), width=10,
data=dat, chromosome=1, genome="mm9",
name="random data")

## from GRanges
library(GenomicRanges)
gr <- GRanges(seqnames="chr1", ranges=IRanges(seq(1,1000, len=100),
width=10))
values(gr) <- t(dat)
dtTrack <- DataTrack(range=gr, genome="mm9", name="random data")

## from IRanges
dtTrack <- DataTrack(range=ranges(gr), data=dat, genome="mm9",
name="random data", chromosome=1)

## from a data.frame
df <- as.data.frame(gr)
colnames(df)[1] <- "chromosome"
dtTrack <- DataTrack(range=df,  genome="mm9", name="random data")

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(dtTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
plotTracks(dtTrack)

## Track names
names(dtTrack)
names(dtTrack) <- "foo"
plotTracks(dtTrack)

## Subsetting and splitting
subTrack <- subset(dtTrack, from=100, to=300)
length(subTrack)
subTrack[1:2,]
subTrack[,1:2]
split(dtTrack, rep(1:2, each=50))

## Accessors
start(dtTrack)
end(dtTrack)
width(dtTrack)
position(dtTrack)
width(subTrack) <- width(subTrack)-5

strand(dtTrack)
strand(subTrack) <- "-"

chromosome(dtTrack)
chromosome(subTrack) <- "chrX"

genome(dtTrack)
genome(subTrack) <- "mm9"

range(dtTrack)
ranges(dtTrack)

## Data
values(dtTrack)
score(dtTrack)

## coercion
as(dtTrack, "data.frame")

}

\keyword{classes}






























