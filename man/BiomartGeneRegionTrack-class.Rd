\name{BiomartGeneRegionTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{BiomartGeneRegionTrack-class}
\alias{BiomartGeneRegionTrack}
\alias{initialize,BiomartGeneRegionTrack-method}

\title{BiomartGeneRegionTrack class and methods}

\description{
  
  A class to hold gene model data for a genomic region fetched
  dynamically from EBI's Biomart Ensembl data source.
  
}

\section{Objects from the class}{

  Objects can be created using the constructor function
  \code{BiomartGeneRegionTrack}.

}

\usage{

BiomartGeneRegionTrack(start, end, biomart, chromosome, strand="*",
genome=NULL, stacking="squish", filters=list(), name="BiomartGeneRegionTrack", ...)

}

\arguments{
  
  \item{start}{An integer scalar with the genomic start coordinates for
    the gene model range.}
  
  \item{end}{An integer scalar with the genomic end coordinates for
    the gene model range.}

  \item{biomart}{An optional \code{\linkS4class{Mart}} object
    providing access to the EBI Biomart webservice. As default the
    appropriate Ensembl data source is selected based on the provided
    genome and chromosome.}
  
  \item{strand}{Character scalar, the strand for which to fetch gene
    information from Biomart. One in \code{+}, \code{-}, or \code{+-}.}
  
  \item{chromosome}{The chromosome on which the track's genomic ranges
    are defined. A valid UCSC chromosome identifier. Please note that at
    this stage only syntactic checking takes place, i.e., the argument
    value needs to be a single integer, numeric character or a character
    of the form \code{chrx}, where \code{x} may be any possible
    string. The user has to make sure that the respective chromosome is
    indeed defined for the the track's genome.}
  
  \item{genome}{The genome on which the track's ranges are
    defined. Usually this is a valid UCSC genome identifier, however
    this is not being formally checked at this point. If no mapping from
    genome to Biomart Ensembl data source is possible, the
    \code{biomart} argument needs to be provided by the user.}
  
  \item{stacking}{The stacking type for overlapping items of the
    track. One in \code{c(hide, dense, squish, pack,full)}. Currently,
    only hide (don't show the track items, squish (make best use of the
    available space) and dense (no stacking at all) are implemented.}
  
  \item{filters}{A list of additional filters to be applied in the
    Biomart query. See \code{\link{getBM}} for details.}
  
  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}
  
  \item{\dots}{Additional items which will all be interpreted as further
    display parameters. See \code{\link{settings}} and the "Display
    Parameters" section below for details.}
  
}

\details{
  
  A track containing all gene models in a particular region as fetched
  from EBI's Biomart service. Usually the user does not have to take
  care of the Biomart connection, which will be established
  automatically based on the provided genome and chromosome
  information. However, for full flexibility a valid
  \code{\linkS4class{Mart}} object may be passed on to the
  constructor. Please note that this assumes a connection to one of the
  Ensembl gene data sources, mapping the available query data back to
  the internal object slots.

}

\value{

  The return value of the constructor function is a new object of class
  \code{BiomartGeneRegionTrack}.
  
}

\section{Slots}{
  \describe{
    \item{\code{biomart}:}{Object of class \code{"MartOrNULL"}, the
      connection to the Ensembl Biomart webservice. }
    
    \item{\code{filter}:}{Object of class \code{"list"}, additional
      filters for the data base query. }
    
    \item{\code{start}:}{Object of class \code{"numeric"}, inherited
      from class \code{\linkS4class{GeneRegionTrack}}. The start
      coordinates of the annotation range. The coorrdinates for the
      individual gene model items are stored in the \code{range} slot.}
    
    \item{\code{end}:}{Object of class \code{"numeric"}, inherited from
      class \code{\linkS4class{GeneRegionTrack}}. The end
      coordinates of the annotation range. The corrdinates for the
      individual gene model items are stored in the \code{range} slot.}

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
      \code{\linkS4class{GdObject}}}

    \item{\code{name}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{GdObject}}}
    
     \item{\code{imageMap}:}{Object of class
      \code{\linkS4class{ImageMap}}, inherited from class
      \code{\linkS4class{GdObject}}}
    
  }
}

\section{Extends}{
  
  Class \code{"\linkS4class{GeneRegionTrack}"}, directly.
  
  Class \code{"\linkS4class{AnnotationTrack}"}, by class "GeneRegionTrack", distance 2.

  Class \code{"\linkS4class{StackedTrack}"}, by class "GeneRegionTrack", distance 3.
  
  Class \code{"\linkS4class{RangeTrack}"}, by class "GeneRegionTrack", distance 4.
  
  Class \code{"\linkS4class{GdObject}"}, by class "GeneRegionTrack", distance 5.
  
}

\section{Methods}{

  In the following code chunks, \code{obj} is considered to be an
  object of class \code{BiomartGeneRegionTrack}.
   
  \describe{

    \bold{\emph{Internal methods:}}
    
    \item{initialize}{\code{signature(.Object =
    "BiomartGeneRegionTrack")}: initialize the object. }
    
  }
  
  \bold{\emph{Inherited methods:}}

  \describe{

     \item{group}{\code{signature(gdObject="BiomartGeneRegionTrack")}: extract
      the group membership for all track items.
		
      \emph{Usage:}
        
      \code{group(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{group(obj)}}

      }
    }
		
    \item{group<-}{\code{signature(gdObject="BiomartGeneRegionTrack",
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
	    
    \item{identifier}{\code{signature(gdObject="BiomartGeneRegionTrack")}:
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

	\item{}{\code{identifier(obj, lowest=FALSE)}}

      }
    }
	    
    \item{identifier<-}{\code{signature(gdObject="BiomartGeneRegionTrack",
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

    \item{exon}{\code{signature(GdObject="BiomartGeneRegionTrack")}: Extract
      the exon identifiers for all exons in the gene models.
      
      \emph{Usage:}
        
      \code{exon(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{exon(obj)}}

      }
    }
    
    \item{exon<-}{\code{signature(GdObject="BiomartGeneRegionTrack",
      value="character")}: replace the exon identifiers for all exons in
      the gene model.  The replacement value must be a character of
      appropriate length or another vector that can be coerced into
      such.
      
      \emph{Usage:}
      
      \code{exon<-(GdObject, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{exon(obj) <- paste("Exon", 1:5)}}

      }
    }
      
    \item{gene}{\code{signature(GdObject="BiomartGeneRegionTrack")}: Extract
      the gene identifiers for all gene models.
      
      \emph{Usage:}
        
      \code{gene(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{gene(obj)}}

      }
    }
    
    \item{gene<-}{\code{signature(GdObject="BiomartGeneRegionTrack",
	value="character")}: replace the gene identifiers for all gene
      models.  The replacement value must be a character of
      appropriate length or another vector that can be coerced into
      such.
      
      \emph{Usage:}
      
      \code{gene<-(GdObject, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{gene(obj) <- paste("Gene", LETTERS[1:5])}}

      }
    }

    \item{symbol}{\code{signature(GdObject="BiomartGeneRegionTrack")}: Extract
      the human-readble gene symbol for all gene models.
	
      \emph{Usage:}
        
      \code{symbol(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{symbol(obj)}}

      }
    }
      
    \item{symbol<-}{\code{signature(GdObject="BiomartGeneRegionTrack",
	value="character")}: replace the human-readable gene symbol for
      all gene models.  The replacement value must be a character of
      appropriate length or another vector that can be coerced into
      such.

      \emph{Usage:}
      
      \code{gene<-(GdObject, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{symbol(obj) <- letters[1:5]}}

      }
    }

    \item{transcript}{\code{signature(GdObject="BiomartGeneRegionTrack")}:
      Extract the transcript identifiers for all transcripts in the
      gene models.
	
      \emph{Usage:}
        
      \code{transcript(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{transcript(obj)}}

      }
    }
    
    \item{transcript<-}{\code{signature(GdObject="BiomartGeneRegionTrack",
	value="character")}: replace the transcript identifiers for
      all transcripts in the gene model.  The replacement value must
      be a character of appropriate length or another vector that
      can be coerced into such.
      
      \emph{Usage:}
      
      \code{transcript<-(GdObject, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{transcript(obj) <- paste("Exon", 1:5)}}
	
      }
    }
  }
	
  \emph{Internal methods:}

  \describe{
	    
    \item{coerce}{\code{signature(from="BiomartGeneRegionTrack",
	to="UCSCData")}: coerce to a \code{UCSCData} object for export to
      the UCSC genome browser.
  
       \emph{Examples:}
      
      \describe{
	
	\item{}{\code{as(obj, "UCSCData")}}

      }
    }

    \item{collapseTrack}{\code{signature(GdObject="BiomartGeneRegionTrack")}:
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
    
    \item{show}{\code{signature(object="BiomartGeneRegionTrack")}: show a
      human-readable summary of the object }
    
    \item{drawGD}{\code{signature(GdObject="BiomartGeneRegionTrack")}: plot the
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

    \item{drawGrid}{\code{signature(GdObject="BiomartGeneRegionTrack")}: superpose a grid on top of a track.

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

    \item{setStacks}{\code{signature(GdObject="BiomartGeneRegionTrack")}:
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
	  
    \item{stacking}{\code{signature(GdObject="BiomartGeneRegionTrack")}: return
      the current stacking type.
    	
      \emph{Usage:}
      
      \code{stacking(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{stacking(obj)}}

      }
    }
    	
    \item{stacking<-}{\code{signature(GdObject="BiomartGeneRegionTrack",
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

    \item{stacks}{\code{signature(GdObject="BiomartGeneRegionTrack")}: return
      the stack indices for each track item.
      
      \emph{Usage:}
      
      \code{stacks(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::stacks(obj)}}
	
      }
    }
    
    \item{[}{\code{signature(x="BiomartGeneRegionTrack", i="ANY", j="ANY",
	drop="ANY")}: subset the items in the \code{BiomartGeneRegionTrack}
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
      
      
    \item{chromosome}{\code{signature(GdObject="BiomartGeneRegionTrack")}:
      return the chromosome for which the track is defined.
      
      \emph{Usage:}
      
      \code{chromosome(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj)}}
	
      }
    }
    
    \item{chromosome<-}{\code{signature(GdObject="BiomartGeneRegionTrack")}:
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
    
    \item{start, end, width}{\code{signature(x="BiomartGeneRegionTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="BiomartGeneRegionTrack")}:
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
    
    \item{position}{\code{signature(GdObject="BiomartGeneRegionTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{feature}{\code{signature(GdObject="BiomartGeneRegionTrack")}: return the
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
    
    \item{feature<-}{\code{signature(gdObject="BiomartGeneRegionTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{BiomartGeneRegionTrack}. See \code{\link{grouping}}
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
    
    \item{genome}{\code{signature(x="BiomartGeneRegionTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="BiomartGeneRegionTrack")}: set the track's
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

    \item{length}{\code{signature(x="BiomartGeneRegionTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="BiomartGeneRegionTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="BiomartGeneRegionTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{split}{\code{signature(x="BiomartGeneRegionTrack")}: split a
      \code{BiomartGeneRegionTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of objects of the same class as the input
      object, all inheriting from class \code{BiomartGeneRegionTrack}.
      
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
    
    \item{strand}{\code{signature(x="BiomartGeneRegionTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="BiomartGeneRegionTrack")}: replace the
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
    
    \item{values}{\code{signature(x="BiomartGeneRegionTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="BiomartGeneRegionTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{subset}{\code{signature(x="BiomartGeneRegionTrack")}: subset a
      \code{BiomartGeneRegionTrack} by coordinates and sort if necessary.
      
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
    
    \item{displayPars}{\code{signature(x="BiomartGeneRegionTrack",
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
      
    
    \item{displayPars}{\code{signature(x="BiomartGeneRegionTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="BiomartGeneRegionTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="BiomartGeneRegionTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="BiomartGeneRegionTrack",
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
    
    \item{setPar}{\code{signature(x="BiomartGeneRegionTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{BiomartGeneRegionTrack} class are
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
    
    \item{setPar}{\code{signature(x="BiomartGeneRegionTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{BiomartGeneRegionTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{names}{\code{signature(x="BiomartGeneRegionTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="BiomartGeneRegionTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="BiomartGeneRegionTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="BiomartGeneRegionTrack")}: return the tags from the
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
  \code{BiomartGeneRegionTrack} upon instantiation, unless one or more
  of them have already been set by one of the optional sub-class
  initializers, which always get precedence over these global defaults.
  See \code{\link{settings}} for details on setting graphical parameters
  for tracks. \describe{ 

    \item{}{\code{C_segment="burlywood4"}: Character or integer scalar.
      Fill color for annotation objects of type 'C_segment'.} 

    \item{}{\code{D_segment="lightblue"}: Character or integer scalar.
      Fill color for annotation objects of type 'C_segment'.} 

    \item{}{\code{J_segment="dodgerblue2"}: Character or integer scalar.
      Fill color for annotation objects of type 'C_segment'.} 

    \item{}{\code{miRNA="cornflowerblue"}: Character or integer scalar.
      Fill color for annotation objects of type 'L_segment'.} 

    \item{}{\code{miRNA_pseudogene="cornsilk"}: Character or integer
      scalar. Fill color for annotation objects of type 'miRNA_pseudogene'.} 

    \item{}{\code{misc_RNA="cornsilk3"}: Character or integer scalar.
      Fill color for annotation objects of type 'misc_RNA'.} 

    \item{}{\code{misc_RNA_pseudogene="cornsilk4"}: Character or integer
      scalar. Fill color for annotation objects of type 'misc_RNA_pseudogene'.} 

    \item{}{\code{Mt_rRNA="yellow"}: Character or integer scalar. Fill
      color for annotation objects of type 'Mt_rRNA'.} 

    \item{}{\code{Mt_tRNA="darkgoldenrod"}: Character or integer scalar.
      Fill color for annotation objects of type 'Mt_tRNA'.} 

    \item{}{\code{Mt_tRNA_pseudogene="darkgoldenrod1"}: Character or
      integer scalar. Fill color for annotation objects of type
      'Mt_tRNA_pseudogene'.} 

    \item{}{\code{protein_coding="gold4"}: Character or integer scalar.
      Fill color for annotation objects of type 'protein_coding'.} 

    \item{}{\code{pseudogene="brown1"}: Character or integer scalar.
      Fill color for annotation objects of type 'pseudogene'.} 

    \item{}{\code{retrotransposed="blueviolet"}: Character or integer
      scalar. Fill color for annotation objects of type 'retrotransposed'.} 

    \item{}{\code{rRNA="darkolivegreen1"}: Character or integer scalar.
      Fill color for annotation objects of type 'rRNA'.} 

    \item{}{\code{rRNA_pseudogene="darkolivegreen"}: Character or
      integer scalar. Fill color for annotation objects of type
      'rRNA_pseudogene'.} 

    \item{}{\code{scRNA="darkorange"}: Character or integer scalar.
      Fill color for annotation objects of type 'scRNA'.} 

    \item{}{\code{scRNA_pseudogene="darkorange2"}: Character or integer
      scalar. Fill color for annotation objects of type 'scRNA_pseudogene'.} 

    \item{}{\code{snoRNA="cyan"}: Character or integer scalar. Fill
      color for annotation objects of type 'snoRNA'.} 

    \item{}{\code{snoRNA_pseudogene="cyan2"}: Character or integer
      scalar. Fill color for annotation objects of type 'snoRNA_pseudogene'.} 

    \item{}{\code{snRNA="coral"}: Character or integer scalar. Fill
      color for annotation objects of type 'snRNA'.} 

    \item{}{\code{snRNA_pseudogene="coral3"}: Character or integer
      scalar. Fill color for annotation objects of type 'snRNA_pseudogene'.} 

    \item{}{\code{tRNA_pseudogene="antiquewhite3"}: Character or integer
      scalar. Fill color for annotation objects of type 'tRNA_pseudogene'.} 

    \item{}{\code{V_segment="aquamarine"}: Character or integer scalar.
      Fill color for annotation objects of type 'V_segment'.} 

  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{BiomartGeneRegionTrack} objects. 

    \describe{ 

    \item{}{\code{\linkS4class{GeneRegionTrack}}: 

      \describe{ 

        \item{}{\code{fill="orange"}: Character or integer scalar. The
          fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{geneSymbols=TRUE}: Logical scalar. Use human-readable
          gene symbols or gene IDs for the transcript annotation.} 

        \item{}{\code{shape=c("smallArrow", "box")}: Character scalar.
          The shape in which to display the track items. Currently only
          \code{box}, \code{arrow}, \code{ellipse}, and \code{smallArrow}
          are implemented.} 

        \item{}{\code{showExonId=FALSE}: Logical scalar. Control whether
          to plot the individual exon identifiers.}

	\item{}{\code{mergeGroups=FALSE}: Logical scalar. Merge fully
	  overlapping groups if \code{collapse==TRUE}.}

      } 

    } 

    \item{}{\code{\linkS4class{AnnotationTrack}}: 

      \describe{ 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{col="transparent"}: Character or integer scalar.
          The border color for all track items.} 

        \item{}{\code{fontcolor="white"}: Character or integer scalar.
          The font color for item identifiers.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{fontface=1}: Integer scalar. The font face for
          item identifiers.} 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{fontfamily="sans"}: Character scalar. The font
          family for item identifiers.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          item identifiers.} 

        \item{}{\code{lex=1}: Numeric scalar. The line expansion factor
          for all track items. This is also used to connect grouped
          items. See \code{\link{grouping}} for details.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for item identifiers.} 

        \item{}{\code{lty="solid"}: Character or integer scalar. The
          line type for all track items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{lwd=1}: Integer scalar. The line width for all
          track items. This is also used to connect grouped items. See
          \code{\link{grouping}} for details.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

        \item{}{\code{showFeatureId=FALSE}: Logical scalar. Control
          whether to plot the individual track item identifiers.} 

        \item{}{\code{showId=FALSE}: Logical scalar. Control whether
          to annotate individual groups.} 

        \item{}{\code{showOverplotting=FALSE}: Logical scalar. Use a
          color gradient to show the amount of overplotting for collapsed
          items. This implies that \code{collapse==TRUE}} 

        \item{}{\code{size=1}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function.}

	\item{}{\code{mergeGroups=FALSE}: Logical scalar. Merge fully
	  overlapping groups if \code{collapse==TRUE}.}

      } 

    }

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

\references{
  
  EBI Biomart webservice at \url{http://www.biomart.org}.

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

  \code{\linkS4class{Mart}}

  \code{\linkS4class{RangeTrack}}

  \code{\linkS4class{StackedTrack}}

  \code{\link{collapsing}}

  \code{\link{DataTrack}}

  \code{\link{getBM}}

  \code{\link{grouping}}

  \code{\link{panel.grid}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

  \code{\link{useMart}}

}

\examples{

\dontshow{
## Load some sample data
data(bmTrack)
}

## Construct the object
\dontrun{
bmTrack <- BiomartGeneRegionTrack(start=26682683, end=26711643,
chromosome=7, genome="mm9")
}

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(bmTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
plotTracks(bmTrack)

## Track names
names(bmTrack)
names(bmTrack) <- "foo"
plotTracks(bmTrack)

## Subsetting and splitting
subTrack <- subset(bmTrack, from=26700000, to=26705000)
length(subTrack)
subTrack <- bmTrack[transcript(bmTrack)=="ENSMUST00000144140"]
split(bmTrack, transcript(bmTrack))

## Accessors
start(bmTrack)
end(bmTrack)
width(bmTrack)
position(bmTrack)
width(subTrack) <- width(subTrack)+100

strand(bmTrack)
strand(subTrack) <- "-"

chromosome(bmTrack)
chromosome(subTrack) <- "chrX"

genome(bmTrack)
genome(subTrack) <- "hg19"

range(bmTrack)
ranges(bmTrack)

## Annotation
identifier(bmTrack)
identifier(bmTrack, lowest=TRUE)
identifier(subTrack) <- "bar"

feature(bmTrack)
feature(subTrack) <- "foo"

exon(bmTrack)
exon(subTrack) <- letters[1:2]

gene(bmTrack)
gene(subTrack) <- "bar"

symbol(bmTrack)
symbol(subTrack) <- "foo"

transcript(bmTrack)
transcript(subTrack) <- c("foo", "bar")
chromosome(subTrack) <- "chr7"
plotTracks(subTrack)

values(bmTrack)

## Grouping
group(bmTrack)
group(subTrack) <- "Group 1"
transcript(subTrack)
plotTracks(subTrack)

## Stacking
stacking(bmTrack)
stacking(bmTrack) <- "dense"
plotTracks(bmTrack)

## coercion
as(bmTrack, "data.frame")
as(bmTrack, "UCSCData")

## HTML image map
coords(bmTrack)
tags(bmTrack)
bmTrack <- plotTracks(bmTrack)$foo
coords(bmTrack)
tags(bmTrack)

}

\keyword{classes}











































