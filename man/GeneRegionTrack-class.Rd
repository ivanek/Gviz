\name{GeneRegionTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{GeneRegionTrack-class}
\alias{GeneRegionTrack}
\alias{coerce,GeneRegionTrack,UCSCData-method}
\alias{collapseTrack,GeneRegionTrack-method}
\alias{exon}
\alias{exon,GeneRegionTrack-method}
\alias{exon<-}
\alias{exon<-,GeneRegionTrack,character-method}
\alias{gene,GeneRegionTrack-method}
\alias{gene}
\alias{gene<-,GeneRegionTrack,character-method}
\alias{gene<-}
\alias{group,GeneRegionTrack-method}
\alias{group<-,GeneRegionTrack,character-method}
\alias{identifier,GeneRegionTrack-method}
\alias{identifier<-,GeneRegionTrack,character-method}
\alias{initialize,GeneRegionTrack-method}
\alias{show,GeneRegionTrack-method}
\alias{symbol,GeneRegionTrack-method}
\alias{symbol<-,GeneRegionTrack,character-method}
\alias{symbol}
\alias{symbol<-}
\alias{transcript,GeneRegionTrack-method}
\alias{transcript<-,GeneRegionTrack,character-method}
\alias{transcript}
\alias{transcript<-}
\alias{coerce,GRanges,GeneRegionTrack-method}
\alias{coerce,GRangesList,GeneRegionTrack-method}
\alias{coerce,TranscriptDb,GeneRegionTrack-method}

\title{GeneRegionTrack class and methods}

\description{
  
  A class to hold gene model data for a genomic region.

}

\section{Objects from the class}{
  
  Objects can be created using the constructor function
  \code{GeneRegionTrack}.

}


\usage{

GeneRegionTrack(range=NULL, rstarts=NULL, rends=NULL, rwidths=NULL, strand="*",
               feature="unknown", exon=feature, transcript=feature, gene=feature,
               symbol=feature, chromosome=NULL, genome, stacking="squish",
               name="GeneRegionTrack", start=NULL, end=NULL, ...)

}

\arguments{

  Since \code{GeneRegionTrack} object are essentially just a specific
  type of \code{\linkS4class{AnnotationTrack}}, their constructors are
  quite similar. However, in the case of the \code{GeneRegionTrack}
  certain assumptions are made about the type of grouping on different
  levels (see the \code{Details} section for more information). The
  natural representation for gene model in the Bioconductor world are
  \code{\linkS4class{TranscriptDb}} objects, and we tried to make it as
  straight forward as possible to create \code{GeneRegionTracks}
  starting from those. Building the object from individual function
  arguments is of course still possible.
  
  \item{range}{
    
    An optional meta argument to handle the different input types. If
    the \code{range} argument is missing, all the relevant information
    to create the object has to be provided as individual function
    arguments (see below).

    The different input options for \code{range} are:

    \describe{
      
      \item{}{A \code{TranscriptDb} object: all the necessary gene model
	information including exon locations, transcript groupings and
	associated gene ids are contained in \code{TranscriptDb}
	objects, and the coercion between the two is almost completely
	automated. If desired, the data to be fetched from the
	\code{TranscriptDb} object can be restricted using the
	constructor's \code{chromosome}, \code{start} and \code{end}
	arguments. See below for details. A direct coercion method
	\code{as(obj, "GeneRegionTrack")} is also available. A nice
	added benefit of this input option is that the UTR and coding
	region information that is part of the original
	\code{TranscriptDb} object is retained in the
	\code{GeneRegionTrack}.}
      
      \item{}{A \code{GRanges} object: the genomic ranges for the
	\code{GeneRegion} track as well as the optional additional
	\code{elementMetadata} columns \code{feature},
	\code{transcript}, \code{gene}, \code{exon} and \code{symbol}
	(see description of the individual function parameters below for
	details). Calling the constructor on a \code{GRanges} object
	without further arguments, e.g.
	\code{GeneRegionTrack(range=obj)} is equivalent to calling the
	coerce method \code{as(obj, "GeneRegionTrack")}.}

      \item{}{A \code{GRangesList} object: this is very similar to the
	previous case, except that the grouping information that is part
	of the list structure is preserved in the
	\code{GeneRegionTrack}. I.e., all the elements within one list
	item receive the same group id. For consistancy, there is also a
	coercion method from \code{GRangesLists} \code{as(obj,
	"GeneRegionTrack")}. Please note that unless the necessary
	information about gene ids, symbols, etc. is present in the
	individual \code{GRanges} meta data slots, the object will not
	be particularly useful, because all the identifiers will be set
	to a common default value.}

      \item{}{An \code{\linkS4class{IRanges}} object: almost identical
	to the \code{GRanges} case, except that the chromosome and
	strand information as well as all additional data has to be
	provided in the separate \code{chromosome}, \code{strand},
	\code{feature}, \code{transcript},\ code{symbol}, \code{exon} or
	\code{gene} arguments, because it can not be directly encoded in
	an \code{IRange} object. Note that only the former two are
	mandatory (if not provided explicitely the more or less
	reasonable default values \code{chromosome=NA} and
	\code{strand=*} are used, but not providing information about
	the gene-to-transcript relationship or the human-readble symbols
	renders a lot of the class' functionality useles.}

      \item{}{A \code{data.frame} object: the \code{data.frame} needs to
	contain at least the two mandatory columns \code{start} and
	\code{end} with the range coordinates. It may also contain a
	\code{chromosome} and a \code{strand} column with the chromosome
	and strand information for each range. If missing it will be
	drawn from the separate \code{chromosome} or \code{strand}
	arguments. In addition, the \code{feature}, \code{exon},
	\code{transcript}, \code{gene} and \code{symbol} data can be
	provided as additional columns. The above comments about
	potential default values also apply here.}

    }
  }	
  
  \item{start, end}{An integer scalar with the genomic start or end
    coordinate for the gene model range. If those is missing, the
    default value will automatically be the smallest (or largest value,
    respectively) in \code{rstarts} for the currently active
    chromosome. When building a \code{GeneRegionTrack} from a
    \code{TranscriptDb} object, these arguments can be used to subset
    the desired annotation data by genomic coordinates. Please note that
    in that case the \code{chromosome} parameter must also be set.}

  \item{rstarts}{An integer vector of the start coordinates for the
    actual gene model items, i.e., for the individual exons. The
    relationship between exons is handled via the \code{gene} and
    \code{transcript} factors. Alternatively, this can be a vector of
    comma-separated lists of integer coordinates, one vector item for
    each transcript, and each comma-separated element being the start
    location of a single exon within that transcript. Those lists will
    be exploded upon object instantiation and all other annotation
    arguments will be recycled accordingly to regenerate the
    exon/transcript/gene relationship structure. This implies the
    approriate number of items in all annotation and coordinates
    arguments.}
  
  \item{rends}{An integer vector of the end coordinates for the actual
    gene model items. Both \code{rstarts} and \code{rends} have to be of
    equal length.}

   \item{rwidths}{An integer vector of widths for the actual gene model
    items. This can be used instead of either \code{rstarts} or
    \code{rends} to specify the range coordinates.}

   \item{feature}{Factor (or other vector that can be coerced into one),
    giving the feature types for the individual track exons.  When
    plotting the track to the device, if a display parameter with the
    same name as the value of \code{feature} is set, this will be used
    as the track item's fill color. See \code{\link{grouping}} for
    details.}

  \item{exon}{Character vector of exon identifiers. It's values will be
    used as the identifier tag when plotting to the device if the
    display parameter \code{showExonId=TRUE}.}

  \item{strand}{Character vector, the strand information for the
    individual track exons. It may be provided in the form \code{+} for
    the Watson strand, \code{-} for the Crick strand or \code{*} for
    either one of the two. Please note that all items within a single
    gene or transcript model need to be on the same strand, and
    erroneous entries will result in casting of an error.}

  \item{transcript}{Factor (or other vector that can be coerced into
    one), giving the transcript memberships for the individual track
    exons. All items with the same transcript identifier will be
    visually connected when plotting to the device.  See
    \code{\link{grouping}} for details. Will be used as labels when
    \code{showId=TRUE}, and \code{geneSymbol=FALSE}.}

  \item{gene}{Factor (or other vector that can be coerced into one),
    giving the gene memberships for the individual track exons.}

  \item{symbol}{ A factor with human-readable gene name aliases which
    will be used as labels when \code{showId=TRUE}, and
    \code{geneSymbol=TRUE}.}

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
    different active chromosome. When creating a \code{GeneRegionTrack}
    from a \code{TranscriptDb} object, the value of this parameter can
    be used to subset the data to fetch by chromosome.}

  \item{genome}{The genome on which the track's ranges are
    defined. Usually this is a valid UCSC genome identifier, however
    this is not being formally checked at this point. If not provided
    here the constructor will try to extract this information from the
    provided input, and eventually will fall back to the default value
    of \code{NA}.}
  
  \item{stacking}{The stacking type for overlapping items of the
    track. One in \code{c(hide, dense, squish, pack,full)}. Currently,
    only hide (don't show the track items, squish (make best use of the
    available space) and dense (no stacking at all) are implemented.}
  
  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

  \item{\dots}{Additional items which will all be interpreted as further
    display parameters. See \code{\link{settings}} and the "Display
    Parameters" section below for details.}
  
}

\details{

  A track containing all gene models in a particular region. The data
  are usually fetched dynamially from an online data store, but it is
  also possible to manully construct objects from local
  data. Connections to particular online data sources should be
  implemented as sub-classes, and \code{GeneRegionTrack} is just the
  commone denominator that is being used for plotting later on. There
  are several levels of data associated to a \code{GeneRegionTrack}:

 \describe{

   \item{exon level:}{identifiers are stored in the exon column of the
     \code{\linkS4class{GRanges}} object in the \code{range} slot. Data
     may be extracted using the \code{exon} method.}
   
    \item{transcript level:}{identifiers are stored in the transcript
      column of the \code{\linkS4class{GRanges}} object. Data may be
      extracted using the \code{transcript} method.}
      
    \item{gene level:}{identifiers are stored in the gene column of the
      \code{\linkS4class{GRanges}} object, more human-readable versions
      in the symbol column. Data may be extracted using the \code{gene}
      or the \code{symbol} methods.}
    
    \item{transcript-type level:}{information is stored in the feature
      column of the \code{\linkS4class{GRanges}} object.  If a display
      parameter of the same name is specified, the software will use its
      value for the coloring.}
  
  }

  \code{GeneRegionTrack} objects also know about coding regions and
  non-coding regions (e.g., UTRs) in a transctipt, and will indicate
  those by using different shapes (wide boxes for all coding regions,
  thinner boxes for non-coding regions). This is archived by setting the
  \code{feature} values of the object for non-coding elements to one of
  the options that are provided in the \code{thinBoxFeature} display
  parameters. All other elements are expected to be coding elements.
  
}

\value{

  The return value of the constructor function is a new object of class
  \code{GeneRegionTrack}.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{start}:}{Object of class \code{"numeric"}, the start
      coordinates of the annotation range. The coorrdinates for the
      individual gene model items are stored in the \code{range} slot.}
    
    \item{\code{end}:}{Object of class \code{"numeric"}, the end
      coordinates of the annotation range. The corrdinates for the
      individual gene model items are stored in the \code{range} slot. }

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
  
  Class \code{"\linkS4class{AnnotationTrack}"}, directly.
  
  Class \code{"\linkS4class{StackedTrack}"}, by class "AnnotationTrack",
  distance2.
  
  Class \code{"\linkS4class{RangeTrack}"}, by class "AnnotationTrack",
  distance3.
  
  Class \code{"\linkS4class{GdObject}"}, by class "AnnotationTrack",
  distance4.
  
}

\section{Methods}{

  In the following code chunks, \code{obj} is considered to be an object of class \code{GeneRegionTrack}.
	
  \bold{\emph{Exported in the name space:}}

  \describe{

     \item{group}{\code{signature(gdObject="GeneRegionTrack")}: extract
      the group membership for all track items.
		
      \emph{Usage:}
        
      \code{group(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{group(obj)}}

      }
    }
		
    \item{group<-}{\code{signature(gdObject="GeneRegionTrack",
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
	    
    \item{identifier}{\code{signature(gdObject="GeneRegionTrack")}:
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
	    
    \item{identifier<-}{\code{signature(gdObject="GeneRegionTrack",
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

    \item{exon}{\code{signature(GdObject="GeneRegionTrack")}: Extract
      the exon identifiers for all exons in the gene models.
      
      \emph{Usage:}
        
      \code{exon(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{exon(obj)}}

      }
    }
    
    \item{exon<-}{\code{signature(GdObject="GeneRegionTrack",
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
      
    \item{gene}{\code{signature(GdObject="GeneRegionTrack")}: Extract
      the gene identifiers for all gene models.
      
      \emph{Usage:}
        
      \code{gene(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{gene(obj)}}

      }
    }
    
    \item{gene<-}{\code{signature(GdObject="GeneRegionTrack",
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

    \item{symbol}{\code{signature(GdObject="GeneRegionTrack")}: Extract
      the human-readble gene symbol for all gene models.
	
      \emph{Usage:}
        
      \code{symbol(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{symbol(obj)}}

      }
    }
      
    \item{symbol<-}{\code{signature(GdObject="GeneRegionTrack",
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

    \item{transcript}{\code{signature(GdObject="GeneRegionTrack")}:
      Extract the transcript identifiers for all transcripts in the
      gene models.
	
      \emph{Usage:}
        
      \code{transcript(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{transcript(obj)}}

      }
    }
    
    \item{transcript<-}{\code{signature(GdObject="GeneRegionTrack",
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
	    
    \item{coerce}{\code{signature(from="GeneRegionTrack",
	to="UCSCData")}: coerce to a \code{UCSCData} object for export to
      the UCSC genome browser.
  
       \emph{Examples:}
      
      \describe{
	
	\item{}{\code{as(obj, "UCSCData")}}

      }
    }

    \item{collapseTrack}{\code{signature(GdObject="GeneRegionTrack")}:
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
    
    \item{initialize}{\code{signature(.Object="GeneRegionTrack")}:
      initialize the object }

    \item{show}{\code{signature(object="GeneRegionTrack")}: show a
      human-readable summary of the object }
    
  }

   \bold{\emph{Inherited methods:}}

   \describe{

     \item{drawGD}{\code{signature(GdObject="GeneRegionTrack")}: plot the
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

    \item{drawGrid}{\code{signature(GdObject="GeneRegionTrack")}: superpose a grid on top of a track.

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

    \item{setStacks}{\code{signature(GdObject="GeneRegionTrack")}:
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
	  
    \item{stacking}{\code{signature(GdObject="GeneRegionTrack")}: return
      the current stacking type.
    	
      \emph{Usage:}
      
      \code{stacking(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{stacking(obj)}}

      }
    }
    	
    \item{stacking<-}{\code{signature(GdObject="GeneRegionTrack",
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

    \item{stacks}{\code{signature(GdObject="GeneRegionTrack")}: return
      the stack indices for each track item.
      
      \emph{Usage:}
      
      \code{stacks(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::stacks(obj)}}
	
      }
    }
    
    \item{[}{\code{signature(x="GeneRegionTrack", i="ANY", j="ANY",
	drop="ANY")}: subset the items in the \code{GeneRegionTrack}
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
      
      
    \item{chromosome}{\code{signature(GdObject="GeneRegionTrack")}:
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
    
    \item{chromosome<-}{\code{signature(GdObject="GeneRegionTrack")}:
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
    
    \item{start, end, width}{\code{signature(x="GeneRegionTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="GeneRegionTrack")}:
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
    
    \item{position}{\code{signature(GdObject="GeneRegionTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{feature}{\code{signature(GdObject="GeneRegionTrack")}: return the
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
    
    \item{feature<-}{\code{signature(gdObject="GeneRegionTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{GeneRegionTrack}. See \code{\link{grouping}}
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
    
    \item{genome}{\code{signature(x="GeneRegionTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="GeneRegionTrack")}: set the track's
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

    \item{length}{\code{signature(x="GeneRegionTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="GeneRegionTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="GeneRegionTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{split}{\code{signature(x="GeneRegionTrack")}: split a
      \code{GeneRegionTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of objects of the same class as the input
      object, all inheriting from class \code{GeneRegionTrack}.
      
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
    
    \item{strand}{\code{signature(x="GeneRegionTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="GeneRegionTrack")}: replace the
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
    
    \item{values}{\code{signature(x="GeneRegionTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="GeneRegionTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{subset}{\code{signature(x="GeneRegionTrack")}: subset a
      \code{GeneRegionTrack} by coordinates and sort if necessary.
      
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
    
    \item{displayPars}{\code{signature(x="GeneRegionTrack",
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
      
    
    \item{displayPars}{\code{signature(x="GeneRegionTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="GeneRegionTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="GeneRegionTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="GeneRegionTrack",
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
    
    \item{setPar}{\code{signature(x="GeneRegionTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{GeneRegionTrack} class are
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
    
    \item{setPar}{\code{signature(x="GeneRegionTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{GeneRegionTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{names}{\code{signature(x="GeneRegionTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="GeneRegionTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="GeneRegionTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="GeneRegionTrack")}: return the tags from the
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
  \code{GeneRegionTrack} upon instantiation, unless one or more of them
  have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters
  for tracks. \describe{

    \item{}{\code{min.distance=0}: Numeric scalar. The minimum pixel
      distance before collapsing range items, only if
      \\code{collapse==TRUE}. See \\code{\\link{collapsing}} for
      details. Note that a value larger than 0 may lead to UTR regions
      being merged to CDS regions, which in most cases is not particularly
      useful.}

    \item{}{\code{col=NULL}: Character or integer scalar. The border
      color for all items. Defaults to using the same color as in
      \code{fill}, also taking into account different track
      \code{features}.}

    \item{}{\code{fill="orange"}: Character or integer scalar. The fill
      color for untyped items. This is also used to connect grouped
      items. See \code{\link{grouping}} for details.} 

    \item{}{\code{geneSymbols=TRUE}: Logical scalar. Use human-readable
      gene symbols or gene IDs for the transcript annotation.} 

    \item{}{\code{shape=c("smallArrow", "box")}: Character scalar. The
      shape in which to display the track items. Currently only
      \code{box}, \code{arrow}, \code{ellipse}, and \code{smallArrow}
      are implemented.} 

    \item{}{\code{showExonId=FALSE}: Logical scalar. Control whether
      to plot the individual exon identifiers.}

    \item{}{\code{collapseTranscripts=FALSE}: Logical scalar. Merge all
      transcripts of the same gene into one single gene
      model. Essentially, this will only keep the start location of the
      first exon and the end location of the last exon from all
      transcripts of a gene.}

      \item{}{\code{thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5",
	  "miRNA", "lincRNA")}: Character vector. A listing of feature
	  types that should be drawn with thin boxes. Typically those
	  are non-coding elements.}

  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{GeneRegionTrack} objects. 

    \describe{ 

    \item{}{\code{\linkS4class{AnnotationTrack}}: 

      \describe{ 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{col="transparent"}: Character or integer scalar.
          The border color for all track items.}
	
	\item{}{\code{col.line="darkgray"}: Character scalar. The color used for
	  connecting lines between grouped items. Defaults to a dark gray,
	  but if set to \code{NULL} the same color as for the first item in
	  the group is used.}

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


\author{Florian Hahne, Steve Lianoglou}

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

## The empty object
GeneRegionTrack()

## Load some sample data
data(cyp2b10)

## Construct the object
grTrack <- GeneRegionTrack(start=26682683, end=26711643,
rstart=cyp2b10$start, rends=cyp2b10$end, chromosome=7, genome="mm9",
transcript=cyp2b10$transcript, gene=cyp2b10$gene, symbol=cyp2b10$symbol,
name="Cyp2b10", strand=cyp2b10$strand)

## From a TranscriptDb object
if(require(GenomicFeatures)){
samplefile <- system.file("extdata", "UCSC_knownGene_sample.sqlite", package="GenomicFeatures")
txdb <- loadDb(samplefile)
GeneRegionTrack(txdb)
GeneRegionTrack(txdb, chromosome="chr6", start=300000,
end=350000)
}



\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(grTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
plotTracks(grTrack)

## Track names
names(grTrack)
names(grTrack) <- "foo"
plotTracks(grTrack)

## Subsetting and splitting
subTrack <- subset(grTrack, from=26700000, to=26705000)
length(subTrack)
subTrack <- grTrack[transcript(grTrack)=="ENSMUST00000144140"]
split(grTrack, transcript(grTrack))

## Accessors
start(grTrack)
end(grTrack)
width(grTrack)
position(grTrack)
width(subTrack) <- width(subTrack)+100

strand(grTrack)
strand(subTrack) <- "-"

chromosome(grTrack)
chromosome(subTrack) <- "chrX"

genome(grTrack)
genome(subTrack) <- "hg19"

range(grTrack)
ranges(grTrack)

## Annotation
identifier(grTrack)
identifier(grTrack, lowest=TRUE)
identifier(subTrack) <- "bar"

feature(grTrack)
feature(subTrack) <- "foo"

exon(grTrack)
exon(subTrack) <- letters[1:2]

gene(grTrack)
gene(subTrack) <- "bar"

symbol(grTrack)
symbol(subTrack) <- "foo"

transcript(grTrack)
transcript(subTrack) <- c("foo", "bar")
chromosome(subTrack) <- "chr7"
plotTracks(subTrack)

values(grTrack)

## Grouping
group(grTrack)
group(subTrack) <- "Group 1"
transcript(subTrack)
plotTracks(subTrack)

## Collapsing transcripts
plotTracks(grTrack, collapseTranscripts=TRUE, showId=TRUE,
extend.left=10000, shape="arrow")

## Stacking
stacking(grTrack)
stacking(grTrack) <- "dense"
plotTracks(grTrack)

## coercion
as(grTrack, "data.frame")
as(grTrack, "UCSCData")

## HTML image map
coords(grTrack)
tags(grTrack)
grTrack <- plotTracks(grTrack)$foo
coords(grTrack)
tags(grTrack)

}

\keyword{classes}






























