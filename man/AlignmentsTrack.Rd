\name{AlignmentsTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{AlignmentsTrack-class}
\alias{AlignmentsTrack}
\alias{chromosome<-,AlignmentsTrack-method}
\alias{drawAxis,AlignmentsTrack-method}
\alias{drawGD,AlignmentsTrack-method}

\alias{[,AlignmentsTrack,ANY,ANY-method}
\alias{coverage,AlignmentsTrack-method}
\alias{initialize,AlignmentsTrack-method}
\alias{initialize,ReferenceAlignmentsTrack-method}
\alias{subset,AlignmentsTrack-method}
\alias{subset,ReferenceAlignmentsTrack-method}
\alias{stacks,AlignmentsTrack-method}
\alias{values,AlignmentsTrack-method}
\alias{show,AlignmentsTrack-method}

\title{AlignmentsTrack class and methods}

\description{
  
  A class to represent short sequences that have been aligned to a
  reference genome as they are typically generated in next generation
  sequencing experiments.
  
}

\section{Objects from the Class}{

  Objects can be created using the constructor function \code{AlignmentsTrack}.

}


\usage{

AlignmentsTrack(range=NULL, start=NULL, end=NULL, width=NULL, strand, chromosome, genome,
                stacking="squish", id, cigar, mapq, flag, isize, groupid, status, md, seqs,
                name="AlignmentsTrack", isPaired=TRUE, importFunction, referenceSequence, ...)

}

\arguments{

  In the most common case \code{AlignmentsTrack} objects will be created
  directly from \code{BAM} files, and we strongly recommend to do
  this. A \code{BAM} file contains all the information that is needed to
  properly display the aligned reads, but more importantly, it allows to
  dynamically stream the data for the desired plotting range off the
  disk rather than having to load potentially gigantic amounts of data
  into memory upon object instantiation. That being said, there are
  other starting points to build \code{AlignmentsTracks}.

   \item{range}{
    
    An optional meta argument to handle the different input types. If
    the \code{range} argument is missing, all the relevant information
    to create the object has to be provided as individual function
    arguments (see below).

    The different input options for \code{range} are:

    \describe{

      \item{}{A \code{character} string: the path to a \code{BAM} file
	containing the read alignments. To be precise, this will result
	in the instantiation of a \code{ReferenceAlignmentsTrack}
	object, but for the user this implementation detail should be of
	no concern.}

      \item{}{A \code{GRanges} object: the genomic ranges of the
	individual reads as well as the optional additional
	metadata columns \code{id}, \code{cigar},
	\code{mapq}, \code{flag}, \code{isize}, \code{groupid},
	\code{status}, \code{md} and \code{seqs} (see description of the
	individual function parameters below for details). Calling the
	constructor on a \code{GRanges} object without further
	arguments, e.g. \code{AlignmentsTrack(range=obj)} is equivalent
	to calling the coerce method \code{as(obj, "AlignmentsTrack")}.}

      \item{}{An \code{\linkS4class{IRanges}} object: almost identical
	to the \code{GRanges} case, except that the chromosome and
	strand information as well as all additional metadata has to be
	provided in the separate \code{chromosome}, \code{strand},
	\code{feature}, \code{group} or \code{id} arguments, because it
	can not be directly encoded in an \code{IRanges} object. Note
	that none of those inputs are mandatory, and if not provided
	explicitely the more or less reasonable default values
	\code{chromosome=NA} and \code{strand="*"} are used. }

      \item{}{A \code{data.frame} object: the \code{data.frame} needs to
	contain at least the two mandatory columns \code{start} and
	\code{end} with the range coordinates. It may also contain a
	\code{chromosome} and a \code{strand} column with the chromosome
	and strand information for each range. If missing it will be
	drawn from the separate \code{chromosome} or \code{strand}
	arguments. In addition, the \code{id}, \code{cigar},
	\code{mapq}, \code{flag}, \code{isize}, \code{groupid},
	\code{status}, \code{md} and \code{seqs} data can be provided as
	additional columns. The above comments about potential default
	values also apply here.}

    }
  }

  \item{start, end, width}{Integer vectors, giving the start and the
    end coordinates for the individual track items, or their
    width. Two of the three need to be specified, and have to be of
    equal length or of length one, in which case this single value
    will be recycled. Otherwise, the usual R recycling rules for
    vectors do not apply here.}

  \item{id}{Character vector of read identifiers. Those identifiers have
    to be unique, i.e., each range representing a read needs to have a
    unique \code{id}.}

  \item{cigar}{A character vector of valid CIGAR strings describing
    details of the alignment. Typically those include alignemnts gaps or
    insertions and deletions, but also hard and soft clipped read
    regions. If missing, a fully mapped read without gaps or indels is
    assumed. Needs to be of equal length as the provided genomic
    coordinates, or of length 1.}

  \item{mapq}{A numeric vector of read mapping qualities. Needs to be of
    equal length as the provided genomic coordinates, or of length 1.}

  \item{flag}{A numeric vector of flag values. Needs to be of equal
    length as the provided genomic coordinates, or of length
    1. Currently not used.}

  \item{isize}{A numeric vector of empirical insert sizes. This only
    applies if the reads are paired. Needs to be of equal length as the
    provided genomic coordinates, or of length 1. Currently not used.}
  
  \item{groupid}{A factor (or vector than can be coerced into one)
    defining the read pairs. Reads with the same \code{groupid} are
    considered to be mates. Please note that each read group may only
    have one or two members. Needs to be of equal length as the provided
    genomic coordinates, or of length 1.}

  \item{status}{A factor describing the mapping status of a read. Has to
    be one in \code{mated}, \code{unmated} or \code{ambiguous}. Needs to
    be of equal length as the provided genomic coordinates, or of length
    1.}
  
  \item{md}{A character vector describing the mapping details. This is
    effectively and alternative to the CIGAR encoding and it removes the
    dependency on a reference sequence to figure out read
    mismatches. Needs to be of equal length as the provided genomic
    coordinates, or of length 1. Currently not used.}
  
  \item{seqs}{\code{DNAStringSet} of read sequences.}


  \item{strand}{Character vector, the strand information for the
    reads. It may be provided in the form \code{+} for the Watson
    strand, \code{-} for the Crick strand or \code{*} for either one of
    the two. Needs to be of equal length as the provided genomic
    coordinates, or of length 1. Please note that paired reads need to
    be on opposite strands, and erroneous entries will result in casting
    of an error.}

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
    track. One in \code{c(hide, dense, squish, pack, full)}. Currently,
    only squish (make best use of the available space), dense (no
    stacking, collapse overlapping ranges), and hide (do not show any
    track items at all) are implemented.}

  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}

  \item{isPaired}{A logical scalar to determine whether the reads are
    paired or not. While this may be used to render paired-end data as
    single-end, the oppsite will typically not have any effect because
    the appropriate \code{groupid} settings will not be present.  Thus
    setting \code{isPaired} to \code{TRUE} can usually be used to
    autodetect the pairing state of the input data.}

  \item{importFunction}{A user-defined function to be used to import the
    data from a file. This only applies when the \code{range} argument
    is a character string with the path to the input data file. The
    function needs to accept an argument \code{x} containing the file
    path and a second argument \code{selection} with the desired
    plotting ranges. It has to return a proper \code{GRanges} object
    with all the necessary metadata columns set. A single
    default import function is already implemented in the package for
    \code{BAM} files.}

  \item{referenceSequence}{An optional
    \code{\linkS4class{SequenceTrack}} object containing the reference
    sequence against which the reads have been aligned. This is only
    needed when mismatch information has to be added to the plot (i.e.,
    the \code{showMismatchs} display parameter is \code{TRUE}) because
    this is normally not encoded in the \code{BAM} file. If not provided
    through this argument, the \code{\link{plotTracks}} function is
    smart enough to detect the presence of a
    \code{\linkS4class{SequenceTrack}} object in the track list and will
    use that as a reference sequence.}
  
  \item{\dots}{Additional items which will all be interpreted as further
    display parameters. See \code{\link{settings}} and the "Display
    Parameters" section below for details.}
}


\value{

  The return value of the constructor function is a new object of class
  \code{AlignmentsTrack} or \code{ReferenceAlignmentsTrack}.
  
}

\section{details}{
  \code{AlignmentTracks} usually have two section: the coverage section
  on top showing a histogram of the read coverage, and the pile-up
  section below with the individual reads. Both can be toggled on or off
  using the \code{type} display parameter. If reference sequence has
  been provided either during object instantiation or with the track
  list to the call to \code{plotTracks}, sequence mismatch information
  will be shown in both sections: as a stacked histogram in the coverage
  plot and as colored boxes or characters (depending on available space)
  in for the pile-ups.

}


\section{Slots}{
  
  \describe{
    
    \item{\code{stackRanges}:}{Object of class \code{"GRanges"}, the
      ranges of the precomputed mate or gaps stacks that should remain
      conmnected.}

    \item{\code{sequences}:}{Object of class \code{"DNAStringSet"}, the
      processed read sequences.}

    \item{\code{referenceSequence}:}{Object of class \code{"SequenceTrack"}, the
      reference sequence to which the reads have been aligned to.}
    
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
  of class \code{AlignmentsTrack}.
  
  \bold{\emph{Exported in the name space:}}
  
  \describe{
    
    \item{[}{\code{signature(x="AlignmentsTrack")}: subset the items in
      the \code{AlignmentsTrack}. This is essentially similar to
      subsetting of the \code{\linkS4class{GRanges}} object in the
      \code{range} slot. For most applications, the \code{subset} method
      may be more appropriate.

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{i}: subsetting indices}
	
      }
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{obj[1:5]}}
	
      }
    }

    \item{subset}{\code{signature(x="AlignmentsTrack")}: subset a
      \code{AlignmentsTrack} by coordinates and sort if necessary.
      
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
      
    \item{split}{\code{signature(x="AlignmentsTrack")}: split an
      \code{AlignmentsTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of \code{AlignmentsTrack} objects.
      
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
      
  }

  \bold{\emph{Internal methods:}}
  
  \describe{
    
    
    \item{drawAxis}{\code{signature(GdObject="AlignmentsTrack")}: add a y-axis
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
   
    \item{drawGD}{\code{signature(gdObject="AlignmentsTrack")}: plot the
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

     
    \item{initialize}{\code{signature(.Object="AligendReadTrack")}:
      initialize the object.}

    \item{show}{\code{signature(object="AlignmentsTrack")}: show a
      human-readable summary of the object.}
        
  }

  \bold{\emph{Inherited methods:}}

  \describe{
	  
    \item{stacking}{\code{signature(GdObject="AlignmentsTrack")}: return
      the current stacking type.
    	
      \emph{Usage:}
      
      \code{stacking(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{stacking(obj)}}

      }
    }
    	
    \item{stacking<-}{\code{signature(GdObject="AlignmentsTrack",
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
    
    \item{setStacks}{\code{signature(GdObject="AlignmentsTrack")}:
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
    
    \item{stacks}{\code{signature(GdObject="AlignmentsTrack")}: return
      the stack indices for each track item.

      \emph{Usage:}
      
      \code{stacks(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
      	\item{}{\code{Gviz:::stacks(obj)}}

      }
    }
    
    \item{chromosome}{\code{signature(GdObject="AlignmentsTrack")}:
      return the chromosome for which the track is defined.
      
      \emph{Usage:}
      
      \code{chromosome(GdObject)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{chromosome(obj)}}
	
      }
    }
    
    \item{chromosome<-}{\code{signature(GdObject="AlignmentsTrack")}:
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
    
    \item{start, end, width}{\code{signature(x="AlignmentsTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="AlignmentsTrack")}:
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
    
    \item{position}{\code{signature(GdObject="AlignmentsTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}
      
      \code{position(GdObject)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{position(obj)}}

      }
    }
    
    \item{feature}{\code{signature(GdObject="AlignmentsTrack")}: return the
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
    
    \item{feature<-}{\code{signature(gdObject="AlignmentsTrack",
      value="character")}: set the grouping information for track
      items. This has to be a factor vector (or another type of vector
      that can be coerced into one) of the same length as the number of
      items in the \code{AlignmentsTrack}. See \code{\link{grouping}} or
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
    
    \item{genome}{\code{signature(x="AlignmentsTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="AlignmentsTrack")}: set the track's
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

    \item{length}{\code{signature(x="AlignmentsTrack")}: return the number
      of items in the track.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
    
    \item{range}{\code{signature(x="AlignmentsTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.
      
      \emph{Usage:}
      
      \code{range(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{range(obj)}}

      }
    }
    
    \item{ranges}{\code{signature(x="AlignmentsTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.
      
      \emph{Usage:}
      
      \code{ranges(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{ranges(obj)}}

      }
    }
    
    \item{strand}{\code{signature(x="AlignmentsTrack")}: return a vector of
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
    
    \item{strand<-}{\code{signature(x="AlignmentsTrack")}: replace the
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
    
    \item{values}{\code{signature(x="AlignmentsTrack")}: return all
      additional annotation information except for the genomic coordinates
      for the track items as a data.frame.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{values(obj)}}

      }
    }
    
    \item{coerce}{\code{signature(from="AlignmentsTrack",
	to="data.frame")}: coerce the \code{\linkS4class{GRanges}}
	object in the \code{range} slot into a regular data.frame.
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{as(obj, "data.frame")}}

      }
    }

    \item{displayPars}{\code{signature(x="AlignmentsTrack",
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
      
    
    \item{displayPars}{\code{signature(x="AlignmentsTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="AlignmentsTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="AlignmentsTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="AlignmentsTrack",
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
    
    \item{setPar}{\code{signature(x="AlignmentsTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{AlignmentsTrack} class are
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
    
    \item{setPar}{\code{signature(x="AlignmentsTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{AlignmentsTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{group}{\code{signature(GdObject="AlignmentsTrack")}: return
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
    
    \item{names}{\code{signature(x="AlignmentsTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="AlignmentsTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="AlignmentsTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="AlignmentsTrack")}: return the tags from the
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
  \code{AlignmentsTrack} upon instantiation, unless one or more of
  them have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters
  for tracks.

  \describe{

    \item{}{\code{alpha.mismatch=NULL}: Numeric scalar between 0 and 1. The
      transparency of the mismatch base information.}

    \item{}{\code{alpha.reads=0.5}: Numeric scalar between 0 and 1. The
      transparency of the individual read icons. Can be used to indicate
      overlapping regions in read pairs. Only on supported devices.}

    \item{}{\code{cex.mismatch=0.7}: Numeric Scalar. The character
      expansion factor for the mismatch base letters.}

    \item{}{\code{col.coverage=NULL}: Integer or character
      scalar. The line color for the coverage profile.}
    
    \item{}{\code{col.gap=.DEFAULT_SHADED_COL}: Integer or character
      scalar. The color of the line that is bridging the gap regions in
      gapped alignments.}
    
    \item{}{\code{col.mates=.DEFAULT_BRIGHT_SHADED_COL}: Integer or
      character scalar. The color of the line that is connecting two
      paired reads.}

    \item{}{\code{col.mismatch=.DEFAULT_SHADED_COL}: Integer or
      character scalar. The box color around mismatch bases.}

    \item{}{\code{col.reads=NULL}: Integer or
      character scalar. The box color around reads.}
    
    \item{}{\code{col=.DEFAULT_SHADED_COL}: Integer or character
      scalar. The default color of all line elements.}
    
    \item{}{\code{collapse=FALSE}: Logical scalar. Do not perform any
      collapsing of overlapping elements. Currently not supported.}
    
    \item{}{\code{coverageHeight=0.1}: Numeric scalar. The height of the
      coverage region of the track. Can either be a value between 0 and
      1 in which case it is taken as a relative height, or a positive
      value greater 1 in which case it is interpreted as pixels.}

    \item{}{\code{fill.coverage=NULL}: Integer or character
      scalar. The fill color for the coverage profile.}

    \item{}{\code{fill.reads=NULL}: Integer or character
      scalar. The fill color for the read icons.}
      
    \item{}{\code{fill="#BABABA"}: Integer or character
      scalar. The default fill color of all plot elements.}
    
    \item{}{\code{fontface.mismatch=2}: Integer scalar. The font face
      for mismatch bases.}

    \item{}{\code{lty.coverage=NULL}: Integer or character scalar. The line
      type of the coverage profile.}
    
    \item{}{\code{lty.gap=NULL}: Integer or character scalar. The type
      of the line that is bridging the gap regions in gapped
      alignments.}
    
    \item{}{\code{lty.mates=NULL}: Integer or character scalar. The
      type of the line that is connecting two paired reads.}

    \item{}{\code{lty.mismatch=NULL}: Integer or character scalar. The box
      line type around mismatch bases.}

    \item{}{\code{lty.reads=NULL}: Integer or character scalar. The box
      line type around reads.}
    
    \item{}{\code{lty=1}: Integer or character scalar. The default type
      of all line elements.}
    
    \item{}{\code{lwd.coverage=NULL}: Integer or character scalar. The line
      width of the coverage profile.}
    
    \item{}{\code{lwd.gap=NULL}: Integer scalar. The width of the line
      that is bridging the gap regions in gapped alignments.}
    
    \item{}{\code{lwd.mates=NULL}: Integer scalar. The width of the
      line that is connecting two paired reads.}

    \item{}{\code{lwd.mismatch=NULL}: Integer scalar.The box
      line width around mismatch bases.}

    \item{}{\code{lwd.reads=NULL}: Integer scalar.The box
      line width around reads.}
    
    \item{}{\code{lwd=1}: Integer scalar. The default width of all line
      elements.}

    \item{}{\code{noLetters=FALSE}: Logical scalar. Always plot colored
      boxes for mismatch bases regardles of the available space.}
    
    \item{}{\code{max.height=10}: Integer scalar. The maximum height
      of an individual read in pixels. Can be used in combination with
      \code{min.height} to control the read and stacking appearance. }
    
    \item{}{\code{min.height=5}: Integer scalar. The minimum height of
      an individual read in pixels. Can be used in combination with
      \code{max.height} to control the read and stacking appearance.}
    
    \item{}{\code{minCoverageHeight=50}: Integer scalar. The minimum
      height of the coverage section. Uselful in combination with a
      relative setting of \code{coverageHeight}.}
    
    \item{}{\code{showMismatches=TRUE}: Logical scalar. Add mismatch
      information, either as individual base letters or using color
      coded bars. This implies that the reference sequence has been
      provided, either to the class constructor or as part of the
      track list.}
    
    \item{}{\code{size=NULL}: Numeric scalar. The size of the
      track. Defaults to automatic sizing.}
    
    \item{}{\code{type=c("coverage", "pileup")}: Character vactor. The
      type of information to plot. For \code{coverage} a coverage plot,
      potentially augmented by base mismatch information, and for
      \code{"pileup"} the pileups of the individual reads. The two can
      be combined.}
    
  } 

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{AlignmentsTrack} objects. 

  \describe{

    \item{}{\code{\linkS4class{StackedTrack}}:
      
      \describe{ 

	\item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse the
	  y-ordering of stacked items. I.e., features that are plotted on
	  the bottom-most stacks will be moved to the top-most stack and
	  vice versa.}

	\item{}{\code{stackHeight=0.75}: Numeric between 0 and 1. Controls
	  the vertical size and spacing between stacked elements. The number
	  defines the proportion of the total available space for the stack
	  that is used to draw the glyphs. E.g., a value of 0.5 means that
	  half of the available vertical drawing space (for each stacking
	  line) is used for the glyphs, and thus one quarter of the available
	  space each is used for spacing above and below the glyph. Defaults
	  to 0.75.}

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

        \item{}{\code{col.border.title="transparent"}: Integer or character
          scalar. The border color for the title panels.}  

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.}

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
## Creating objects
afrom <- 2960000
ato <- 3160000
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata",
"gapped.bam"), isPaired=TRUE)
plotTracks(alTrack, from=afrom, to=ato, chromosome="chr12")

## Omit the coverage or the pile-ups part
plotTracks(alTrack, from=afrom, to=ato, chromosome="chr12",
type="coverage")
plotTracks(alTrack, from=afrom, to=ato, chromosome="chr12",
type="pileup")

## Including sequence information with the constructor
if(require(BSgenome.Hsapiens.UCSC.hg19)){
strack <- SequenceTrack(Hsapiens, chromosome="chr21")
afrom <- 44945200
ato <- 44947200
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata",
"snps.bam"), isPaired=TRUE, referenceSequence=strack)
plotTracks(alTrack, chromosome="chr21", from=afrom, to=ato)

## Including sequence information in the track list
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata",
"snps.bam"), isPaired=TRUE)
plotTracks(c(alTrack, strack), chromosome="chr21", from=44946590,
to=44946660)
}

}
\keyword{classes}































