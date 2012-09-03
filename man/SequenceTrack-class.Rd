\name{SequenceTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{SequenceTrack-class}
\alias{SequenceTrack}
\alias{chromosome,SequenceTrack-method}
\alias{chromosome<-,SequenceTrack-method}
\alias{end,SequenceTrack-method}
\alias{genome,SequenceTrack-method}
\alias{initialize,SequenceTrack-method}
\alias{initialize,SequenceBSgenomeTrack-method}
\alias{initialize,SequenceDNAStringSetTrack-method}
\alias{seqlevels,SequenceBSgenomeTrack-method}
\alias{seqnames,SequenceBSgenomeTrack-method}
\alias{seqlevels,SequenceDNAStringSetTrack-method}
\alias{seqnames,SequenceDNAStringSetTrack-method}
\alias{length,SequenceTrack-method}
\alias{start,SequenceTrack-method}
\alias{subseq,SequenceTrack-method}
\alias{coerce,DNAString,Rle-method}
\alias{width,SequenceTrack-method}


\title{SequenceTrack class and methods}

\description{

  A track class to represent genomic sequences. The two child classes
  \code{SequenceDNAStringSetTrack} and \code{SequenceBSgenomeTrack} do
  most of the work, however in practise they are of no particular
  relevance to the user. 
  
}

\section{Objects from the class}{
  
  Objects can be created using the constructor function
  \code{SequenceTrack}.
  
}

\usage{

SequenceTrack(sequence = NULL, chromosome = NULL, genome, name = "Sequence", ...)

}

\arguments{
  
  \item{sequence}{Either an object of class \code{\linkS4class{DNAStringSet}} or an
    object of class \code{\linkS4class{BSgenome}}. In the former case,
    the individual \code{\linkS4class{DNAString}}s are considered to be
    the different chromosome sequences. In the latter case, the
    \code{Gviz} package tries to follow the \code{BSgenome} philosophy
    in that the respective chromosome sequences are only realized once
    they are first accessed.}

  \item{chromosome}{The currently active chromosome of the track. A
    valid UCSC chromosome identifier if
    \code{options(ucscChromosomeNames=TRUE)}. Please note that in this
    case only syntactic checking takes place, i.e., the argument value
    needs to be an integer, numeric character or a character of the form
    \code{chrx}, where \code{x} may be any possible string. The user has
    to make sure that sequences for the respective chromosomes are
    indeed part of the object. If not provided here, the constructor
    will set it to the first available sequence. Please note that by
    definition all objects in the \code{Gviz} package can only have a
    single active chromosome at a time (although internally the
    information for more than one chromosome may be present), and the
    user has to call the \code{chromosome<-} replacement method in order
    to change to a different active chromosome.}
  
  \item{genome}{The genome on which the track's ranges are
    defined. Usually this is a valid UCSC genome identifier, however
    this is not being formally checked at this point. For a
    \code{SequenceBSgenomeTrack} object, the genome information is
    extracted from the input \code{BSgenome} package. For a
    \code{DNAStringSet} it has too be provided or the constructor will
    fall back to the default value of \code{NA}.}
  
  \item{name}{Character scalar of the track's name used in the title
    panel when plotting.}
  
  \item{\dots}{Additional items which will all be interpreted as further
    display parameters. See \code{\link{settings}} and the "Display
    Parameters" section below for details.}

} 
    
\value{

  The return value of the constructor function is a new object of class
  \code{SequenceDNAStringSetTrack} or of class
  \code{SequenceBSgenomeTrack} depending on the constructor
  arguments. Typically the user will not have to be troubled with this
  distinction and can rely on the constructor to make the right choice.
  
}

\section{details}{

  Depending on the available space the class will use different options
  to plot a sequence. If single letters can be accomodated without
  overplotting those will be show. Otherwise, colored boxes will be used
  to indicate letters, and if there is not enough horizontal room to
  show those, a simple line will indicate presence of a sequence. The
  \code{min.width} and \code{fontsize} display parameters directly
  control this behaviour. Each of the five possible nucleotides (G, A, T,
  C, and N) will be endoded in a separate color. As default we use the
  colors suggested in the \code{biovizBase} package, but a user is
  free to set their own color scheme by providing a named character
  vector with color as display parameter \code{fontcolor}, with names
  equal to the five possible bases.

}

\section{Slots}{

  \describe{
    	
    \item{\code{chromosome}:}{Object of class \code{"character"}, the
      chromosome on which the track is defined. There can only be a
      single chromosome for one track. Throughout the
      package, chromosome name have to be entered either as a single
      integer scalar or as a character scalar of the form
      \code{chrXYZ}, where \emph{XYZ} may be an arbitrary
      character string.}
    
    \item{\code{genome}:}{Object of class \code{"character"}, the genome
      for which the track is defined. This should be a valid UCSC genome
      identifier, however this may not always be formally checked upon
      object instantiation.}
    
    \item{\code{dp}:}{Object of class
      \code{\linkS4class{DisplayPars}}, inherited from class
      \code{\linkS4class{GdObject}}.}
    
    \item{\code{name}:}{Object of class \code{"character"}, inherited
      from class \code{\linkS4class{GdObject}} }

    \item{\code{imageMap}:}{Object of class
      \code{\linkS4class{ImageMap}}, inherited from class
      \code{\linkS4class{GdObject}}}
  
  }
  
}

\section{Extends}{
  
  Class \code{"\linkS4class{GdObject}"}, directly.
  
}

\section{Methods}{
    
  In the following code chunks, \code{obj} is considered to be an object
  inheriting from class \code{SequenceTrack}.
  
  \bold{\emph{Exported in the name space:}}

  \describe{
      
      \item{chromosome}{\code{signature(GdObject="SequenceTrack")}:
	return the chromosome for which the track is defined.
      
	\emph{Usage:}
      
	\code{chromosome(GdObject)}

        \emph{Examples:}
      
	\describe{

	  \item{}{\code{chromosome(obj)}}

	}
      }
    
      \item{chromosome<-}{\code{signature(GdObject="SequenceTrack")}:
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
    
    \item{genome}{\code{signature(x="SequenceTrack")}: return the track's genome. 
      
      \emph{Usage:}
      
      \code{genome(x)}
      
      \emph{Examples:}
      
      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="SequenceTrack")}: set the track's
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

    \item{length}{\code{signature(x="SequenceTrack")}: return the number
      of nucleotides in the track's sequence.
      
      \emph{Usage:}
      
      \code{length(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{length(obj)}}

      }
    }
   
    \item{seqnames}{\code{signature(x="SequenceTrack")}: return the
    names (i.e., the chromosome) of the sequences contained in the object.
      
      \emph{Usage:}
      
      \code{values(x)}

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{seqnames(obj)}}

      }
    }

    \item{subseq}{\code{signature(x="SequenceTrack")}: Extract a
      sub-sequence from the track.
      
      \emph{Usage:}
	
      \code{subseq(x, start=NA, end=NA, width=NA)}

      \emph{Additional Arguments:}
      
      \describe{
	
	\item{}{\code{start}: the start coordinate for the sub-sequence.}

	\item{}{\code{end}: the end coordinate for the sub-sequence.}

	\item{}{\code{width}: the width of the sub-sequence.}
	  
      }
      
       \emph{Examples:}
      
      \describe{
  
	\item{}{\code{subseq(obj, 1, 10)}}

      }
    }
  }
	
  \bold{\emph{Internal methods:}}
  
  \describe{
    
    \item{initialize}{\code{signature(.Object="SequenceTrack")}:
      initialize the object. }
  }

   \bold{\emph{Inherited methods:}}
  
   \describe{

     \item{displayPars}{\code{signature(x="SequenceTrack",
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
      
    
    \item{displayPars}{\code{signature(x="SequenceTrack", name="missing")}:
      list the value of all available display parameters. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{displayPars(obj)}}

      }
    }
    
    \item{getPar}{\code{signature(x="SequenceTrack", name="character")}:
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
    
    \item{getPar}{\code{signature(x="SequenceTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.
  
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{getPar(obj)}}

      }
    }
    
    \item{displayPars<-}{\code{signature(x="SequenceTrack", value="list")}:
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
    
    \item{setPar}{\code{signature(x="SequenceTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{SequenceTrack} class are
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
    
    \item{setPar}{\code{signature(x="SequenceTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{SequenceTrack} class are pass-by-reference, so no re-assignmnet
      to the symbol \code{obj} is necessary. See \code{\link{settings}}
      for details on display parameters and customization.
      
       
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }
    
    \item{names}{\code{signature(x="SequenceTrack")}: return the value of
      the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj)}}

      }
    }
    
    \item{names<-}{\code{signature(x="SequenceTrack", value="character")}:
      set the value of the \code{name} slot.
      
      \emph{Usage:}
      
      \code{names<-(x, value)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="SequenceTrack")}: return the
      coordinates from the internal image map.
      
      \emph{Usage:}
      
      \code{coords(ImageMap)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="SequenceTrack")}: return the tags from the
      internal image map.
      
      \emph{Usage:}
      
      \code{tags(x)}
      
      \emph{Examples:}
      
      \describe{
	
	\item{}{\code{tags(obj)}}

      }
    }
    
    \item{drawAxis}{\code{signature(GdObject="SequenceTrack")}: add a
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
    
    \item{drawGrid}{\code{signature(GdObject="SequenceTrack")}: superpose a
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
  \code{SequenceTrack} upon instantiation
  
  \describe{ 
   
     \item{}{\code{size=null}: Numeric scalar. The size of
       the track item. Defaults to auto-detect the size based on the
       other parameter settings.}

     \item{}{\code{fontcolor=getBioColor("DNA_BASES_N")}: Character
       vector. The colors
       used for the 5 possible nucleotides (G, A, T, C, N). Defaults to
       use colors as defined in the \code{biovizBase} package.}

     \item{}{\code{fontsize=10}: Numeric scalar. Controls the size of
       the sequence and thus also the level of plotable details. }

     \item{}{\code{fontface=2}: Numeric scalar. The face of the font. }

     \item{}{\code{lwd=2}: Numeric scalar. The width of the line when no
       indiviual letters can be plotted due to size limitations.}

     \item{}{\code{col="darkgray"}: Character scalar. The color of the
       line when no indiviual letters can be plotted due to size
       limitations.}

     \item{}{\code{min.width=2}: Numeric scalar. The minimum width of
       the colored boxes that are drawn when no indiviual letters can be
       plotted due to size limitations.}

     \item{}{\code{showTitle=FALSE}: Logical scalar. Do not show a title
       panel by default.}
       
     \item{}{\code{background.title="transparent"}: Character
       scalar. Make the title panel transparent by default.}

     \item{}{\code{noLetters=FALSE}: Logical scalar. Always plot colored
       boxes (or a line) regardles of the available space.}
     
     \item{}{\code{add53=FALSE}: Logical scalar. Add a direction
       indicator.}
     
     \item{}{\code{add53=FALSE}: Logical scalar. Plot the sequence complement.}
     
   }

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{SequenceTrack} objects. 

    \describe{ 

    \item{}{\code{\linkS4class{GdObject}}: 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

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

  \code{\linkS4class{DataTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\linkS4class{GeneRegionTrack}}

  \code{\linkS4class{GRanges}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{IRanges}}

  \code{\linkS4class{BSgenome}}

  \code{\linkS4class{DNAStringSet}}

  \code{\link{plotTracks}}

  \code{\link{settings}}

}


\examples{
## An empty object
SequenceTrack()

## Construct from DNAStringSet
library(Biostrings)
letters <- c("A", "C", "T", "G", "N")
set.seed(999)
seqs <- DNAStringSet(c(chr1=paste(sample(letters, 100000, TRUE),
collapse=""), chr2=paste(sample(letters, 200000, TRUE), collapse="")))
sTrack <- SequenceTrack(seqs, genome="hg19")
sTrack

## Construct from BSGenome object
if(require(BSgenome.Hsapiens.UCSC.hg19)){
sTrack <- SequenceTrack(Hsapiens)
sTrack
}


## Set active chromosome
chromosome(sTrack)
chromosome(sTrack) <- "chr2"
head(seqnames(sTrack))



\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(sTrack) <- list(fontfamily=font, fontfamily.title=font)
}
}

## Plotting
## Sequences
plotTracks(sTrack, from=199970, to=200000)
## Boxes
plotTracks(sTrack, from=199800, to=200000)
## Line
plotTracks(sTrack, from=1, to=200000)
## Force boxes
plotTracks(sTrack, from=199970, to=200000, noLetters=TRUE)
## Direction indicator
plotTracks(sTrack, from=199970, to=200000, add53=TRUE)
## Sequence complement
plotTracks(sTrack, from=199970, to=200000, add53=TRUE, complement=TRUE)
## Colors
plotTracks(sTrack, from=199970, to=200000, add53=TRUE, fontcolor=c(A=1,
C=1, G=1, T=1, N=1))

## Track names
names(sTrack)
names(sTrack) <- "foo"

## Accessors
genome(sTrack)
genome(sTrack) <- "mm9"
length(sTrack)

## Sequence extraction
subseq(sTrack, start=100000, width=20)
## beyond the stored sequence range
subseq(sTrack, start=length(sTrack), width=20)

}




\keyword{classes}































