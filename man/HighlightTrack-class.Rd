\name{HighlightTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{HighlightTrack-class}
\alias{HighlightTrack}
\alias{initialize,HighlightTrack-method}
\alias{setStacks,HighlightTrack-method}
\alias{subset,HighlightTrack-method}
\alias{consolidateTrack,HighlightTrack-method}
\alias{show,HighlightTrack-method}
\alias{length,HighlightTrack-method}
\alias{chromosome<-,HighlightTrack-method}
\alias{displayPars<-,HighlightTrack,list-method}

\title{HighlightTrack class and methods}

\description{

  A container for other track objects from the Gviz package that allows
  for the addition of a common highlighting area across tracks.

}

\section{Objects from the Class}{

 Objects can be created using the constructor function
  \code{HighlightTrack}.

}

\usage{

HighlightTrack(trackList=list(), range=NULL, start=NULL, end=NULL, width=NULL, chromosome, genome,
                           name="HighlightTrack",  ...)


}

\arguments{

  \item{trackList}{A list of Gviz track objects that all have to inherit
    from class \code{GdObject}.}

  \item{range}{

    An optional meta argument to handle the different input types. If
    the \code{range} argument is missing, all the relevant information
    to create the object has to be provided as individual function
    arguments (see below).

    The different input options for \code{range} are:

    \describe{

      \item{}{A \code{GRanges} object: the genomic ranges for the
	highlighting regions.}

      \item{}{An \code{\linkS4class{IRanges}} object: almost identical
	to the \code{GRanges} case, except that the chromosome
	information has to be provided in the separate \code{chromosome}
	argument, because it can not be directly encoded in an
	\code{IRanges} object.}

      \item{}{A \code{data.frame} object: the \code{data.frame} needs to
	contain at least the two mandatory columns \code{start} and
	\code{end} with the range coordinates. It may also contain a
	\code{chromosome} column with the chromosome information for
	each range. If missing, this information will be drawn from the
	constructor's \code{chromosome} argument.}

    }
  }

  \item{start, end}{An integer scalar with the genomic start or end
    coordinates for the highlighting range. Can also be supplied as
    part of the \code{range} argument.}

  \item{width}{An integer vector of widths for highlighting
    ranges. This can be used instead of either \code{start} or
    \code{end} to specify the range coordinates.}

   \item{chromosome}{The chromosome on which the track's genomic ranges
    are defined. A valid UCSC chromosome identifier if
    \code{options(ucscChromosomeNames=TRUE)}. Please note that in this
    case only syntactic checking takes place, i.e., the argument value
    needs to be an integer, numeric character or a character of the form
    \code{chrx}, where \code{x} may be any possible string. The user has
    to make sure that the respective chromosome is indeed defined for
    the the track's genome. If not provided here, the constructor will
    try to build the chromosome information based on the available
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
    provided inputs, and eventually will fall back to the default value
    of \code{NA}.}

  \item{name}{Character scalar of the track's name. This is not really
    used and only exists fro completeness.}

   \item{\dots}{All additional parameters are ignored.}
}

\details{

  A track to conceptionally group other Gviz track objects into a meta
  track for the sole purpose of overlaying all the contained tracks with
  the same highlighting region as defined by the objects genomic
  ranges. During rendering the contained tracks will be treated as if
  they had been provided to the \code{plotTracks} function as individual
  objects.

}


\section{Slots}{

  \describe{

    \item{\code{trackList}:}{Object of class \code{"list"}, holding the
      subtrack objects.}

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
  of class \code{HighlightTrack}.


  \bold{\emph{Internal methods:}}

  \describe{

    \item{setStacks}{\code{signature(GdObject="HighlightTrack")}:
      recompute the stacks based on the available space and on the
      object's track items and stacking settings. This really just calls
      the \code{setStacks} methods for the contained tracks and only
      exists for dispatching reasons.

      \emph{Usage:}

      \code{setStacks(GdObject, ...)}

      \emph{Examples:}

      \describe{

      	\item{}{\code{Gviz:::setStacks(obj)}}

      }
    }


    \item{initialize}{\code{signature(.Object="HighlightTrack")}:
      initialize the object. }

     \item{subset}{\code{signature(x="HighlightTrack")}: subset all the
      contained tracks in an \code{HighlightTrack} by coordinates and sort
      if necessary.

      \emph{Usage:}

      \code{subset(x, ...)}

      \emph{Additional Arguments:}

      \describe{

	\item{}{\code{\dots}: additional arguments are passed on to the
	  next methods.}

      }

      \emph{Examples:}

      \describe{

	\item{}{\code{subset(obj)}}

      }
    }

    \item{length}{\code{signature(x="HighlightTrack")}: return the number
      of items in the track.

      \emph{Usage:}

      \code{length(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{length(obj)}}

      }
    }

  }

  \bold{\emph{Inherited methods:}}

  \describe{

    \item{[}{\code{signature(x="HighlightTrack", i="ANY", j="ANY",
	drop="ANY")}: subset the items in the \code{HighlightTrack}
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


    \item{chromosome}{\code{signature(GdObject="HighlightTrack")}:
      return the chromosome for which the track is defined.

      \emph{Usage:}

      \code{chromosome(GdObject)}

      \emph{Examples:}

      \describe{

	\item{}{\code{chromosome(obj)}}

      }
    }

    \item{chromosome<-}{\code{signature(GdObject="HighlightTrack")}:
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

    \item{start, end, width}{\code{signature(x="HighlightTrack")}: the
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

    \item{start<-, end<-, width<-}{\code{signature(x="HighlightTrack")}:
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

    \item{position}{\code{signature(GdObject="HighlightTrack")}: the
      arithmetic mean of the track item's coordionates, i.e.,
      \code{(end(obj)-start(obj))/2}.

      \emph{Usage:}

      \code{position(GdObject)}

      \emph{Examples:}

      \describe{

	\item{}{\code{position(obj)}}

      }
    }

    \item{genome}{\code{signature(x="HighlightTrack")}: return the track's genome.

      \emph{Usage:}

      \code{genome(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{genome(obj)}}

      }
    }

    \item{genome<-}{\code{signature(x="HighlightTrack")}: set the track's
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


    \item{range}{\code{signature(x="HighlightTrack")}: return the genomic
      coordinates for the track as an object of class
      \code{\linkS4class{IRanges}}.

      \emph{Usage:}

      \code{range(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{range(obj)}}

      }
    }

    \item{ranges}{\code{signature(x="HighlightTrack")}: return the genomic
      coordinates for the track along with all additional annotation
      information as an object of class \code{\linkS4class{GRanges}}.

      \emph{Usage:}

      \code{ranges(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{ranges(obj)}}

      }
    }

    \item{split}{\code{signature(x="HighlightTrack")}: split a
      \code{HighlightTrack} object by an appropriate factor vector (or
      another vector that can be coerced into one). The output of this
      operation is a list of objects of the same class as the input
      object, all inheriting from class \code{HighlightTrack}.

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


    \item{displayPars}{\code{signature(x="HighlightTrack",
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


    \item{displayPars}{\code{signature(x="HighlightTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.

      \emph{Examples:}

      \describe{

	\item{}{\code{displayPars(obj)}}

      }
    }

    \item{getPar}{\code{signature(x="HighlightTrack", name="character")}:
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

    \item{getPar}{\code{signature(x="HighlightTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Examples:}

      \describe{

	\item{}{\code{getPar(obj)}}

      }
    }

    \item{displayPars<-}{\code{signature(x="HighlightTrack",
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

    \item{setPar}{\code{signature(x="HighlightTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{HighlightTrack} class are
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

    \item{setPar}{\code{signature(x="HighlightTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{HighlightTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.


      \emph{Examples:}

      \describe{

	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }


    \item{names}{\code{signature(x="HighlightTrack")}: return the value of
      the \code{name} slot.

      \emph{Usage:}

      \code{names(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{names(obj)}}

      }
    }

    \item{names<-}{\code{signature(x="HighlightTrack", value="character")}:
      set the value of the \code{name} slot.

      \emph{Usage:}

      \code{names<-(x, value)}

      \emph{Examples:}

      \describe{

	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="HighlightTrack")}: return the
      coordinates from the internal image map.

      \emph{Usage:}

      \code{coords(ImageMap)}

      \emph{Examples:}

      \describe{

	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="HighlightTrack")}: return the tags from the
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
  \code{HighlightTrack} upon instantiation, unless one or more of them
  have already been set by one of the optional sub-class initializers,
  which always get precedence over these global defaults. See
  \code{\link{settings}} for details on setting graphical parameters for
  tracks.

  \describe{

    \item{}{\code{col="red"}: Integer or character scalar. The boder color for the highlighting regions.}

    \item{}{\code{fill="#FFE3E6"}: Integer or character scalar. The fill color for the highlighting regions.}

  }

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{HighlightTrack} objects.

    \describe{

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

  \code{\linkS4class{OverlayTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

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






























