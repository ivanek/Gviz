\name{OverlayTrack-class}
\Rdversion{1.1}
\docType{class}
\alias{OverlayTrack-class}
\alias{OverlayTrack}
\alias{drawGD,OverlayTrack-method}
\alias{initialize,OverlayTrack-method}
\alias{setStacks,OverlayTrack-method}
\alias{subset,OverlayTrack-method}
\alias{show,OverlayTrack-method}
\alias{length,OverlayTrack-method}
\alias{consolidateTrack,OverlayTrack-method}
\alias{chromosome<-,OverlayTrack-method}
\alias{displayPars<-,OverlayTrack,list-method}

\title{OverlayTrack class and methods}

\description{

  A container for other track objects from the Gviz package that allows
  for overlays of their content on the same region of the plot.

}

\section{Objects from the Class}{

 Objects can be created using the constructor function
  \code{OverlayTrack}.

}

\usage{

OverlayTrack(trackList=list(), name="OverlayTrack", ...)

}

\arguments{

  \item{trackList}{A list of Gviz track objects that all have to inherit
    from class \code{GdObject}.}

  \item{name}{Character scalar of the track's name. This is not really
    used and only exists fro completeness.}

  \item{\dots}{All additional parameters are ignored.}


}

\details{

  A track to conceptionally group other Gviz track objects into a meta
  track in order to merge them into a single overlay visualization. Only
  the first track in the supplied list will be inferred when setting up
  the track title and axis, for all the other tracks only the panel
  content is plotted.

}


\section{Slots}{

  \describe{

    \item{\code{trackList}:}{Object of class \code{"list"}, holding the
      subtrack objects.}

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

  Class \code{"\linkS4class{GdObject}"}, directly.

}

\section{Methods}{

  In the following code chunks, \code{obj} is considered to be an object
  of class \code{OverlayTrack}.


  \bold{\emph{Internal methods:}}

  \describe{

    \item{setStacks}{\code{signature(GdObject="OverlayTrack")}:
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


    \item{initialize}{\code{signature(.Object="OverlayTrack")}:
      initialize the object. }


    \item{length}{\code{signature(x="OverlayTrack")}: return the number
      of items in the track.

      \emph{Usage:}

      \code{length(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{length(obj)}}

      }
    }


    \item{subset}{\code{signature(x="OverlayTrack")}: subset all the
      contained tracks in an \code{OverlayTrack} by coordinates and sort
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

  }

  \bold{\emph{Inherited methods:}}

  \describe{

    \item{displayPars}{\code{signature(x="OverlayTrack",
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


    \item{displayPars}{\code{signature(x="OverlayTrack",
      name="missing")}: list the value of all available display
      parameters. See \code{\link{settings}} for details on display
      parameters and customization.

      \emph{Examples:}

      \describe{

	\item{}{\code{displayPars(obj)}}

      }
    }

    \item{getPar}{\code{signature(x="OverlayTrack", name="character")}:
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

    \item{getPar}{\code{signature(x="OverlayTrack", name="missing")}:
      alias for the \code{displayPars} method. See
      \code{\link{settings}} for details on display parameters and
      customization.

      \emph{Examples:}

      \describe{

	\item{}{\code{getPar(obj)}}

      }
    }

    \item{displayPars<-}{\code{signature(x="OverlayTrack",
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

    \item{setPar}{\code{signature(x="OverlayTrack", value="character")}:
      set the single display parameter \code{name} to \code{value}. Note
      that display parameters in the \code{OverlayTrack} class are
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

    \item{setPar}{\code{signature(x="OverlayTrack", value="list")}: set
      display parameters by the values of the named list in
      \code{value}.  Note that display parameters in the
      \code{OverlayTrack} class are pass-by-reference, so no
      re-assignmnet to the symbol \code{obj} is necessary. See
      \code{\link{settings}} for details on display parameters and
      customization.


      \emph{Examples:}

      \describe{

	\item{}{\code{setPar(obj, list(col="red", lwd=2))}}

      }
    }


    \item{names}{\code{signature(x="OverlayTrack")}: return the value of
      the \code{name} slot.

      \emph{Usage:}

      \code{names(x)}

      \emph{Examples:}

      \describe{

	\item{}{\code{names(obj)}}

      }
    }

    \item{names<-}{\code{signature(x="OverlayTrack", value="character")}:
      set the value of the \code{name} slot.

      \emph{Usage:}

      \code{names<-(x, value)}

      \emph{Examples:}

      \describe{

	\item{}{\code{names(obj) <- "foo"}}

      }
    }

    \item{coords}{\code{signature(ImageMap="OverlayTrack")}: return the
      coordinates from the internal image map.

      \emph{Usage:}

      \code{coords(ImageMap)}

      \emph{Examples:}

      \describe{

	\item{}{\code{coords(obj)}}

      }
    }

    \item{tags}{\code{signature(x="OverlayTrack")}: return the tags from the
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

  No formal display parameters are defined for objects of class
  \code{OverlayTrack}.

  Additional display parameters are being inherited from the respective
  parent classes. Note that not all of them may have an effect on the
  plotting of \code{OverlayTrack} objects.

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

  \code{\linkS4class{HighlightTrack}}

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






























