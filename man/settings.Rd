\name{settings}
\Rdversion{1.1}
\alias{settings}

\title{Setting display parameters to control the look and feel of the plots}

\description{
  
  The genome track plots in this package are all highly customizable by
  means of so called 'display parameters'. This page highlights the use
  of these parameters and list all available settings for the different
  track classes.
  
}

\details{
  
  All of the package's track objects inherit the \code{dp} slot from the
  \code{\linkS4class{GdObject}} parent class, which is the main
  container to store an object's display parameters. Internally, the
  content of this slot has to be an object of class
  \code{\linkS4class{DisplayPars}}, but the user is usually not exposed
  to this low level implementation. Instead, there are two main
  interaction points, namely the individual object constructor functions
  and the final \code{\link{plotTracks}} function. In both cases, all
  additional arguments that are not caught by any of the formally
  defined function parameters are being interpreted as additional
  display parameters and are automatically added to the aforementioned
  slot. The main difference here is that display parameters that are
  passed on to the constructor function are specific for an individual
  track object, whereas those supplied to the \code{plotTracks} function
  will be applied to all the objects in the plotting list. Not all
  display parameters have an effect on the plotting of all track
  classes, and those will be silently ignored.
  
  One can query the available display parameters for a given class as
  well as their default values by calling the
  \code{\link{availableDisplayPars}} function, or by inspecting the man
  pages of the individual track classes. The structure of the classes
  defined in this package is hierarchical, and so are the available
  display parameters, i.e., all objects inherit the parameters defined
  in the commom \code{GdObject} parent class, and so on.
  
  Once a track object has been created, the display parameters are still
  open for modification. To this end, the \code{\link{displayPars}}
  replacement method is available for all objects inheriting from class
  \code{GdObject}. The method takes a named list of parameters as input,
  e.g.:
  
  \code{displayPars(foo) <- list(col="red", lwd=2)}
  
  In the same spirit, the currently set display parameters for the
  object \code{foo} can be inferred using the \code{displayPars} method
  directly, e.g.:
  
  \code{displayPars(foo)}
  
  For track objects inheriting from class
  \code{\linkS4class{AnnotationTrack}}, display parameters that are not
  formally defined in the class definition or in any of the parent
  classes are considered to be valid R color identifiers that are used
  to distinguish between different types of annotation features. For
  instance, the parameter 'miRNA' will be used to color all annotation
  features of class miRNA. The annotation types can be set in the
  constructor function of the track object via the \code{feature}
  argument. For most of the tracks that have been inferred from one of
  the online repositories, this classification will usually be
  downloaded along with the actual annotation data.
  
}


\section{Display Parameters}{ 

  \describe{ 

    \item{GenomeAxisTrack}{: 

      \describe{ 

        \item{}{\code{cex.id=0.7}: Numeric scalar. The text size for
          the optional range annotation.} 

        \item{}{\code{fontsize=10}: Numeric scalar. Font size for the
          axis annotation text in points.} 

        \item{}{\code{background.title="transparent"}: Character scalar.
          The background color for the title panel. Defaults to omit
          the background.} 

        \item{}{\code{fontcolor="#808080"}: Character scalar. The font
          color for the axis annotation text.} 

        \item{}{\code{size=NULL}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function. Defaults to the ideal size based on the other track
          settings.} 

        \item{}{\code{lwd=2}: Numeric scalar. The line width for the
          axis elementes.} 

        \item{}{\code{add35=FALSE}: Logical scalar. Add 3' to 5'
          direction indicators.} 

        \item{}{\code{showTitle=FALSE}: Logical scalar. Plot a title
          panel. Defaults to omit the title panel.} 

        \item{}{\code{add53=FALSE}: Logical scalar. Add 5' to 3'
          direction indicators.} 

        \item{}{\code{littleTicks=FALSE}: Logical scalar. Add more
          fine-grained tick marks.} 

        \item{}{\code{distFromAxis=1}: Numeric scalar. Control the
          distance of the axis annotation from the tick marks.} 

        \item{}{\code{col.id="white"}: Character scalar. The text color
          for the optional range annotation.} 

        \item{}{\code{labelPos="alternating"}: Character vector, one
          in "alternating", "revAlternating", "above" or "below". The
          vertical positioning of the axis labels.} 

        \item{}{\code{exponent=NULL}: Numeric scalar. The exponent for
          the axis coordinates, e.g., 3 means mb, 6 means gb, etc. The
          default is to automatically determine the optimal exponent.} 

        \item{}{\code{cex=0.8}: Numeric scalar. The overall font
          expansion factor for the axis annotation text.} 

        \item{}{\code{fill.range="cornsilk3"}: Character scalar. The
          fill color for highlighted regions on the axis.} 

        \item{}{\code{showId=FALSE}: Logical scalar. Show the optional
          range highlighting annotation.} 

        \item{}{\code{col.range="cornsilk4"}: Character scalar. The
          border color for highlighted regions on the axis.} 

        \item{}{\code{col="darkgray"}: Character scalar. The color for
          the axis lines and tickmarks.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text.} 

        \item{}{\code{fill="lightgray"}: Integer or character scalar.
          Default fill color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text.} 

      } 

    } 

    \item{DataTrack}{: 

      \describe{ 

        \item{}{\code{pch=20}: Integer scalar. The type of glyph used
          for plotting symbols.} 

        \item{}{\code{notch=FALSE}: Logical scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{cex=0.7}: Numeric scalar. The default pixel size
          for plotting symbols.} 

        \item{}{\code{box.width=NULL}: Numeric scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{stackedBars=TRUE}: Logical scalar. When there
          are several data groups, draw the histogram-type plots as
          stacked barplots or grouped side by side.} 

        \item{}{\code{levels.fos=NULL}: Numeric scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{gradient=c("#F7FBFF", "#DEEBF7", "#C6DBEF",
          "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C",
          "#08306B")}: Character vector. The base colors for the
          'gradient' plotting type.} 

        \item{}{\code{min.distance=0}: Numeric scalar. The mimimum
          distance in pixel below which to collapse ranges.} 

        \item{}{\code{degree=1}: Numeric scalar. Parameter controlling
          the loess calculation for smooth and mountain-type plots.
          See \code{\link{panel.loess}} for details.} 

        \item{}{\code{type="p"}: Character vector. The plot type, one
          or several in \code{c("p","l", "b", "a", "s", "g", "r", "S",
          "smooth", "histogram", "mountain", "h", "boxplot", "gradient",
          "heatmap")}. See 'Details' section in \code{\linkS4class{DataTrack}}
          for more information on the individual plotting types.} 

        \item{}{\code{notch.frac=0.5}: Numeric scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{ylim=NULL}: Numeric vector of length 2. The range
          of the y-axis scale.} 

        \item{}{\code{fill.histogram="lightgray"}: Character scalar.
          Fill color in histogram-type plots, defaults to the setting
          of \code{fill}.} 

        \item{}{\code{baseline=NULL}: Numeric scalar. Y-axis position
          of an optional baseline. This parameter has a special meaning
          for mountain-type plots, see the 'Details' section in
          \code{\linkS4class{DataTrack}} for more information.} 

        \item{}{\code{collapse=FALSE}: Logical scalar. Collapse
          overlapping ranges and aggregate the underlying data.} 

        \item{}{\code{size=NULL}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function. By default the size will be set automatically based
          on the selected plotting type.} 

        \item{}{\code{col.mountain=NULL}: Character scalar. Line color
          in mountain-type plots, defaults to the setting of \code{col}.} 

        \item{}{\code{span=0.2}: Numeric scalar. Parameter controlling
          the loess calculation for smooth and mountain-type plots.
          See \code{\link{panel.loess}} for details.} 

        \item{}{\code{transformation=NULL}: Function. Applied to the
          data matrix prior to plotting or when calling the \code{score}
          method. The function should accept exactly one input argument
          and its return value needs to be a numeric vector which can
          be coerced back into a data matrix of identical dimensionality
          as the input data.} 

        \item{}{\code{box.ratio=1}: Numeric scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{varwidth=FALSE}: Logical scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{family="symmetric"}: Character scalar. Parameter
          controlling the loess calculation for smooth and mountain-type
          plots. See \code{\link{panel.loess}} for details.} 

        \item{}{\code{separator=0}: Numeric scalar. Number of pixels
          used to separate individual samples in heatmap-type plots.} 

        \item{}{\code{ncolor=100}: Integer scalar. The number of colors
          for the 'gradient' plotting type} 

        \item{}{\code{stats=X[[44]]}: Function. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{amount=NULL}: Numeric scalar. Amount of jittering
          in xy-type plots. See \code{\link{panel.xyplot}} for details.} 

        \item{}{\code{na.rm=FALSE}: Boolean controlling whether to
          discard all NA values when plotting or to keep empty spaces for NAs} 

        \item{}{\code{lwd.baseline=NULL}: Numeric scalar. Line width
          of the optional baseline, defaults to the setting of \code{lwd}.} 

        \item{}{\code{windowSize=NULL}: Numeric scalar. The size of
          the running window when the value of \code{window} is negative.} 

        \item{}{\code{lty.mountain=NULL}: Character or numeric scalar.
          Line type in mountain-type plots, defaults to the setting of
          \code{lty}.} 

        \item{}{\code{jitter.x=FALSE}: Logical scalar. Toggle on
          jittering on the x axis in xy-type plots. See \code{\link{panel.xyplot}}
          for details.} 

        \item{}{\code{jitter.y=FALSE}: Logical scalar. Toggle off
          jittering on the y axis in xy-type plots. See \code{\link{panel.xyplot}}
          for details.} 

        \item{}{\code{groups=NULL}: Vector coercable to a factor.
          Optional sample grouping. See 'Details' section in
          \code{\linkS4class{DataTrack}} for further information.} 

        \item{}{\code{lty.baseline=NULL}: Character or numeric scalar.
          Line type of the optional baseline, defaults to the setting
          of \code{lty}.} 

        \item{}{\code{col.baseline=NULL}: Character scalar. Color for
          the optional baseline, defaults to the setting of \code{col}.} 

        \item{}{\code{evaluation=50}: Numeric scalar. Parameter
          controlling the loess calculation for smooth and mountain-type
          plots. See \code{\link{panel.loess}} for details.} 

        \item{}{\code{fill.mountain=c("#CCFFFF", "#FFCCFF")}: Character
          vector of length 2. Fill color in mountain-type plots.} 

        \item{}{\code{factor=0.5}: Numeric scalar. Factor to control
          amount of jittering in xy-type plots. See \code{\link{panel.xyplot}}
          for details.} 

        \item{}{\code{do.out=TRUE}: Logical scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{col.histogram="#808080"}: Character scalar. Line
          color in histogram-type plots.} 

        \item{}{\code{aggregation="mean"}: Function or character scalar.
          Used to aggregate values in windows or for collapsing
          overlapping items. The function has to accept a numeric vector
          as a single input parameter and has to return a numeric scalar
          with the aggregated value. Alternatively, one of the predefined
          options \code{mean}, \code{median} \code{sum}, \code{min},
          \code{max} or \code{extreme} can be supplied as a character
          scalar. Defaults to \code{mean}.} 

        \item{}{\code{coef=1.5}: Numeric scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{lwd.mountain=NULL}: Numeric scalar. Line width
          in mountain-type plots, defaults to the setting of \code{lwd}.} 

        \item{}{\code{window=NULL}: Numeric or character scalar.
          Aggregate the rows values of the data matrix to \code{window}
          equally sized slices on the data range using the method
          defined in \code{aggregation}. If negative, apply a running
          window of size \code{windowSize} using the same aggregation
          method. Alternatively, the special value \code{auto} causes
          the function to determine the optimal window size to avoid
          overplotting.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text.} 

        \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
          for all plotting elements, unless there is a more specific
          control defined elsewhere.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{showTitle=TRUE}: Boolean controlling whether to
          plot a title panel. Although this can be set individually
          for each track, in multi-track plots as created by
          \code{\link{plotTracks}} there will still be an empty
          placeholder in case any of the other tracks include a title.
          The same holds true for axes. Note that the the title panel
          background color could be set to transparent in order to
          completely hide it.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text.} 

        \item{}{\code{fill="lightgray"}: Integer or character scalar.
          Default fill color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{col="#0080FF"}: Integer or character scalar.
          Default line color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text.} 

      } 

    } 

    \item{IdeogramTrack}{: 

      \describe{ 

        \item{}{\code{fontsize=10}: Numeric scalar. The font size for
          the chromosome name text.} 

        \item{}{\code{background.title="transparent"}: Character scalar.
          The background color for the title panel. Defaults to omit
          the background.} 

        \item{}{\code{fontcolor="#808080"}: Character scalar. The font
          color for the chromosome name text.} 

        \item{}{\code{size=NULL}: Numeric scalar. The relative size of
          the track. Defaults to automatic size setting. Can be overridden
          in the \code{\link{plotTracks}} function.} 

        \item{}{\code{showTitle=FALSE}: Logical scalar. Plot a title
          panel. Defaults to omit the title panel.} 

        \item{}{\code{cex=0.8}: Numeric scalar. The overall font
          expansion factor for the chromosome name text.} 

        \item{}{\code{fill="#FFE3E6"}: Character scalar. The fill color
          used for the highlighting of the currently displayed genomic
          region.} 

        \item{}{\code{showId=TRUE}: Logical scalar. Indicate the
          chromosome name next to the ideogram.} 

        \item{}{\code{col="red"}: Character scalar. The border color
          used for the highlighting of the currently displayed genomic
          region.} 

        \item{}{\code{bevel=0.45}: Numeric scalar, between 0 and 1.
          The level of smoothness for the two ends of the ideogram.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
          for all plotting elements, unless there is a more specific
          control defined elsewhere.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text.} 

      } 

    } 

    \item{AnnotationTrack}{: 

      \describe{ 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{showFeatureId=FALSE}: Logical scalar. Control
          whether to plot the individual track item identifiers.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          item identifiers.} 

        \item{}{\code{showOverplotting=FALSE}: Logical scalar. Use a
          color gradient to show the amount of overplotting for collapsed
          items. This implies that \code{collapse==TRUE}} 

        \item{}{\code{lex=1}: Numeric scalar. The line expansion factor
          for all track items. This is also used to connect grouped
          items. See \code{\link{grouping}} for details.} 

        \item{}{\code{fontcolor="white"}: Character or integer scalar.
          The font color for item identifiers.} 

        \item{}{\code{fontfamily="sans"}: Character scalar. The font
          family for item identifiers.} 

        \item{}{\code{size=1}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function.} 

        \item{}{\code{lwd=1}: Integer scalar. The line width for all
          track items. This is also used to connect grouped items. See
          \code{\link{grouping}} for details.} 

        \item{}{\code{lty="solid"}: Character or integer scalar. The
          line type for all track items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{shape="arrow"}: Character scalar. The shape in
          which to display the track items. Currently only \code{box},
          \code{arrow}, \code{ellipse}, and \code{smallArrow} are
          implemented.} 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{fontface=1}: Integer scalar. The font face for
          item identifiers.} 

        \item{}{\code{fill="lightblue"}: Character or integer scalar.
          The fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{showId=FALSE}: Logical scalar. Control whether
          to annotate individual groups.} 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for item identifiers.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{col="transparent"}: Character or integer scalar.
          The border color for all track items.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{showTitle=TRUE}: Boolean controlling whether to
          plot a title panel. Although this can be set individually
          for each track, in multi-track plots as created by
          \code{\link{plotTracks}} there will still be an empty
          placeholder in case any of the other tracks include a title.
          The same holds true for axes. Note that the the title panel
          background color could be set to transparent in order to
          completely hide it.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

      } 

    } 

    \item{GeneRegionTrack}{: 

      \describe{ 

        \item{}{\code{showExonId=FALSE}: Logical scalar. Control whether
          to plot the individual exon identifiers.} 

        \item{}{\code{shape=c("smallArrow", "box")}: Character scalar.
          The shape in which to display the track items. Currently only
          \code{box}, \code{arrow}, \code{ellipse}, and \code{smallArrow}
          are implemented.} 

        \item{}{\code{geneSymbols=TRUE}: Logical scalar. Use human-readable
          gene symbols or gene IDs for the transcript annotation.} 

        \item{}{\code{fill="orange"}: Character or integer scalar. The
          fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

      } 

      \bold{\emph{Inherited from class AnnotationTrack:}} 

      \describe{ 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{showFeatureId=FALSE}: Logical scalar. Control
          whether to plot the individual track item identifiers.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          item identifiers.} 

        \item{}{\code{showOverplotting=FALSE}: Logical scalar. Use a
          color gradient to show the amount of overplotting for collapsed
          items. This implies that \code{collapse==TRUE}} 

        \item{}{\code{lex=1}: Numeric scalar. The line expansion factor
          for all track items. This is also used to connect grouped
          items. See \code{\link{grouping}} for details.} 

        \item{}{\code{fontcolor="white"}: Character or integer scalar.
          The font color for item identifiers.} 

        \item{}{\code{fontfamily="sans"}: Character scalar. The font
          family for item identifiers.} 

        \item{}{\code{size=1}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function.} 

        \item{}{\code{lwd=1}: Integer scalar. The line width for all
          track items. This is also used to connect grouped items. See
          \code{\link{grouping}} for details.} 

        \item{}{\code{lty="solid"}: Character or integer scalar. The
          line type for all track items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{fontface=1}: Integer scalar. The font face for
          item identifiers.} 

        \item{}{\code{showId=FALSE}: Logical scalar. Control whether
          to annotate individual groups.} 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for item identifiers.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{col="transparent"}: Character or integer scalar.
          The border color for all track items.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{showTitle=TRUE}: Boolean controlling whether to
          plot a title panel. Although this can be set individually
          for each track, in multi-track plots as created by
          \code{\link{plotTracks}} there will still be an empty
          placeholder in case any of the other tracks include a title.
          The same holds true for axes. Note that the the title panel
          background color could be set to transparent in order to
          completely hide it.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

      } 

    } 

    \item{BiomartGeneRegionTrack}{: 

      \describe{ 

        \item{}{\code{scRNA_pseudogene="darkorange2"}: Character or
          integer scalar. Fill color for annotation objects of type
          'scRNA_pseudogene'.} 

        \item{}{\code{miRNA_pseudogene="cornsilk"}: Character or integer
          scalar. Fill color for annotation objects of type
          'miRNA_pseudogene'.} 

        \item{}{\code{D_segment="lightblue"}: Character or integer
          scalar. Fill color for annotation objects of type 'C_segment'.} 

        \item{}{\code{scRNA="darkorange"}: Character or integer scalar.
          Fill color for annotation objects of type 'scRNA'.} 

        \item{}{\code{J_segment="dodgerblue2"}: Character or integer
          scalar. Fill color for annotation objects of type 'C_segment'.} 

        \item{}{\code{Mt_tRNA_pseudogene="darkgoldenrod1"}: Character
          or integer scalar. Fill color for annotation objects of type
          'Mt_tRNA_pseudogene'.} 

        \item{}{\code{Mt_rRNA="yellow"}: Character or integer scalar.
          Fill color for annotation objects of type 'Mt_rRNA'.} 

        \item{}{\code{tRNA_pseudogene="antiquewhite3"}: Character or
          integer scalar. Fill color for annotation objects of type
          'tRNA_pseudogene'.} 

        \item{}{\code{rRNA_pseudogene="darkolivegreen"}: Character or
          integer scalar. Fill color for annotation objects of type
          'rRNA_pseudogene'.} 

        \item{}{\code{miRNA="cornflowerblue"}: Character or integer
          scalar. Fill color for annotation objects of type 'L_segment'.} 

        \item{}{\code{pseudogene="brown1"}: Character or integer scalar.
          Fill color for annotation objects of type 'pseudogene'.} 

        \item{}{\code{rRNA="darkolivegreen1"}: Character or integer
          scalar. Fill color for annotation objects of type 'rRNA'.} 

        \item{}{\code{protein_coding="gold4"}: Character or integer
          scalar. Fill color for annotation objects of type 'protein_coding'.} 

        \item{}{\code{Mt_tRNA="darkgoldenrod"}: Character or integer
          scalar. Fill color for annotation objects of type 'Mt_tRNA'.} 

        \item{}{\code{retrotransposed="blueviolet"}: Character or
          integer scalar. Fill color for annotation objects of type
          'retrotransposed'.} 

        \item{}{\code{snRNA="coral"}: Character or integer scalar. Fill
          color for annotation objects of type 'snRNA'.} 

        \item{}{\code{V_segment="aquamarine"}: Character or integer
          scalar. Fill color for annotation objects of type 'V_segment'.} 

        \item{}{\code{snRNA_pseudogene="coral3"}: Character or integer
          scalar. Fill color for annotation objects of type
          'snRNA_pseudogene'.} 

        \item{}{\code{misc_RNA="cornsilk3"}: Character or integer
          scalar. Fill color for annotation objects of type 'misc_RNA'.} 

        \item{}{\code{misc_RNA_pseudogene="cornsilk4"}: Character or
          integer scalar. Fill color for annotation objects of type
          'misc_RNA_pseudogene'.} 

        \item{}{\code{snoRNA_pseudogene="cyan2"}: Character or integer
          scalar. Fill color for annotation objects of type
          'snoRNA_pseudogene'.} 

        \item{}{\code{snoRNA="cyan"}: Character or integer scalar. Fill
          color for annotation objects of type 'snoRNA'.} 

        \item{}{\code{C_segment="burlywood4"}: Character or integer
          scalar. Fill color for annotation objects of type 'C_segment'.} 

      } 

      \bold{\emph{Inherited from class GeneRegionTrack:}} 

      \describe{ 

        \item{}{\code{showExonId=FALSE}: Logical scalar. Control whether
          to plot the individual exon identifiers.} 

        \item{}{\code{shape=c("smallArrow", "box")}: Character scalar.
          The shape in which to display the track items. Currently only
          \code{box}, \code{arrow}, \code{ellipse}, and \code{smallArrow}
          are implemented.} 

        \item{}{\code{geneSymbols=TRUE}: Logical scalar. Use human-readable
          gene symbols or gene IDs for the transcript annotation.} 

        \item{}{\code{fill="orange"}: Character or integer scalar. The
          fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

      } 

      \bold{\emph{Inherited from class AnnotationTrack:}} 

      \describe{ 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{showFeatureId=FALSE}: Logical scalar. Control
          whether to plot the individual track item identifiers.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          item identifiers.} 

        \item{}{\code{showOverplotting=FALSE}: Logical scalar. Use a
          color gradient to show the amount of overplotting for collapsed
          items. This implies that \code{collapse==TRUE}} 

        \item{}{\code{lex=1}: Numeric scalar. The line expansion factor
          for all track items. This is also used to connect grouped
          items. See \code{\link{grouping}} for details.} 

        \item{}{\code{fontcolor="white"}: Character or integer scalar.
          The font color for item identifiers.} 

        \item{}{\code{fontfamily="sans"}: Character scalar. The font
          family for item identifiers.} 

        \item{}{\code{size=1}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function.} 

        \item{}{\code{lwd=1}: Integer scalar. The line width for all
          track items. This is also used to connect grouped items. See
          \code{\link{grouping}} for details.} 

        \item{}{\code{lty="solid"}: Character or integer scalar. The
          line type for all track items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{fontface=1}: Integer scalar. The font face for
          item identifiers.} 

        \item{}{\code{showId=FALSE}: Logical scalar. Control whether
          to annotate individual groups.} 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for item identifiers.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{col="transparent"}: Character or integer scalar.
          The border color for all track items.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{showTitle=TRUE}: Boolean controlling whether to
          plot a title panel. Although this can be set individually
          for each track, in multi-track plots as created by
          \code{\link{plotTracks}} there will still be an empty
          placeholder in case any of the other tracks include a title.
          The same holds true for axes. Note that the the title panel
          background color could be set to transparent in order to
          completely hide it.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{collapse=TRUE}: Boolean controlling wether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

      } 

    } 

    \item{AlignedReadTrack}{: 

      \describe{ 

        \item{}{\code{detail="coverage"}: the amount of detail to plot
          the data. Either \code{coverage} to show the coverage only,
          or \code{reads} to show individual reads. For large data sets
          the latter can be very inefficient. Please note that \code{reads}
          is only available when the object has been created with option
          \code{coverageOnly=FALSE}.} 

        \item{}{\code{type="histogram"}: the plot type, one or several
          in \code{c("p","l", "b", "a", "s", "g", "r", "S", "smooth",
          "histogram", "mountain", "h", "boxplot", "gradient", "heatmap")}.
          See	the 'Details' section in \code{\linkS4class{DataTrack}}
          for more information on the individual plotting types.} 

        \item{}{\code{size=NULL}: the relative size of the track.
          Defaults to size selection based on the underlying data. Can
          be overridden in the \code{\link{plotTracks}} function.} 

        \item{}{\code{collapse=FALSE}: collapse overlapping ranges and
          aggregate the underlying data.} 

        \item{}{\code{fill="#0080ff"}: the fill color for the coverage
          indicator.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text.} 

        \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
          for all plotting elements, unless there is a more specific
          control defined elsewhere.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panels.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{col.grid="#808080"}: Integer or character scalar.
          Default line color for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{showTitle=TRUE}: Boolean controlling whether to
          plot a title panel. Although this can be set individually
          for each track, in multi-track plots as created by
          \code{\link{plotTracks}} there will still be an empty
          placeholder in case any of the other tracks include a title.
          The same holds true for axes. Note that the the title panel
          background color could be set to transparent in order to
          completely hide it.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{col.line=NULL}: Integer or character scalar.
          Default colors for plot lines. Usually the same as the global
          \code{col} parameter.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{col.symbol=NULL}: Integer or character scalar.
          Default colors for plot symbols. Usually the same as the
          global \code{col} parameter.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{cex=1}: Numeric scalar. The overall font expansion
          factor for all text.} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is computed based on the available space.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text.} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.frame="lightgray"}: Integer or character
          scalar. The line color used for the panel frame, if
          \code{frame==TRUE}} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{col="#0080FF"}: Integer or character scalar.
          Default line color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{col.title="white"}: Integer or character scalar.
          The font color for the title panels.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text.} 

      } 

    } 

  } 

}


\author{
  Florian Hahne
}
\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{DataTrack}}

  \code{\linkS4class{DisplayPars}}

  \code{\linkS4class{GdObject}}

  \code{\link{availableDisplayPars}}

  \code{\link{collapsing}}

  \code{\link{displayPars}}

  \code{\link{grouping}}

  \code{\link{panel.bwplot}}

  \code{\link{panel.grid}}

  \code{\link{panel.loess}}

  \code{\link{panel.xyplot}}

  \code{\link{plotTracks}}

}












