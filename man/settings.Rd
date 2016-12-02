\name{settings}
\Rdversion{1.1}
\alias{settings}
\alias{addScheme}
\alias{getScheme}

\title{Setting display parameters to control the look and feel of the plots}

\description{
  
  The genome track plots in this package are all highly customizable by
  means of so called 'display parameters'. This page highlights the use
  of these parameters and list all available settings for the different
  track classes.
  
}

\usage{

addScheme(scheme, name)

getScheme(name=getOption("Gviz.scheme"))

}

\arguments{

  \item{scheme}{A named nested list of display parameters, where the
    first level of nesting represents Gviz track object classes, and the
    second level of nesting represents parameters.}
  
  \item{name}{A character scalar with the scheme name.}
  
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

  Users might find themselves changing the same parameters over and over
  again, and it would make sense to register these modifications in a
  central location once and for all. To this end the Gviz package
  supports display parameter schemes. A scheme is essentially just a
  bunch of nested named lists, where the names on the first level of
  nesting should correspond to track class names, and the names on the
  second level to the display parameters to set. The currently active
  schmeme can be changed by setting the global option
  \code{Gviz.scheme}, and a new scheme can be registered by using the
  \code{addScheme} function, providing both the list and the name for
  the new scheme. The \code{getScheme} function is useful to get the
  current scheme as a list structure, for instance to use as a skeleton
  for your own custom scheme.

  In order to make these settings persitant across R sessions one can
  create one or several schemes in the global environment in the special
  object \code{.GvizSchemes}, for instance by putting the necessary code
  in the \code{.Rprofile} file. This object needs to be a named list of
  schemes, and it will be collected when the Givz package loads. Its
  content is then automatically added to the collection of available
  schemes.
  
  Please note that because display parameters are stored with the track
  objects, a scheme change only has an effect on those objects that are
  created after the change has taken place.

  
}


\section{Display Parameters}{ 

  \describe{ 

    \item{GenomeAxisTrack}{: 

      \describe{ 

        \item{}{\code{add35=FALSE}: Logical scalar. Add 3' to 5'
          direction indicators.} 

        \item{}{\code{add53=FALSE}: Logical scalar. Add 5' to 3'
          direction indicators.} 

        \item{}{\code{background.title="transparent"}: Character scalar.
          The background color for the title panel. Defaults to omit
          the background.} 

        \item{}{\code{cex.id=0.7}: Numeric scalar. The text size for
          the optional range annotation.} 

        \item{}{\code{cex=0.8}: Numeric scalar. The overall font
          expansion factor for the axis annotation text.} 

        \item{}{\code{col.border.title="transparent"}: Integer or
          character scalar. The border color for the title panels.} 

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{col.id="white"}: Character scalar. The text color
          for the optional range annotation.} 

        \item{}{\code{col.range="cornsilk4"}: Character scalar. The
          border color for highlighted regions on the axis.} 

        \item{}{\code{distFromAxis=1}: Numeric scalar. Control the
          distance of the axis annotation from the tick marks.} 

        \item{}{\code{exponent=NULL}: Numeric scalar. The exponent for
          the axis coordinates, e.g., 3 means mb, 6 means gb, etc. The
          default is to automatically determine the optimal exponent.} 

        \item{}{\code{fill.range="cornsilk3"}: Character scalar. The
          fill color for highlighted regions on the axis.} 

        \item{}{\code{fontcolor="#808080"}: Character scalar. The font
          color for the axis annotation text.} 

        \item{}{\code{fontsize=10}: Numeric scalar. Font size for the
          axis annotation text in points.} 

        \item{}{\code{labelPos="alternating"}: Character vector, one
          in "alternating", "revAlternating", "above" or "below". The
          vertical positioning of the axis labels. If \code{scale} is
          not \code{NULL}, the possible values are "above", "below"
          and "beside".} 

        \item{}{\code{littleTicks=FALSE}: Logical scalar. Add more
          fine-grained tick marks.} 

        \item{}{\code{lwd=2}: Numeric scalar. The line width for the
          axis elementes.} 

        \item{}{\code{scale=NULL}: Numeric scalar. If not \code{NULL}
          a small scale is drawn instead of the full axis, if the value
          is between 0 and 1 it is interpreted as a fraction of the
          current plotting region, otherwise as an absolute length
          value in genomic coordinates.} 

        \item{}{\code{showId=FALSE}: Logical scalar. Show the optional
          range highlighting annotation.} 

        \item{}{\code{showTitle=FALSE}: Logical scalar. Plot a title
          panel. Defaults to omit the title panel.} 

        \item{}{\code{size=NULL}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function. Defaults to the ideal size based on the other track
          settings.} 

        \item{}{\code{col="darkgray"}: Character scalar. The color for
          the axis lines and tickmarks.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

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
          The border color for the title panels} 

        \item{}{\code{collapse=TRUE}: Boolean controlling whether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{fill="lightgray"}: Integer or character scalar.
          Default fill color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text, unless a more specific definition exists.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text, unless a more specific definition
          exists.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text, unless a more specific definition exists.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

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

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

        \item{}{\code{rotation=0}: The rotation angle for all text
          unless a more specific definiton exists.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

      } 

    } 

    \item{DataTrack}{: 

      \describe{ 

        \item{}{\code{aggregateGroups=FALSE}: Logical scalar. Aggregate
          the values within a sample group using the aggregation
          funnction specified in the \code{aggregation} parameter.} 

        \item{}{\code{aggregation="mean"}: Function or character scalar.
          Used to aggregate values in windows or for collapsing
          overlapping items. The function has to accept a numeric vector
          as a single input parameter and has to return a numeric scalar
          with the aggregated value. Alternatively, one of the predefined
          options \code{mean}, \code{median} \code{sum}, \code{min},
          \code{max} or \code{extreme} can be supplied as a character
          scalar. Defaults to \code{mean}.} 

        \item{}{\code{alpha.confint=0.3}: Numeric scalar. The transparency
          for the confidence intervalls in confint-type plots.} 

        \item{}{\code{amount=NULL}: Numeric scalar. Amount of jittering
          in xy-type plots. See \code{\link{panel.xyplot}} for details.} 

        \item{}{\code{baseline=NULL}: Numeric scalar. Y-axis position
          of an optional baseline. This parameter has a special meaning
          for mountain-type and polygon-type plots, see the 'Details'
          section in \code{\linkS4class{DataTrack}} for more information.} 

        \item{}{\code{box.legend=FALSE}: Logical scalar. Draw a box
          around a legend.} 

        \item{}{\code{box.ratio=1}: Numeric scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{box.width=NULL}: Numeric scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{grid=FALSE}: Logical vector. Draw a line grid
          under the track content.} 

        \item{}{\code{cex.legend=0.8}: Numeric scalar. The size factor
          for the legend text.} 

        \item{}{\code{cex.sampleNames=NULL}: Numeric scalar. The size
          factor for the sample names text in heatmap or horizon plots.
          Defaults to an automatic setting.} 

        \item{}{\code{cex=0.7}: Numeric scalar. The default pixel size
          for plotting symbols.} 

        \item{}{\code{coef=1.5}: Numeric scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{col.baseline=NULL}: Character scalar. Color for
          the optional baseline, defaults to the setting of \code{col}.} 

        \item{}{\code{col.confint=NA}: Character vector. Border colors
          for the confidence intervals for confint-type plots.} 

        \item{}{\code{col.histogram="#808080"}: Character scalar. Line
          color in histogram-type plots.} 

        \item{}{\code{col.horizon=NA}: The line color for the segments
          in the \code{horizon}-type plot. See \code{\link{horizonplot}}
          for details.} 

        \item{}{\code{col.mountain=NULL}: Character scalar. Line color
          in mountain-type and polygon-type plots, defaults to the
          setting of \code{col}.} 

        \item{}{\code{col.sampleNames="white"}: Character or integer
          scalar. The color used for the sample names in heatmap plots.} 

        \item{}{\code{col=c("#0080ff", "#ff00ff", "darkgreen", "#ff0000",
          "orange", "#00ff00", "brown")}: Character or integer vector.
          The color used for all line and symbol elements, unless there
          is a more specific control defined elsewhere. Unless
          \code{groups} are specified, only the first color in the
          vector is usually regarded.} 

        \item{}{\code{collapse=FALSE}: Logical scalar. Collapse
          overlapping ranges and aggregate the underlying data.} 

        \item{}{\code{degree=1}: Numeric scalar. Parameter controlling
          the loess calculation for smooth and mountain-type plots.
          See \code{\link{panel.loess}} for details.} 

        \item{}{\code{do.out=TRUE}: Logical scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{evaluation=50}: Numeric scalar. Parameter
          controlling the loess calculation for smooth and mountain-type
          plots. See \code{\link{panel.loess}} for details.} 

        \item{}{\code{factor=0.5}: Numeric scalar. Factor to control
          amount of jittering in xy-type plots. See \code{\link{panel.xyplot}}
          for details.} 

        \item{}{\code{family="symmetric"}: Character scalar. Parameter
          controlling the loess calculation for smooth and mountain-type
          plots. See \code{\link{panel.loess}} for details.} 

        \item{}{\code{fill.confint=NULL}: Character vector. Fill colors
          for the confidence intervals for confint-type plots.} 

        \item{}{\code{fill.histogram=NULL}: Character scalar. Fill
          color in histogram-type plots, defaults to the setting of
          \code{fill}.} 

        \item{}{\code{fill.horizon=c("#B41414", "#E03231", "#F7A99C",
          "#9FC8DC", "#468CC8", "#0165B3")}: The fill colors for the
          segments in the \code{horizon}-type plot. This should be a
          vector of length six, where the first three entries are the
          colors for positive changes, and the latter three entries
          are the colors for negative changes. Defaults to a red-blue
          color scheme. See \code{\link{horizonplot}} for details.} 

        \item{}{\code{fill.mountain=c("#CCFFFF", "#FFCCFF")}: Character
          vector of length 2. Fill color in mountain-type and polygon-type
          plots.} 

        \item{}{\code{fontface.legend=NULL}: Integer or character
          scalar. The font face for the legend text.} 

        \item{}{\code{fontfamily.legend=NULL}: Integer or character
          scalar. The font family for the legend text.} 

        \item{}{\code{fontsize.legend=NULL}: Numeric scalar. The pixel
          size for the legend text.} 

        \item{}{\code{fontcolor.legend="#808080"}: Integer or character
          scalar. The font color for the legend text.} 

        \item{}{\code{gradient=c("#F7FBFF", "#DEEBF7", "#C6DBEF",
          "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C",
          "#08306B")}: Character vector. The base colors for the
          \code{gradient} plotting type or the \code{heatmap} type with
          a single group. When plotting heatmaps with more than one
          group, the \code{col} parameter can be used to control the
          group color scheme, however the gradient will always be from
          white to 'col' and thus does not offer as much flexibility
          as this \code{gradient} parameter.} 

        \item{}{\code{groups=NULL}: Vector coercable to a factor.
          Optional sample grouping. See 'Details' section in
          \code{\linkS4class{DataTrack}} for further information.} 

        \item{}{\code{horizon.origin=0}: The baseline relative to which
          changes are indicated on the \code{horizon}-type plot. See
          \code{\link{horizonplot}} for details.} 

        \item{}{\code{horizon.scale=NULL}: The scale for each of the
          segments in the \code{horizon}-type plot. Defaults to 1/3 of
          the absolute data range. See \code{\link{horizonplot}} for details.} 

        \item{}{\code{jitter.x=FALSE}: Logical scalar. Toggle on
          jittering on the x axis in xy-type plots. See \code{\link{panel.xyplot}}
          for details.} 

        \item{}{\code{jitter.y=FALSE}: Logical scalar. Toggle off
          jittering on the y axis in xy-type plots. See \code{\link{panel.xyplot}}
          for details.} 

        \item{}{\code{levels.fos=NULL}: Numeric scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{legend=TRUE}: Boolean triggering the addition of
          a legend to the track to indicate groups. This only has an
          effect if at least two groups are present.} 

        \item{}{\code{lineheight.legend=NULL}: Numeric scalar. The line
          height for the legend text.} 

        \item{}{\code{lty.baseline=NULL}: Character or numeric scalar.
          Line type of the optional baseline, defaults to the setting
          of \code{lty}.} 

        \item{}{\code{lty.mountain=NULL}: Character or numeric scalar.
          Line type in mountain-type and polygon-type plots, defaults
          to the setting of \code{lty}.} 

        \item{}{\code{lwd.baseline=NULL}: Numeric scalar. Line width
          of the optional baseline, defaults to the setting of \code{lwd}.} 

        \item{}{\code{lwd.mountain=NULL}: Numeric scalar. Line width
          in mountain-type and polygon-type plots, defaults to the
          setting of \code{lwd}.} 

        \item{}{\code{min.distance=0}: Numeric scalar. The mimimum
          distance in pixel below which to collapse ranges.} 

        \item{}{\code{na.rm=FALSE}: Boolean controlling whether to
          discard all NA values when plotting or to keep empty spaces for NAs} 

        \item{}{\code{ncolor=100}: Integer scalar. The number of colors
          for the 'gradient' plotting type} 

        \item{}{\code{notch.frac=0.5}: Numeric scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{notch=FALSE}: Logical scalar. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{pch=20}: Integer scalar. The type of glyph used
          for plotting symbols.} 

        \item{}{\code{separator=0}: Numeric scalar. Number of pixels
          used to separate individual samples in heatmap- and horizon-type
          plots.} 

        \item{}{\code{showColorBar=TRUE}: Boolean. Indicate the data
          range color mapping in the axis for 'heatmap' or 'gradient' types.} 

        \item{}{\code{showSampleNames=FALSE}: Boolean. Display the
          names of the individual samples in a heatmap or a horizon plot.} 

        \item{}{\code{size=NULL}: Numeric scalar. The relative size of
          the track. Can be overridden in the \code{\link{plotTracks}}
          function. By default the size will be set automatically based
          on the selected plotting type.} 

        \item{}{\code{span=0.2}: Numeric scalar. Parameter controlling
          the loess calculation for smooth and mountain-type plots.
          See \code{\link{panel.loess}} for details.} 

        \item{}{\code{stackedBars=TRUE}: Logical scalar. When there
          are several data groups, draw the histogram-type plots as
          stacked barplots or grouped side by side.} 

        \item{}{\code{stats=X[[i]]}: Function. Parameter controlling
          the boxplot appearance. See \code{\link{panel.bwplot}} for details.} 

        \item{}{\code{transformation=NULL}: Function. Applied to the
          data matrix prior to plotting or when calling the \code{score}
          method. The function should accept exactly one input argument
          and its return value needs to be a numeric vector which can
          be coerced back into a data matrix of identical dimensionality
          as the input data.} 

        \item{}{\code{type="p"}: Character vector. The plot type, one
          or several in \code{c("p","l", "b", "a", "a_confint", "s",
          "g", "r", "S", "confint", "smooth", "histogram", "mountain",
          "polygon", "h", "boxplot", "gradient", "heatmap", "horizon")}.
          See 'Details' section in \code{\linkS4class{DataTrack}} for
          more information on the individual plotting types.} 

        \item{}{\code{varwidth=FALSE}: Logical scalar. Parameter
          controlling the boxplot appearance. See \code{\link{panel.bwplot}}
          for details.} 

        \item{}{\code{window=NULL}: Numeric or character scalar.
          Aggregate the rows values of the data matrix to \code{window}
          equally sized slices on the data range using the method
          defined in \code{aggregation}. If negative, apply a running
          window of size \code{windowSize} using the same aggregation
          method. Alternatively, the special value \code{auto} causes
          the function to determine the optimal window size to avoid
          overplotting, and \code{fixed} uses fixed-size windows of
          size \code{windowSize}.} 

        \item{}{\code{windowSize=NULL}: Numeric scalar. The size of
          the running window when the value of \code{window} is negative.} 

        \item{}{\code{ylim=NULL}: Numeric vector of length 2. The range
          of the y-axis scale.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{col.border.title="white"}: Integer or character
          scalar. The border color for the title panels.} 

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
          The border color for the title panels} 

        \item{}{\code{fill="lightgray"}: Integer or character scalar.
          Default fill color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text, unless a more specific definition exists.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text, unless a more specific definition exists.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text, unless a more specific definition exists.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
          for all plotting elements, unless there is a more specific
          control defined elsewhere.} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

        \item{}{\code{rotation=0}: The rotation angle for all text
          unless a more specific definiton exists.} 

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

    \item{IdeogramTrack}{: 

      \describe{ 

        \item{}{\code{background.title="transparent"}: Character scalar.
          The background color for the title panel. Defaults to omit
          the background.} 

        \item{}{\code{bevel=0.45}: Numeric scalar, between 0 and 1.
          The level of smoothness for the two ends of the ideogram.} 

        \item{}{\code{cex.bands=0.7}: Numeric scalar. The font expansion
          factor for the chromosome band identifier text.} 

        \item{}{\code{cex=0.8}: Numeric scalar. The overall font
          expansion factor for the chromosome name text.} 

        \item{}{\code{col="red"}: Character scalar. The border color
          used for the highlighting of the currently displayed genomic
          region.} 

        \item{}{\code{col.border.title="transparent"}: Integer or
          character scalar. The border color for the title panels.} 

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{fill="#FFE3E6"}: Character scalar. The fill color
          used for the highlighting of the currently displayed genomic
          region.} 

        \item{}{\code{fontface=1}: Character scalar. The font face for
          the chromosome name text.} 

        \item{}{\code{fontfamily="sans"}: Character scalar. The font
          family for the chromosome name text.} 

        \item{}{\code{fontcolor="#808080"}: Character scalar. The font
          color for the chromosome name text.} 

        \item{}{\code{fontsize=10}: Numeric scalar. The font size for
          the chromosome name text.} 

        \item{}{\code{outline=FALSE}: Logical scalar. Add borders to
          the individual chromosome staining bands.} 

        \item{}{\code{showBandId=FALSE}: Logical scalar. Show the
          identifier for the chromosome bands if there is space for it.} 

        \item{}{\code{lty=1}: Character or integer scalar. The line
          type used for the highlighting of the currently displayed
          genomic region.} 

        \item{}{\code{lwd=1}: Numeric scalar. The line width used for
          the highlighting of the currently displayed genomic region.} 

        \item{}{\code{showId=TRUE}: Logical scalar. Indicate the
          chromosome name next to the ideogram.} 

        \item{}{\code{showTitle=FALSE}: Logical scalar. Plot a title
          panel. Defaults to omit the title panel.} 

        \item{}{\code{size=NULL}: Numeric scalar. The relative size of
          the track. Defaults to automatic size setting. Can also be
          overridden in the \code{\link{plotTracks}} function.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

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
          The border color for the title panels} 

        \item{}{\code{collapse=TRUE}: Boolean controlling whether to
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

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text, unless a more specific definition exists.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

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

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

        \item{}{\code{rotation=0}: The rotation angle for all text
          unless a more specific definiton exists.} 

        \item{}{\code{showAxis=TRUE}: Boolean controlling whether to
          plot a y axis (only applies to track types where axes are
          implemented).} 

        \item{}{\code{v=-1}: Integer scalar. Parameter controlling the
          number of vertical grid lines, see \code{\link{panel.grid}}
          for details.} 

      } 

    } 

    \item{AnnotationTrack}{: 

      \describe{ 

        \item{}{\code{arrowHeadWidth=30}: Numeric scalar. The width of
          the arrow head in pixels if \code{shape} is \code{fixedArrow}.} 

        \item{}{\code{arrowHeadMaxWidth=40}: Numeric scalar. The maximum
          width of the arrow head in pixels if \code{shape} is \code{arrow}.} 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{col.line="darkgray"}: Character scalar. The color
          used for connecting lines between grouped items. Defaults to
          a light gray, but if set to \code{NULL} the same color as
          for the first item in the group is used.} 

        \item{}{\code{col="transparent"}: Character or integer scalar.
          The border color for all track items.} 

        \item{}{\code{featureAnnotation=NULL}: Character scalar. Add
          annotation information to the individual track elements. This
          can be a value in \code{id}, \code{group} or \code{feature}.
          Defaults to \code{id}. Only works if \code{showFeatureId} is
          not \code{FALSE}.} 

        \item{}{\code{fill="lightblue"}: Character or integer scalar.
          The fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{fontfamily.group="sans"}: Character scalar. The
          font family for the group-level annotation.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{fontcolor.item="white"}: Character or integer
          scalar. The font color for item identifiers.} 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{fontsize.group=12}: Numeric scalar. The font size
          for the group-level annotation.} 

        \item{}{\code{groupAnnotation=NULL}: Character scalar. Add
          annotation information as group labels. This can be a value
          in \code{id}, \code{group} or \code{feature}. Defaults to
          \code{group}. Only works if \code{showId} is not \code{FALSE}.} 

        \item{}{\code{just.group="left"}: Character scalar. the
          justification of group labels. Either \code{left}, \code{right},
          \code{above} or \code{below}.} 

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

        \item{}{\code{mergeGroups=FALSE}: Logical scalar. Merge fully
          overlapping groups if \code{collapse==TRUE}.} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details. For feathered bars indicating the strandedness
          of grouped items this also controls the height of the arrow
          feathers.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

        \item{}{\code{rotation.group=0}: Numeric scalar. The degree of
          text rotation for group labels.} 

        \item{}{\code{rotation.item=0}: Numeric scalar. The degree of
          text rotation for item identifiers.} 

        \item{}{\code{shape="arrow"}: Character scalar. The shape in
          which to display the track items. Currently only \code{box},
          \code{arrow}, \code{fixedArrow}, \code{ellipse}, and
          \code{smallArrow} are implemented.} 

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

      } 

      \bold{\emph{Inherited from class StackedTrack:}} 

      \describe{ 

        \item{}{\code{stackHeight=0.75}: Numeric between 0 and 1.
          Controls the vertical size and spacing between stacked
          elements. The number defines the proportion of the total
          available space for the stack that is used to draw the glyphs.
          E.g., a value of 0.5 means that half of the available vertical
          drawing space (for each stacking line) is used for the glyphs,
          and thus one quarter of the available space each is used for
          spacing above and below the glyph. Defaults to 0.75.} 

        \item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse
          the y-ordering of stacked items. I.e., features that are
          plotted on the bottom-most stacks will be moved to the top-most
          stack and vice versa.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{col.border.title="white"}: Integer or character
          scalar. The border color for the title panels.} 

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
          The border color for the title panels} 

        \item{}{\code{collapse=TRUE}: Boolean controlling whether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text, unless a more specific definition exists.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text, unless a more specific definition exists.} 

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

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{min.distance=1}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for details.} 

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

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

    \item{GeneRegionTrack}{: 

      \describe{ 

        \item{}{\code{arrowHeadWidth=10}: Numeric scalar. The width of
          the arrow head in pixels if \code{shape} is \code{fixedArrow}.} 

        \item{}{\code{arrowHeadMaxWidth=20}: Numeric scalar. The maximum
          width of the arrow head in pixels if \code{shape} is \code{arrow}.} 

        \item{}{\code{col=NULL}: Character or integer scalar. The border
          color for all track items. Defaults to using the same color
          as in \code{fill}, also taking into account different track
          \code{features}.} 

        \item{}{\code{collapseTranscripts=FALSE}: Logical or character
          scalar. Can be one in \code{gene}, \code{longest}, \code{shortest}
          or \code{meta}. Merge all transcripts of the same gene into
          one single gene model. In the case of \code{gene} (or
          \code{TRUE}), this will only keep the start location of the
          first exon and the end location of the last exon from all
          transcripts of the gene. For \code{shortest} and \code{longest},
          only the longest or shortest transcript model is retained.
          For \code{meta}, a meta-transcript containing the union of
          all exons is formed (essentially identical to the operation
          \code{reduce(geneModel)}).} 

        \item{}{\code{exonAnnotation=NULL}: Character scalar. Add
          annotation information to the individual exon models. This
          can be a value in \code{symbol}, \code{gene}, \code{transcript},
          \code{exon} or \code{feature}. Defaults to \code{exon}. Only
          works if \code{showExonId} is not \code{FALSE}.} 

        \item{}{\code{fill="orange"}: Character or integer scalar. The
          fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{min.distance=0}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for
          details. Note that a value larger than 0 may lead to UTR
          regions being merged to CDS regions, which in most cases is
          not particularly useful.} 

        \item{}{\code{shape=c("smallArrow", "box")}: Character scalar.
          The shape in which to display the track items. Currently only
          \code{box}, \code{arrow}, \code{ellipse}, and \code{smallArrow}
          are implemented.} 

        \item{}{\code{showExonId=NULL}: Logical scalar. Control whether
          to plot the individual exon identifiers.} 

        \item{}{\code{thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5",
          "3UTR", "5UTR", "miRNA", "lincRNA", "three_prime_UTR",
          "five_prime_UTR")}: Character vector. A listing of feature
          types that should be drawn with thin boxes. Typically those
          are non-coding elements.} 

        \item{}{\code{transcriptAnnotation=NULL}: Character scalar.
          Add annotation information as transcript labels. This can be
          a value in \code{symbol}, \code{gene}, \code{transcript},
          \code{exon} or \code{feature}. Defaults to \code{symbol}.
          Only works if \code{showId} is not \code{FALSE}.} 

      } 

      \bold{\emph{Inherited from class AnnotationTrack:}} 

      \describe{ 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{col.line="darkgray"}: Character scalar. The color
          used for connecting lines between grouped items. Defaults to
          a light gray, but if set to \code{NULL} the same color as
          for the first item in the group is used.} 

        \item{}{\code{featureAnnotation=NULL}: Character scalar. Add
          annotation information to the individual track elements. This
          can be a value in \code{id}, \code{group} or \code{feature}.
          Defaults to \code{id}. Only works if \code{showFeatureId} is
          not \code{FALSE}.} 

        \item{}{\code{fontfamily.group="sans"}: Character scalar. The
          font family for the group-level annotation.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{fontcolor.item="white"}: Character or integer
          scalar. The font color for item identifiers.} 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{fontsize.group=12}: Numeric scalar. The font size
          for the group-level annotation.} 

        \item{}{\code{groupAnnotation=NULL}: Character scalar. Add
          annotation information as group labels. This can be a value
          in \code{id}, \code{group} or \code{feature}. Defaults to
          \code{group}. Only works if \code{showId} is not \code{FALSE}.} 

        \item{}{\code{just.group="left"}: Character scalar. the
          justification of group labels. Either \code{left}, \code{right},
          \code{above} or \code{below}.} 

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

        \item{}{\code{mergeGroups=FALSE}: Logical scalar. Merge fully
          overlapping groups if \code{collapse==TRUE}.} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details. For feathered bars indicating the strandedness
          of grouped items this also controls the height of the arrow
          feathers.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

        \item{}{\code{rotation.group=0}: Numeric scalar. The degree of
          text rotation for group labels.} 

        \item{}{\code{rotation.item=0}: Numeric scalar. The degree of
          text rotation for item identifiers.} 

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

      } 

      \bold{\emph{Inherited from class StackedTrack:}} 

      \describe{ 

        \item{}{\code{stackHeight=0.75}: Numeric between 0 and 1.
          Controls the vertical size and spacing between stacked
          elements. The number defines the proportion of the total
          available space for the stack that is used to draw the glyphs.
          E.g., a value of 0.5 means that half of the available vertical
          drawing space (for each stacking line) is used for the glyphs,
          and thus one quarter of the available space each is used for
          spacing above and below the glyph. Defaults to 0.75.} 

        \item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse
          the y-ordering of stacked items. I.e., features that are
          plotted on the bottom-most stacks will be moved to the top-most
          stack and vice versa.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{col.border.title="white"}: Integer or character
          scalar. The border color for the title panels.} 

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
          The border color for the title panels} 

        \item{}{\code{collapse=TRUE}: Boolean controlling whether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text, unless a more specific definition exists.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text, unless a more specific definition exists.} 

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

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

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

    \item{BiomartGeneRegionTrack}{: 

      \describe{ 

        \item{}{\code{C_segment="burlywood4"}: Character or integer
          scalar. Fill color for annotation objects of type 'C_segment'.} 

        \item{}{\code{D_segment="lightblue"}: Character or integer
          scalar. Fill color for annotation objects of type 'C_segment'.} 

        \item{}{\code{J_segment="dodgerblue2"}: Character or integer
          scalar. Fill color for annotation objects of type 'C_segment'.} 

        \item{}{\code{Mt_rRNA="yellow"}: Character or integer scalar.
          Fill color for annotation objects of type 'Mt_rRNA'.} 

        \item{}{\code{Mt_tRNA="darkgoldenrod"}: Character or integer
          scalar. Fill color for annotation objects of type 'Mt_tRNA'.} 

        \item{}{\code{Mt_tRNA_pseudogene="darkgoldenrod1"}: Character
          or integer scalar. Fill color for annotation objects of type
          'Mt_tRNA_pseudogene'.} 

        \item{}{\code{V_segment="aquamarine"}: Character or integer
          scalar. Fill color for annotation objects of type 'V_segment'.} 

        \item{}{\code{miRNA="cornflowerblue"}: Character or integer
          scalar. Fill color for annotation objects of type 'L_segment'.} 

        \item{}{\code{miRNA_pseudogene="cornsilk"}: Character or integer
          scalar. Fill color for annotation objects of type
          'miRNA_pseudogene'.} 

        \item{}{\code{misc_RNA="cornsilk3"}: Character or integer
          scalar. Fill color for annotation objects of type 'misc_RNA'.} 

        \item{}{\code{misc_RNA_pseudogene="cornsilk4"}: Character or
          integer scalar. Fill color for annotation objects of type
          'misc_RNA_pseudogene'.} 

        \item{}{\code{protein_coding="#FFD58A"}: Character or integer
          scalar. Fill color for annotation objects of type 'protein_coding'.} 

        \item{}{\code{pseudogene="brown1"}: Character or integer scalar.
          Fill color for annotation objects of type 'pseudogene'.} 

        \item{}{\code{rRNA="darkolivegreen1"}: Character or integer
          scalar. Fill color for annotation objects of type 'rRNA'.} 

        \item{}{\code{rRNA_pseudogene="darkolivegreen"}: Character or
          integer scalar. Fill color for annotation objects of type
          'rRNA_pseudogene'.} 

        \item{}{\code{retrotransposed="blueviolet"}: Character or
          integer scalar. Fill color for annotation objects of type
          'retrotransposed'.} 

        \item{}{\code{scRNA="gold4"}: Character or integer scalar. Fill
          color for annotation objects of type 'scRNA'.} 

        \item{}{\code{scRNA_pseudogene="darkorange2"}: Character or
          integer scalar. Fill color for annotation objects of type
          'scRNA_pseudogene'.} 

        \item{}{\code{snRNA="coral"}: Character or integer scalar. Fill
          color for annotation objects of type 'snRNA'.} 

        \item{}{\code{snRNA_pseudogene="coral3"}: Character or integer
          scalar. Fill color for annotation objects of type
          'snRNA_pseudogene'.} 

        \item{}{\code{snoRNA="cyan"}: Character or integer scalar. Fill
          color for annotation objects of type 'snoRNA'.} 

        \item{}{\code{snoRNA_pseudogene="cyan2"}: Character or integer
          scalar. Fill color for annotation objects of type
          'snoRNA_pseudogene'.} 

        \item{}{\code{tRNA_pseudogene="antiquewhite3"}: Character or
          integer scalar. Fill color for annotation objects of type
          'tRNA_pseudogene'.} 

        \item{}{\code{utr3="#FFD58A"}: FIXME: PLEASE ADD PARAMETER
          DESCRIPTION.} 

        \item{}{\code{utr5="#FFD58A"}: FIXME: PLEASE ADD PARAMETER
          DESCRIPTION.} 

        \item{}{\code{verbose=FALSE}: Logical scalar. Report data
          loading events from Bioamart or retrieval from cache.} 

      } 

      \bold{\emph{Inherited from class GeneRegionTrack:}} 

      \describe{ 

        \item{}{\code{arrowHeadWidth=10}: Numeric scalar. The width of
          the arrow head in pixels if \code{shape} is \code{fixedArrow}.} 

        \item{}{\code{arrowHeadMaxWidth=20}: Numeric scalar. The maximum
          width of the arrow head in pixels if \code{shape} is \code{arrow}.} 

        \item{}{\code{col=NULL}: Character or integer scalar. The border
          color for all track items. Defaults to using the same color
          as in \code{fill}, also taking into account different track
          \code{features}.} 

        \item{}{\code{collapseTranscripts=FALSE}: Logical or character
          scalar. Can be one in \code{gene}, \code{longest}, \code{shortest}
          or \code{meta}. Merge all transcripts of the same gene into
          one single gene model. In the case of \code{gene} (or
          \code{TRUE}), this will only keep the start location of the
          first exon and the end location of the last exon from all
          transcripts of the gene. For \code{shortest} and \code{longest},
          only the longest or shortest transcript model is retained.
          For \code{meta}, a meta-transcript containing the union of
          all exons is formed (essentially identical to the operation
          \code{reduce(geneModel)}).} 

        \item{}{\code{exonAnnotation=NULL}: Character scalar. Add
          annotation information to the individual exon models. This
          can be a value in \code{symbol}, \code{gene}, \code{transcript},
          \code{exon} or \code{feature}. Defaults to \code{exon}. Only
          works if \code{showExonId} is not \code{FALSE}.} 

        \item{}{\code{fill="orange"}: Character or integer scalar. The
          fill color for untyped items. This is also used to connect
          grouped items. See \code{\link{grouping}} for details.} 

        \item{}{\code{min.distance=0}: Numeric scalar. The minimum
          pixel distance before collapsing range items, only if
          \code{collapse==TRUE}. See \code{\link{collapsing}} for
          details. Note that a value larger than 0 may lead to UTR
          regions being merged to CDS regions, which in most cases is
          not particularly useful.} 

        \item{}{\code{shape=c("smallArrow", "box")}: Character scalar.
          The shape in which to display the track items. Currently only
          \code{box}, \code{arrow}, \code{ellipse}, and \code{smallArrow}
          are implemented.} 

        \item{}{\code{showExonId=NULL}: Logical scalar. Control whether
          to plot the individual exon identifiers.} 

        \item{}{\code{thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5",
          "3UTR", "5UTR", "miRNA", "lincRNA", "three_prime_UTR",
          "five_prime_UTR")}: Character vector. A listing of feature
          types that should be drawn with thin boxes. Typically those
          are non-coding elements.} 

        \item{}{\code{transcriptAnnotation=NULL}: Character scalar.
          Add annotation information as transcript labels. This can be
          a value in \code{symbol}, \code{gene}, \code{transcript},
          \code{exon} or \code{feature}. Defaults to \code{symbol}.
          Only works if \code{showId} is not \code{FALSE}.} 

      } 

      \bold{\emph{Inherited from class AnnotationTrack:}} 

      \describe{ 

        \item{}{\code{cex.group=0.6}: Numeric scalar. The font expansion
          factor for the group-level annotation.} 

        \item{}{\code{cex=1}: Numeric scalar. The font expansion factor
          for item identifiers.} 

        \item{}{\code{col.line="darkgray"}: Character scalar. The color
          used for connecting lines between grouped items. Defaults to
          a light gray, but if set to \code{NULL} the same color as
          for the first item in the group is used.} 

        \item{}{\code{featureAnnotation=NULL}: Character scalar. Add
          annotation information to the individual track elements. This
          can be a value in \code{id}, \code{group} or \code{feature}.
          Defaults to \code{id}. Only works if \code{showFeatureId} is
          not \code{FALSE}.} 

        \item{}{\code{fontfamily.group="sans"}: Character scalar. The
          font family for the group-level annotation.} 

        \item{}{\code{fontcolor.group="#808080"}: Character or integer
          scalar. The font color for the group-level annotation.} 

        \item{}{\code{fontcolor.item="white"}: Character or integer
          scalar. The font color for item identifiers.} 

        \item{}{\code{fontface.group=2}: Numeric scalar. The font face
          for the group-level annotation.} 

        \item{}{\code{fontsize.group=12}: Numeric scalar. The font size
          for the group-level annotation.} 

        \item{}{\code{groupAnnotation=NULL}: Character scalar. Add
          annotation information as group labels. This can be a value
          in \code{id}, \code{group} or \code{feature}. Defaults to
          \code{group}. Only works if \code{showId} is not \code{FALSE}.} 

        \item{}{\code{just.group="left"}: Character scalar. the
          justification of group labels. Either \code{left}, \code{right},
          \code{above} or \code{below}.} 

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

        \item{}{\code{mergeGroups=FALSE}: Logical scalar. Merge fully
          overlapping groups if \code{collapse==TRUE}.} 

        \item{}{\code{min.height=3}: Numeric scalar. The minimum range
          height in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details. For feathered bars indicating the strandedness
          of grouped items this also controls the height of the arrow
          feathers.} 

        \item{}{\code{min.width=1}: Numeric scalar. The minimum range
          width in pixels to display. All ranges are expanded to this
          size in order to avoid rendering issues. See \code{\link{collapsing}}
          for details.} 

        \item{}{\code{rotation=0}: Numeric scalar. The degree of text
          rotation for item identifiers.} 

        \item{}{\code{rotation.group=0}: Numeric scalar. The degree of
          text rotation for group labels.} 

        \item{}{\code{rotation.item=0}: Numeric scalar. The degree of
          text rotation for item identifiers.} 

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

      } 

      \bold{\emph{Inherited from class StackedTrack:}} 

      \describe{ 

        \item{}{\code{stackHeight=0.75}: Numeric between 0 and 1.
          Controls the vertical size and spacing between stacked
          elements. The number defines the proportion of the total
          available space for the stack that is used to draw the glyphs.
          E.g., a value of 0.5 means that half of the available vertical
          drawing space (for each stacking line) is used for the glyphs,
          and thus one quarter of the available space each is used for
          spacing above and below the glyph. Defaults to 0.75.} 

        \item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse
          the y-ordering of stacked items. I.e., features that are
          plotted on the bottom-most stacks will be moved to the top-most
          stack and vice versa.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{col.border.title="white"}: Integer or character
          scalar. The border color for the title panels.} 

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
          The border color for the title panels} 

        \item{}{\code{collapse=TRUE}: Boolean controlling whether to
          collapse the content of the track to accomodate the minimum
          current device resolution. See \code{\link{collapsing}} for
          details.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text, unless a more specific definition exists.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text, unless a more specific definition exists.} 

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

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

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
          "histogram", "mountain", "polygon", "h", "boxplot", "gradient",
          "heatmap", "horizon")}. See the 'Details' section in
          \code{\linkS4class{DataTrack}} for more information on the
          individual plotting types.} 

        \item{}{\code{fill="#0080ff"}: the fill color for the coverage
          indicator.} 

        \item{}{\code{size=NULL}: the relative size of the track.
          Defaults to size selection based on the underlying data. Can
          be overridden in the \code{\link{plotTracks}} function.} 

        \item{}{\code{collapse=FALSE}: collapse overlapping ranges and
          aggregate the underlying data.} 

      } 

      \bold{\emph{Inherited from class StackedTrack:}} 

      \describe{ 

        \item{}{\code{stackHeight=0.75}: Numeric between 0 and 1.
          Controls the vertical size and spacing between stacked
          elements. The number defines the proportion of the total
          available space for the stack that is used to draw the glyphs.
          E.g., a value of 0.5 means that half of the available vertical
          drawing space (for each stacking line) is used for the glyphs,
          and thus one quarter of the available space each is used for
          spacing above and below the glyph. Defaults to 0.75.} 

        \item{}{\code{reverseStacking=FALSE}: Logical flag. Reverse
          the y-ordering of stacked items. I.e., features that are
          plotted on the bottom-most stacks will be moved to the top-most
          stack and vice versa.} 

      } 

      \bold{\emph{Inherited from class GdObject:}} 

      \describe{ 

        \item{}{\code{alpha=1}: Numeric scalar. The transparency for
          all track items.} 

        \item{}{\code{alpha.title=NULL}: Numeric scalar. The transparency
          for the title panel.} 

        \item{}{\code{background.panel="transparent"}: Integer or
          character scalar. The background color of the content panel.} 

        \item{}{\code{background.title="lightgray"}: Integer or character
          scalar. The background color for the title panel.} 

        \item{}{\code{cex.axis=NULL}: Numeric scalar. The expansion
          factor for the axis annotation. Defaults to \code{NULL}, in
          which case it is automatically determined based on the
          available space.} 

        \item{}{\code{cex.title=NULL}: Numeric scalar. The expansion
          factor for the title panel. This effects the fontsize of both
          the title and the axis, if any. Defaults to \code{NULL},
          which means that the text size is automatically adjusted to
          the available space.} 

        \item{}{\code{cex=1}: Numeric scalar. The overall font expansion
          factor for all text and glyphs, unless a more specific
          definition exists.} 

        \item{}{\code{col.axis="white"}: Integer or character scalar.
          The font and line color for the y axis, if any.} 

        \item{}{\code{col.border.title="white"}: Integer or character
          scalar. The border color for the title panels.} 

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
          The border color for the title panels} 

        \item{}{\code{col="#0080FF"}: Integer or character scalar.
          Default line color setting for all plotting elements, unless
          there is a more specific control defined elsewhere.} 

        \item{}{\code{fontcolor="black"}: Integer or character scalar.
          The font color for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontface.title=2}: Integer or character scalar.
          The font face for the title panels.} 

        \item{}{\code{fontface=1}: Integer or character scalar. The
          font face for all text, unless a more specific definition exists.} 

        \item{}{\code{fontfamily.title="sans"}: Integer or character
          scalar. The font family for the title panels.} 

        \item{}{\code{fontfamily="sans"}: Integer or character scalar.
          The font family for all text, unless a more specific definition
          exists.} 

        \item{}{\code{fontsize=12}: Numeric scalar. The font size for
          all text, unless a more specific definition exists.} 

        \item{}{\code{frame=FALSE}: Boolean. Draw a frame around the
          track when plotting.} 

        \item{}{\code{grid=FALSE}: Boolean, switching on/off the plotting
          of a grid.} 

        \item{}{\code{h=-1}: Integer scalar. Parameter controlling the
          number of horizontal grid lines, see \code{\link{panel.grid}}
          for details.} 

        \item{}{\code{lineheight=1}: Numeric scalar. The font line
          height for all text, unless a more specific definition exists.} 

        \item{}{\code{lty.grid="solid"}: Integer or character scalar.
          Default line type for grid lines, both when \code{type=="g"}
          in \code{\link{DataTrack}}s and when display parameter
          \code{grid==TRUE}.} 

        \item{}{\code{lty="solid"}: Numeric scalar. Default line type
          setting for all plotting elements, unless there is a more
          specific control defined elsewhere.} 

        \item{}{\code{lwd.border.title=1}: Integer scalar. The border
          width for the title panels.} 

        \item{}{\code{lwd.title=1}: Integer scalar. The border width
          for the title panels} 

        \item{}{\code{lwd.grid=1}: Numeric scalar. Default line width
          for grid lines, both when \code{type=="g"} in \code{\link{DataTrack}}s
          and when display parameter \code{grid==TRUE}.} 

        \item{}{\code{lwd=1}: Numeric scalar. Default line width setting
          for all plotting elements, unless there is a more specific
          control defined elsewhere.} 

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

        \item{}{\code{reverseStrand=FALSE}: Logical scalar. Set up the
          plotting coordinates in 3' -> 5' direction if \code{TRUE}.
          This will effectively mirror the plot on the vertical axis.} 

        \item{}{\code{rotation.title=90}: The rotation angle for the
          text in the title panel. Even though this can be adjusted,
          the automatic resizing of the title panel will currently not
          work, so use at own risk.} 

        \item{}{\code{rotation=0}: The rotation angle for all text
          unless a more specific definiton exists.} 

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

  \code{\link{horizonplot}}

  \code{\link{panel.bwplot}}

  \code{\link{panel.grid}}

  \code{\link{panel.loess}}

  \code{\link{panel.xyplot}}

  \code{\link{plotTracks}}

}














