#' Setting display parameters to control the look and feel of the plots
#'
#'
#' The genome track plots in this package are all highly customizable by means
#' of so called 'display parameters'. This page highlights the use of these
#' parameters and list all available settings for the different track classes.
#'
#'
#' All of the package's track objects inherit the `dp` slot from the
#' [`GdObject`][GdObject-class] parent class, which is the main container to
#' store an object's display parameters. Internally, the content of this slot
#' has to be an object of class [`DisplayPars`][DisplayPars-class], but the user
#' is usually not exposed to this low level implementation. Instead, there are
#' two main interaction points, namely the individual object constructor
#' functions and the final [`plotTracks`] function. In both cases,
#' all additional arguments that are not caught by any of the formally defined
#' function parameters are being interpreted as additional display parameters
#' and are automatically added to the aforementioned slot. The main difference
#' here is that display parameters that are passed on to the constructor
#' function are specific for an individual track object, whereas those supplied
#' to the [`plotTracks`] function will be applied to all the objects in the
#' plotting list. Not all display parameters have an effect on the plotting of
#' all track classes, and those will be silently ignored.
#'
#' One can query the available display parameters for a given class as well as
#' their default values by calling the [`availableDisplayPars`] function, or by
#' inspecting the man pages of the individual track classes. The structure of
#' the classes defined in this package is hierarchical, and so are the available
#' display parameters, i.e., all objects inherit the parameters defined in the
#' common [`GdObject`][GdObject-class] parent class, and so on.
#'
#' Once a track object has been created, the display parameters are still open
#' for modification. To this end, the [`DisplayPars`][DisplayPars-class]
#' replacement method is available for all objects inheriting from class
#' [`GdObject`][GdObject-class]. The method takes a named list of parameters as
#' input, e.g.:
#'
#' `displayPars(foo) <- list(col="red", lwd=2)`
#'
#' In the same spirit, the currently set display parameters for the object
#' `foo` can be inferred using the `displayPars` method directly,
#' e.g.:
#'
#' `displayPars(foo)`
#'
#' For track objects inheriting from class
#' [`AnnotationTrack`][AnnotationTrack-class], display parameters that are not
#' formally defined in the class definition or in any of the parent classes are
#' considered to be valid R color identifiers that are used to distinguish
#' between different types of annotation features. For instance, the parameter
#' 'miRNA' will be used to color all annotation features of class miRNA. The
#' annotation types can be set in the constructor function of the track object
#' via the `feature` argument. For most of the tracks that have been
#' inferred from one of the online repositories, this classification will
#' usually be downloaded along with the actual annotation data.
#'
#' Users might find themselves changing the same parameters over and over
#' again, and it would make sense to register these modifications in a central
#' location once and for all. To this end the Gviz package supports display
#' parameter schemes. A scheme is essentially just a bunch of nested named
#' lists, where the names on the first level of nesting should correspond to
#' track class names, and the names on the second level to the display
#' parameters to set. The currently active scheme can be changed by setting
#' the global option `Gviz.scheme`, and a new scheme can be registered by
#' using the `addScheme` function, providing both the list and the name
#' for the new scheme. The `getScheme` function is useful to get the
#' current scheme as a list structure, for instance to use as a skeleton for
#' your own custom scheme.
#'
#' In order to make these settings persistent across R sessions one can create
#' one or several schemes in the global environment in the special object
#' `.GvizSchemes`, for instance by putting the necessary code in the
#' `.Rprofile` file. This object needs to be a named list of schemes, and
#' it will be collected when the Gviz package loads. Its content is then
#' automatically added to the collection of available schemes.
#'
#' Please note that because display parameters are stored with the track
#' objects, a scheme change only has an effect on those objects that are
#' created after the change has taken place.
#'
#' @aliases settings addScheme getScheme
#' @param scheme A named nested list of display parameters, where the first
#' level of nesting represents Gviz track object classes, and the second level
#' of nesting represents parameters.
#' @param name A character scalar with the scheme name.
#' @return
#'
#' The `getScheme` function returns current scheme as a list structure.
#' @section Display Parameters:
#'
#' \describe{
#'
#' \item{GenomeAxisTrack:}{
#'
#' \describe{
#'
#' \item{`add35=FALSE`:}{ Logical scalar. Add 3' to 5' direction
#' indicators.}
#'
#' \item{`add53=FALSE`:}{ Logical scalar. Add 5' to 3' direction
#' indicators.}
#'
#' \item{`background.title="transparent"`:}{ Character scalar.  The
#' background color for the title panel. Defaults to omit the background.}
#'
#' \item{`cex.id=0.7`:}{ Numeric scalar. The text size for the optional
#' range annotation.}
#'
#' \item{`cex=0.8`:}{ Numeric scalar. The overall font expansion factor
#' for the axis annotation text.}
#'
#' \item{`col.border.title="transparent"`:}{ Integer or character scalar.
#' The border color for the title panels.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`col.id="white"`:}{ Character scalar. The text color for the
#' optional range annotation.}
#'
#' \item{`col.range="cornsilk4"`:}{ Character scalar. The border color for
#' highlighted regions on the axis.}
#'
#' \item{`distFromAxis=1`:}{ Numeric scalar. Control the distance of the
#' axis annotation from the tick marks.}
#'
#' \item{`exponent=NULL`:}{ Numeric scalar. The exponent for the axis
#' coordinates, e.g., 3 means mb, 6 means gb, etc. The default is to
#' automatically determine the optimal exponent.}
#'
#' \item{`fill.range="cornsilk3"`:}{ Character scalar. The fill color for
#' highlighted regions on the axis.}
#'
#' \item{`fontcolor="#808080"`:}{ Character scalar. The font color for the
#' axis annotation text.}
#'
#' \item{`fontsize=10`:}{ Numeric scalar. Font size for the axis
#' annotation text in points.}
#'
#' \item{`labelPos="alternating"`:}{ Character vector, one in
#' "alternating", "revAlternating", "above" or "below". The vertical
#' positioning of the axis labels. If `scale` is not `NULL`, the
#' possible values are "above", "below" and "beside".}
#'
#' \item{`littleTicks=FALSE`:}{ Logical scalar. Add more fine-grained tick
#' marks.}
#'
#' \item{`lwd=2`:}{ Numeric scalar. The line width for the axis
#' elements.}
#'
#' \item{`scale=NULL`:}{ Numeric scalar. If not `NULL` a small scale
#' is drawn instead of the full axis, if the value is between 0 and 1 it is
#' interpreted as a fraction of the current plotting region, otherwise as an
#' absolute length value in genomic coordinates.}
#'
#' \item{`showId=FALSE`:}{ Logical scalar. Show the optional range
#' highlighting annotation.}
#'
#' \item{`showTitle=FALSE`:}{ Logical scalar. Plot a title panel. Defaults
#' to omit the title panel.}
#'
#' \item{`ticksAt=NULL`:}{ Numeric scalar. The exact x-position for
#' tickmarks (in base-pairs).}
#'
#' \item{`size=NULL`:}{ Numeric scalar. The relative size of the track.
#' Can be overridden in the [`plotTracks`] function. Defaults to the
#' ideal size based on the other track settings.}
#'
#' \item{`col="darkgray"`:}{ Character scalar. The color for the axis
#' lines and tickmarks.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.line=NULL`:}{ Integer or character scalar.  Default colors
#' for plot lines. Usually the same as the global `col` parameter.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`collapse=TRUE`:}{ Boolean controlling whether to collapse the
#' content of the track to accommodate the minimum current device resolution.
#' See [`collapsing`] for details.}
#'
#' \item{`fill="lightgray"`:}{ Integer or character scalar.  Default fill
#' color setting for all plotting elements, unless there is a more specific
#' control defined elsewhere.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontface=1`:}{ Integer or character scalar. The font face for
#' all text, unless a more specific definition exists.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`fontfamily="sans"`:}{ Integer or character scalar.  The font
#' family for all text, unless a more specific definition exists.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`grid=FALSE`:}{ Boolean, switching on/off the plotting of a
#' grid.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for all
#' text, unless a more specific definition exists.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lty="solid"`:}{ Numeric scalar. Default line type setting for
#' all plotting elements, unless there is a more specific control defined
#' elsewhere.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`min.distance=1`:}{ Numeric scalar. The minimum pixel distance
#' before collapsing range items, only if `collapse==TRUE`. See
#' [`collapsing`] for details.}
#'
#' \item{`min.height=3`:}{ Numeric scalar. The minimum range height in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`rotation=0`:}{ The rotation angle for all text unless a more
#' specific definition exists.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' \item{DataTrack:}{
#'
#' \describe{
#'
#' \item{`aggregateGroups=FALSE`:}{ Logical scalar. Aggregate the values
#' within a sample group using the aggregation function specified in the
#' `aggregation` parameter.}
#'
#' \item{`aggregation="mean"`:}{ Function or character scalar.  Used to
#' aggregate values in windows or for collapsing overlapping items. The
#' function has to accept a numeric vector as a single input parameter and has
#' to return a numeric scalar with the aggregated value. Alternatively, one of
#' the predefined options `mean`, `median` `sum`, `min`,
#' `max` or `extreme` can be supplied as a character scalar. Defaults
#' to `mean`.}
#'
#' \item{`missingAsZero=TRUE`:}{ Logical scalar. Defines how the missing
#' values are treated in the aggregation procedure with running window. Setting
#' it to `TRUE` fills empty positions with zeros, which is default.
#' `FALSE` fills empty positions with `NA`.}
#'
#' \item{`alpha.confint=0.3`:}{ Numeric scalar. The transparency for the
#' confidence intervals in confint-type plots.}
#'
#' \item{`amount=NULL`:}{ Numeric scalar. Amount of jittering in xy-type
#' plots. See [lattice::panel.xyplot] for details.}
#'
#' \item{`baseline=NULL`:}{ Numeric scalar. Y-axis position of an optional
#' baseline. This parameter has a special meaning for mountain-type and
#' polygon-type plots, see the 'Details' section in
#' [`DataTrack`][DataTrack-class] for more information.}
#'
#' \item{`box.legend=FALSE`:}{ Logical scalar. Draw a box around a
#' legend.}
#'
#' \item{`box.ratio=1`:}{ Numeric scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`box.width=NULL`:}{ Numeric scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`grid=FALSE`:}{ Logical vector. Draw a line grid under the track
#' content.}
#'
#' \item{`cex.legend=0.8`:}{ Numeric scalar. The size factor for the
#' legend text.}
#'
#' \item{`cex.sampleNames=NULL`:}{ Numeric scalar. The size factor for the
#' sample names text in heatmap or horizon plots.  Defaults to an automatic
#' setting.}
#'
#' \item{`cex=0.7`:}{ Numeric scalar. The default pixel size for plotting
#' symbols.}
#'
#' \item{`coef=1.5`:}{ Numeric scalar. Parameter controlling the boxplot
#' appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`col.baseline=NULL`:}{ Character scalar. Color for the optional
#' baseline, defaults to the setting of `col`.}
#'
#' \item{`col.confint=NA`:}{ Character vector. Border colors for the
#' confidence intervals for confint-type plots.}
#'
#' \item{`col.boxplotFrame="#808080"`:}{ Character scalar.  Line color of
#' the frame around grouped boxplots.}
#'
#' \item{`col.histogram="#808080"`:}{ Character scalar. Line color in
#' histogram-type plots.}
#'
#' \item{`col.horizon=NA`:}{ The line color for the segments in the
#' `horizon`-type plot. See [latticeExtra::panel.horizonplot] for details.}
#'
#' \item{`col.mountain=NULL`:}{ Character scalar. Line color in
#' mountain-type and polygon-type plots, defaults to the setting of
#' `col`.}
#'
#' \item{`col.sampleNames="white"`:}{ Character or integer scalar. The
#' color used for the sample names in heatmap plots.}
#'
#' \item{`col=c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange",
#' "#00ff00", "brown")`:}{ Character or integer vector.  The color used for all
#' line and symbol elements, unless there is a more specific control defined
#' elsewhere. Unless `groups` are specified, only the first color in the
#' vector is usually regarded.}
#'
#' \item{`collapse=FALSE`:}{ Logical scalar. Collapse overlapping ranges
#' and aggregate the underlying data.}
#'
#' \item{`degree=1`:}{ Numeric scalar. Parameter controlling the loess
#' calculation for smooth and mountain-type plots.  See
#' [lattice::panel.loess] for details.}
#'
#' \item{`do.out=TRUE`:}{ Logical scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`evaluation=50`:}{ Numeric scalar. Parameter controlling the
#' loess calculation for smooth and mountain-type plots. See
#' [lattice::panel.loess] for details.}
#'
#' \item{`factor=0.5`:}{ Numeric scalar. Factor to control amount of
#' jittering in xy-type plots. See [lattice::panel.xyplot] for details.}
#'
#' \item{`family="symmetric"`:}{ Character scalar. Parameter controlling
#' the loess calculation for smooth and mountain-type plots. See
#' [lattice::panel.loess] for details.}
#'
#' \item{`fill.confint=NULL`:}{ Character vector. Fill colors for the
#' confidence intervals for confint-type plots.}
#'
#' \item{`fill.histogram=NULL`:}{ Character scalar. Fill color in
#' histogram-type plots, defaults to the setting of `fill`.}
#'
#' \item{`fill.horizon=c("#B41414", "#E03231", "#F7A99C", "#9FC8DC",
#' "#468CC8", "#0165B3")`:}{ The fill colors for the segments in the
#' `horizon`-type plot. This should be a vector of length six, where the
#' first three entries are the colors for positive changes, and the latter
#' three entries are the colors for negative changes. Defaults to a red-blue
#' color scheme. See [latticeExtra::panel.horizonplot] for details.}
#'
#' \item{`fill.mountain=c("#CCFFFF", "#FFCCFF")`:}{ Character vector of
#' length 2. Fill color in mountain-type and polygon-type plots.}
#'
#' \item{`fontface.legend=NULL`:}{ Integer or character scalar. The font
#' face for the legend text.}
#'
#' \item{`fontfamily.legend=NULL`:}{ Integer or character scalar. The font
#' family for the legend text.}
#'
#' \item{`fontsize.legend=NULL`:}{ Numeric scalar. The pixel size for the
#' legend text.}
#'
#' \item{`fontcolor.legend="#808080"`:}{ Integer or character scalar. The
#' font color for the legend text.}
#'
#' \item{`gradient=c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
#' "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")`:}{ Character vector.
#' The base colors for the `gradient` plotting type or the `heatmap`
#' type with a single group. When plotting heatmaps with more than one group,
#' the `col` parameter can be used to control the group color scheme,
#' however the gradient will always be from white to 'col' and thus does not
#' offer as much flexibility as this `gradient` parameter.}
#'
#' \item{`groups=NULL`:}{ Vector coercible to a factor.  Optional sample
#' grouping. See 'Details' section in [`DataTrack`][DataTrack-class] for
#' further information.}
#'
#' \item{`horizon.origin=0`:}{ The baseline relative to which changes are
#' indicated on the `horizon`-type plot. See [latticeExtra::panel.horizonplot]
#' for details.}
#'
#' \item{`horizon.scale=NULL`:}{ The scale for each of the segments in the
#' `horizon`-type plot. Defaults to 1/3 of the absolute data range. See
#' [latticeExtra::panel.horizonplot] for details.}
#'
#' \item{`jitter.x=FALSE`:}{ Logical scalar. Toggle on jittering on the x
#' axis in xy-type plots. See [lattice::panel.xyplot] for details.}
#'
#' \item{`jitter.y=FALSE`:}{ Logical scalar. Toggle off jittering on the y
#' axis in xy-type plots. See [lattice::panel.xyplot] for details.}
#'
#' \item{`levels.fos=NULL`:}{ Numeric scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`legend=TRUE`:}{ Boolean triggering the addition of a legend to
#' the track to indicate groups. This only has an effect if at least two groups
#' are present.}
#'
#' \item{`lineheight.legend=NULL`:}{ Numeric scalar. The line height for
#' the legend text.}
#'
#' \item{`lty.baseline=NULL`:}{ Character or numeric scalar.  Line type of
#' the optional baseline, defaults to the setting of `lty`.}
#'
#' \item{`lty.mountain=NULL`:}{ Character or numeric scalar.  Line type in
#' mountain-type and polygon-type plots, defaults to the setting of
#' `lty`.}
#'
#' \item{`lwd.baseline=NULL`:}{ Numeric scalar. Line width of the optional
#' baseline, defaults to the setting of `lwd`.}
#'
#' \item{`lwd.mountain=NULL`:}{ Numeric scalar. Line width in
#' mountain-type and polygon-type plots, defaults to the setting of
#' `lwd`.}
#'
#' \item{`min.distance=0`:}{ Numeric scalar. The minimum distance in pixel
#' below which to collapse ranges.}
#'
#' \item{`na.rm=FALSE`:}{ Boolean controlling whether to discard all NA
#' values when plotting or to keep empty spaces for NAs}
#'
#' \item{`ncolor=100`:}{ Integer scalar. The number of colors for the
#' 'gradient' plotting type}
#'
#' \item{`notch.frac=0.5`:}{ Numeric scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`notch=FALSE`:}{ Logical scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`pch=20`:}{ Integer scalar. The type of glyph used for plotting
#' symbols.}
#'
#' \item{`separator=0`:}{ Numeric scalar. Number of pixels used to
#' separate individual samples in heatmap- and horizon-type plots.}
#'
#' \item{`showColorBar=TRUE`:}{ Boolean. Indicate the data range color
#' mapping in the axis for 'heatmap' or 'gradient' types.}
#'
#' \item{`showSampleNames=FALSE`:}{ Boolean. Display the names of the
#' individual samples in a heatmap or a horizon plot.}
#'
#' \item{`size=NULL`:}{ Numeric scalar. The relative size of the track.
#' Can be overridden in the [`plotTracks`] function. By default the
#' size will be set automatically based on the selected plotting type.}
#'
#' \item{`span=0.2`:}{ Numeric scalar. Parameter controlling the loess
#' calculation for smooth and mountain-type plots.  See
#' [lattice::panel.loess] for details.}
#'
#' \item{`stackedBars=TRUE`:}{ Logical scalar. When there are several data
#' groups, draw the histogram-type plots as stacked barplots or grouped side by
#' side.}
#'
#' \item{`stats=X[[i]]`:}{ Function. Parameter controlling the boxplot
#' appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`transformation=NULL`:}{ Function. Applied to the data matrix
#' prior to plotting or when calling the `score` method. The function
#' should accept exactly one input argument and its return value needs to be a
#' numeric vector which can be coerced back into a data matrix of identical
#' dimensionality as the input data.}
#'
#' \item{`type="p"`:}{ Character vector. The plot type, one or several in
#' `p`,`l`, `b`, `a`, `a_confint`, `s`, `g`,
#' `r`, `S`, `confint`, `smooth`, `histogram`,
#' `mountain`, `polygon`, `h`, `boxplot`, `gradient`,
#' `heatmap`, `horizon`. See 'Details' section in
#' [`DataTrack`][DataTrack-class] for more information on the individual
#' plotting types.}
#'
#' \item{`varwidth=FALSE`:}{ Logical scalar. Parameter controlling the
#' boxplot appearance. See [lattice::panel.bwplot] for details.}
#'
#' \item{`window=NULL`:}{ Numeric or character scalar.  Aggregate the rows
#' values of the data matrix to `window` equally sized slices on the data
#' range using the method defined in `aggregation`. If negative, apply a
#' running window of size `windowSize` using the same aggregation method.
#' Alternatively, the special value `auto` causes the function to
#' determine the optimal window size to avoid overplotting, and `fixed`
#' uses fixed-size windows of size `windowSize`.}
#'
#' \item{`windowSize=NULL`:}{ Numeric scalar. The size of the running
#' window when the value of `window` is negative.}
#'
#' \item{`ylim=NULL`:}{ Numeric vector of length 2. The range of the
#' y-axis scale.}
#'
#' \item{`yTicksAt=NULL`:}{ Numeric vector. The points at which y-axis
#' tick-marks are to be drawn. By default, when `NULL`, tickmark locations
#' are computed.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.title="lightgray"`:}{ Integer or character scalar.
#' The background color for the title panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.border.title="white"`:}{ Integer or character scalar. The
#' border color for the title panels.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.line=NULL`:}{ Integer or character scalar.  Default colors
#' for plot lines. Usually the same as the global `col` parameter.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`fill="lightgray"`:}{ Integer or character scalar.  Default fill
#' color setting for all plotting elements, unless there is a more specific
#' control defined elsewhere.}
#'
#' \item{`fontcolor="black"`:}{ Integer or character scalar.  The font
#' color for all text, unless a more specific definition exists.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontface=1`:}{ Integer or character scalar. The font face for
#' all text, unless a more specific definition exists.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`fontfamily="sans"`:}{ Integer or character scalar.  The font
#' family for all text, unless a more specific definition exists.}
#'
#' \item{`fontsize=12`:}{ Numeric scalar. The font size for all text,
#' unless a more specific definition exists.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for all
#' text, unless a more specific definition exists.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lty="solid"`:}{ Numeric scalar. Default line type setting for
#' all plotting elements, unless there is a more specific control defined
#' elsewhere.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`lwd=1`:}{ Numeric scalar. Default line width setting for all
#' plotting elements, unless there is a more specific control defined
#' elsewhere.}
#'
#' \item{`min.height=3`:}{ Numeric scalar. The minimum range height in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`rotation=0`:}{ The rotation angle for all text unless a more
#' specific definition exists.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`showTitle=TRUE`:}{ Boolean controlling whether to plot a title
#' panel. Although this can be set individually for each track, in multi-track
#' plots as created by [`plotTracks`] there will still be an empty
#' placeholder in case any of the other tracks include a title.  The same holds
#' true for axes. Note that the title panel background color could be set
#' to transparent in order to completely hide it.}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' \item{IdeogramTrack:}{
#'
#' \describe{
#'
#' \item{`background.title="transparent"`:}{ Character scalar.  The
#' background color for the title panel. Defaults to omit the background.}
#'
#' \item{`bevel=0.45`:}{ Numeric scalar, between 0 and 1.  The level of
#' smoothness for the two ends of the ideogram.}
#'
#' \item{`centromereShape="triangle"`:}{ Character scalar.  The shape of
#' the centromere. Only "triangle" or "circle" is accepted. Default to
#' "triangle"}
#'
#' \item{`cex.bands=0.7`:}{ Numeric scalar. The font expansion factor for
#' the chromosome band identifier text.}
#'
#' \item{`cex=0.8`:}{ Numeric scalar. The overall font expansion factor
#' for the chromosome name text.}
#'
#' \item{`col="red"`:}{ Character scalar. The border color used for the
#' highlighting of the currently displayed genomic region.}
#'
#' \item{`col.border.title="transparent"`:}{ Integer or character scalar.
#' The border color for the title panels.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`fill="#FFE3E6"`:}{ Character scalar. The fill color used for the
#' highlighting of the currently displayed genomic region.}
#'
#' \item{`fontface=1`:}{ Character scalar. The font face for the
#' chromosome name text.}
#'
#' \item{`fontfamily="sans"`:}{ Character scalar. The font family for the
#' chromosome name text.}
#'
#' \item{`fontcolor="#808080"`:}{ Character scalar. The font color for the
#' chromosome name text.}
#'
#' \item{`fontsize=10`:}{ Numeric scalar. The font size for the chromosome
#' name text.}
#'
#' \item{`outline=FALSE`:}{ Logical scalar. Add borders to the individual
#' chromosome staining bands.}
#'
#' \item{`showBandId=FALSE`:}{ Logical scalar. Show the identifier for the
#' chromosome bands if there is space for it.}
#'
#' \item{`lty=1`:}{ Character or integer scalar. The line type used for
#' the highlighting of the currently displayed genomic region.}
#'
#' \item{`lwd=1`:}{ Numeric scalar. The line width used for the
#' highlighting of the currently displayed genomic region.}
#'
#' \item{`showId=TRUE`:}{ Logical scalar. Indicate the chromosome name
#' next to the ideogram.}
#'
#' \item{`showTitle=FALSE`:}{ Logical scalar. Plot a title panel. Defaults
#' to omit the title panel.}
#'
#' \item{`size=NULL`:}{ Numeric scalar. The relative size of the track.
#' Defaults to automatic size setting. Can also be overridden in the
#' [`plotTracks`] function.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.line=NULL`:}{ Integer or character scalar.  Default colors
#' for plot lines. Usually the same as the global `col` parameter.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`collapse=TRUE`:}{ Boolean controlling whether to collapse the
#' content of the track to accommodate the minimum current device resolution.
#' See [`collapsing`] for details.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`grid=FALSE`:}{ Boolean, switching on/off the plotting of a
#' grid.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for all
#' text, unless a more specific definition exists.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`min.distance=1`:}{ Numeric scalar. The minimum pixel distance
#' before collapsing range items, only if `collapse==TRUE`. See
#' [`collapsing`] for details.}
#'
#' \item{`min.height=3`:}{ Numeric scalar. The minimum range height in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`rotation=0`:}{ The rotation angle for all text unless a more
#' specific definition exists.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' \item{AnnotationTrack:}{
#'
#' \describe{
#'
#' \item{`arrowHeadWidth=30`:}{ Numeric scalar. The width of the arrow
#' head in pixels if `shape` is `fixedArrow`.}
#'
#' \item{`arrowHeadMaxWidth=40`:}{ Numeric scalar. The maximum width of
#' the arrow head in pixels if `shape` is `arrow`.}
#'
#' \item{`cex.group=0.6`:}{ Numeric scalar. The font expansion factor for
#' the group-level annotation.}
#'
#' \item{`cex=1`:}{ Numeric scalar. The font expansion factor for item
#' identifiers.}
#'
#' \item{`col.line="darkgray"`:}{ Character scalar. The color used for
#' connecting lines between grouped items. Defaults to a light gray, but if set
#' to `NULL` the same color as for the first item in the group is used.}
#'
#' \item{`col="transparent"`:}{ Character or integer scalar.  The border
#' color for all track items.}
#'
#' \item{`featureAnnotation=NULL`:}{ Character scalar. Add annotation
#' information to the individual track elements. This can be a value in
#' `id`, `group` or `feature`.  Defaults to `id`. Only
#' works if `showFeatureId` is not `FALSE`.}
#'
#' \item{`fill="lightblue"`:}{ Character or integer scalar.  The fill
#' color for untyped items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`fontfamily.group="sans"`:}{ Character scalar. The font family
#' for the group-level annotation.}
#'
#' \item{`fontcolor.group="#808080"`:}{ Character or integer scalar. The
#' font color for the group-level annotation.}
#'
#' \item{`fontcolor.item="white"`:}{ Character or integer scalar. The font
#' color for item identifiers.}
#'
#' \item{`fontface.group=2`:}{ Numeric scalar. The font face for the
#' group-level annotation.}
#'
#' \item{`fontsize.group=12`:}{ Numeric scalar. The font size for the
#' group-level annotation.}
#'
#' \item{`groupAnnotation=NULL`:}{ Character scalar. Add annotation
#' information as group labels. This can be a value in `id`, `group`
#' or `feature`. Defaults to `group`. Only works if `showId` is
#' not `FALSE`.}
#'
#' \item{`just.group="left"`:}{ Character scalar. the justification of
#' group labels. Either `left`, `right`, `above` or
#' `below`.}
#'
#' \item{`lex=1`:}{ Numeric scalar. The line expansion factor for all
#' track items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for item
#' identifiers.}
#'
#' \item{`lty="solid"`:}{ Character or integer scalar. The line type for
#' all track items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`lwd=1`:}{ Integer scalar. The line width for all track items.
#' This is also used to connect grouped items. See [`grouping`] for
#' details.}
#'
#' \item{`mergeGroups=FALSE`:}{ Logical scalar. Merge fully overlapping
#' groups if `collapse==TRUE`.}
#'
#' \item{`min.height=3`:}{ Numeric scalar. The minimum range height in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details. For feathered
#' bars indicating the strandedness of grouped items this also controls the
#' height of the arrow feathers.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`rotation=0`:}{ Numeric scalar. The degree of text rotation for
#' item identifiers.}
#'
#' \item{`rotation.group=0`:}{ Numeric scalar. The degree of text rotation
#' for group labels.}
#'
#' \item{`rotation.item=0`:}{ Numeric scalar. The degree of text rotation
#' for item identifiers.}
#'
#' \item{`shape="arrow"`:}{ Character scalar. The shape in which to
#' display the track items. Currently only `box`, `arrow`,
#' `fixedArrow`, `ellipse`, and `smallArrow` are implemented.}
#'
#' \item{`showFeatureId=FALSE`:}{ Logical scalar. Control whether to plot
#' the individual track item identifiers.}
#'
#' \item{`showId=FALSE`:}{ Logical scalar. Control whether to annotate
#' individual groups.}
#'
#' \item{`showOverplotting=FALSE`:}{ Logical scalar. Use a color gradient
#' to show the amount of overplotting for collapsed items. This implies that
#' `collapse==TRUE`}
#'
#' \item{`size=1`:}{ Numeric scalar. The relative size of the track. Can
#' be overridden in the [`plotTracks`] function.}
#'
#' }
#'
#' **_Inherited from class [`StackedTrack`][StackedTrack-class]:_**
#'
#' \describe{
#'
#' \item{`stackHeight=0.75`:}{ Numeric between 0 and 1.  Controls the
#' vertical size and spacing between stacked elements. The number defines the
#' proportion of the total available space for the stack that is used to draw
#' the glyphs.  E.g., a value of 0.5 means that half of the available vertical
#' drawing space (for each stacking line) is used for the glyphs, and thus one
#' quarter of the available space each is used for spacing above and below the
#' glyph. Defaults to 0.75.}
#'
#' \item{`reverseStacking=FALSE`:}{ Logical flag. Reverse the y-ordering
#' of stacked items. I.e., features that are plotted on the bottom-most stacks
#' will be moved to the top-most stack and vice versa.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.title="lightgray"`:}{ Integer or character scalar.
#' The background color for the title panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.border.title="white"`:}{ Integer or character scalar. The
#' border color for the title panels.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`collapse=TRUE`:}{ Boolean controlling whether to collapse the
#' content of the track to accommodate the minimum current device resolution.
#' See [`collapsing`] for details.}
#'
#' \item{`fontcolor="black"`:}{ Integer or character scalar.  The font
#' color for all text, unless a more specific definition exists.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontface=1`:}{ Integer or character scalar. The font face for
#' all text, unless a more specific definition exists.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`fontfamily="sans"`:}{ Integer or character scalar.  The font
#' family for all text, unless a more specific definition exists.}
#'
#' \item{`fontsize=12`:}{ Numeric scalar. The font size for all text,
#' unless a more specific definition exists.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`grid=FALSE`:}{ Boolean, switching on/off the plotting of a
#' grid.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`min.distance=1`:}{ Numeric scalar. The minimum pixel distance
#' before collapsing range items, only if `collapse==TRUE`. See
#' [`collapsing`] for details.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`showTitle=TRUE`:}{ Boolean controlling whether to plot a title
#' panel. Although this can be set individually for each track, in multi-track
#' plots as created by [`plotTracks`] there will still be an empty
#' placeholder in case any of the other tracks include a title.  The same holds
#' true for axes. Note that the title panel background color could be set
#' to transparent in order to completely hide it.}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' \item{GeneRegionTrack:}{
#'
#' \describe{
#'
#' \item{`arrowHeadWidth=10`:}{ Numeric scalar. The width of the arrow
#' head in pixels if `shape` is `fixedArrow`.}
#'
#' \item{`arrowHeadMaxWidth=20`:}{ Numeric scalar. The maximum width of
#' the arrow head in pixels if `shape` is `arrow`.}
#'
#' \item{`col=NULL`:}{ Character or integer scalar. The border color for
#' all track items. Defaults to using the same color as in `fill`, also
#' taking into account different track `features`.}
#'
#' \item{`collapseTranscripts=FALSE`:}{ Logical or character scalar. Can
#' be one of `gene`, `longest`, `shortest` or `meta`. Merge
#' all transcripts of the same gene into one single gene model. In the case of
#' `gene` (or `TRUE`), this will only keep the start location of the
#' first exon and the end location of the last exon from all transcripts of the
#' gene. For `shortest` and `longest`, only the longest or shortest
#' transcript model is retained.  For `meta`, a meta-transcript containing
#' the union of all exons is formed (essentially identical to the operation
#' `reduce(geneModel)`).}
#'
#' \item{`exonAnnotation=NULL`:}{ Character scalar. Add annotation
#' information to the individual exon models. This can be a value in
#' `symbol`, `gene`, `transcript`, `exon` or
#' `feature`. Defaults to `exon`. Only works if `showExonId` is
#' not `FALSE`.}
#'
#' \item{`fill="orange"`:}{ Character or integer scalar. The fill color
#' for untyped items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`min.distance=0`:}{ Numeric scalar. The minimum pixel distance
#' before collapsing range items, only if `collapse==TRUE`. See
#' [`collapsing`] for details. Note that a value larger than 0 may
#' lead to UTR regions being merged to CDS regions, which in most cases is not
#' particularly useful.}
#'
#' \item{`shape=c("smallArrow", "box")`:}{ Character scalar.  The shape in
#' which to display the track items. Currently only `box`, `arrow`,
#' `ellipse`, and `smallArrow` are implemented.}
#'
#' \item{`showExonId=NULL`:}{ Logical scalar. Control whether to plot the
#' individual exon identifiers.}
#'
#' \item{`thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "3UTR",
#' "5UTR", "miRNA", "lincRNA", "three_prime_UTR", "five_prime_UTR")`:}{
#' Character vector. A listing of feature types that should be drawn with thin
#' boxes. Typically those are non-coding elements.}
#'
#' \item{`transcriptAnnotation=NULL`:}{ Character scalar.  Add annotation
#' information as transcript labels. This can be a value in `symbol`,
#' `gene`, `transcript`, `exon` or `feature`. Defaults to
#' `symbol`.  Only works if `showId` is not `FALSE`.}
#'
#' }
#'
#' **_Inherited from class [`AnnotationTrack`][AnnotationTrack-class]:_**
#'
#' \describe{
#'
#' \item{`cex.group=0.6`:}{ Numeric scalar. The font expansion factor for
#' the group-level annotation.}
#'
#' \item{`cex=1`:}{ Numeric scalar. The font expansion factor for item
#' identifiers.}
#'
#' \item{`col.line="darkgray"`:}{ Character scalar. The color used for
#' connecting lines between grouped items. Defaults to a light gray, but if set
#' to `NULL` the same color as for the first item in the group is used.}
#'
#' \item{`featureAnnotation=NULL`:}{ Character scalar. Add annotation
#' information to the individual track elements. This can be a value in
#' `id`, `group` or `feature`.  Defaults to `id`. Only
#' works if `showFeatureId` is not `FALSE`.}
#'
#' \item{`fontfamily.group="sans"`:}{ Character scalar. The font family
#' for the group-level annotation.}
#'
#' \item{`fontcolor.group="#808080"`:}{ Character or integer scalar. The
#' font color for the group-level annotation.}
#'
#' \item{`fontcolor.item="white"`:}{ Character or integer scalar. The font
#' color for item identifiers.}
#'
#' \item{`fontface.group=2`:}{ Numeric scalar. The font face for the
#' group-level annotation.}
#'
#' \item{`fontsize.group=12`:}{ Numeric scalar. The font size for the
#' group-level annotation.}
#'
#' \item{`groupAnnotation=NULL`:}{ Character scalar. Add annotation
#' information as group labels. This can be a value in `id`, `group`
#' or `feature`. Defaults to `group`. Only works if `showId` is
#' not `FALSE`.}
#'
#' \item{`just.group="left"`:}{ Character scalar. the justification of
#' group labels. Either `left`, `right`, `above` or
#' `below`.}
#'
#' \item{`lex=1`:}{ Numeric scalar. The line expansion factor for all
#' track items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for item
#' identifiers.}
#'
#' \item{`lty="solid"`:}{ Character or integer scalar. The line type for
#' all track items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`lwd=1`:}{ Integer scalar. The line width for all track items.
#' This is also used to connect grouped items. See [`grouping`] for
#' details.}
#'
#' \item{`mergeGroups=FALSE`:}{ Logical scalar. Merge fully overlapping
#' groups if `collapse==TRUE`.}
#'
#' \item{`min.height=3`:}{ Numeric scalar. The minimum range height in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details. For feathered
#' bars indicating the strandedness of grouped items this also controls the
#' height of the arrow feathers.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`rotation=0`:}{ Numeric scalar. The degree of text rotation for
#' item identifiers.}
#'
#' \item{`rotation.group=0`:}{ Numeric scalar. The degree of text rotation
#' for group labels.}
#'
#' \item{`rotation.item=0`:}{ Numeric scalar. The degree of text rotation
#' for item identifiers.}
#'
#' \item{`showFeatureId=FALSE`:}{ Logical scalar. Control whether to plot
#' the individual track item identifiers.}
#'
#' \item{`showId=FALSE`:}{ Logical scalar. Control whether to annotate
#' individual groups.}
#'
#' \item{`showOverplotting=FALSE`:}{ Logical scalar. Use a color gradient
#' to show the amount of overplotting for collapsed items. This implies that
#' `collapse==TRUE`}
#'
#' \item{`size=1`:}{ Numeric scalar. The relative size of the track. Can
#' be overridden in the [`plotTracks`] function.}
#'
#' }
#'
#' **_Inherited from class [`StackedTrack`][StackedTrack-class]:_**
#'
#' \describe{
#'
#' \item{`stackHeight=0.75`:}{ Numeric between 0 and 1.  Controls the
#' vertical size and spacing between stacked elements. The number defines the
#' proportion of the total available space for the stack that is used to draw
#' the glyphs.  E.g., a value of 0.5 means that half of the available vertical
#' drawing space (for each stacking line) is used for the glyphs, and thus one
#' quarter of the available space each is used for spacing above and below the
#' glyph. Defaults to 0.75.}
#'
#' \item{`reverseStacking=FALSE`:}{ Logical flag. Reverse the y-ordering
#' of stacked items. I.e., features that are plotted on the bottom-most stacks
#' will be moved to the top-most stack and vice versa.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.title="lightgray"`:}{ Integer or character scalar.
#' The background color for the title panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.border.title="white"`:}{ Integer or character scalar. The
#' border color for the title panels.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`collapse=TRUE`:}{ Boolean controlling whether to collapse the
#' content of the track to accommodate the minimum current device resolution.
#' See [`collapsing`] for details.}
#'
#' \item{`fontcolor="black"`:}{ Integer or character scalar.  The font
#' color for all text, unless a more specific definition exists.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontface=1`:}{ Integer or character scalar. The font face for
#' all text, unless a more specific definition exists.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`fontfamily="sans"`:}{ Integer or character scalar.  The font
#' family for all text, unless a more specific definition exists.}
#'
#' \item{`fontsize=12`:}{ Numeric scalar. The font size for all text,
#' unless a more specific definition exists.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`grid=FALSE`:}{ Boolean, switching on/off the plotting of a
#' grid.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`showTitle=TRUE`:}{ Boolean controlling whether to plot a title
#' panel. Although this can be set individually for each track, in multi-track
#' plots as created by [`plotTracks`] there will still be an empty
#' placeholder in case any of the other tracks include a title.  The same holds
#' true for axes. Note that the title panel background color could be set
#' to transparent in order to completely hide it.}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' \item{BiomartGeneRegionTrack:}{
#'
#' \describe{
#'
#' \item{`C_segment="burlywood4"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'C_segment'.}
#'
#' \item{`D_segment="lightblue"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'C_segment'.}
#'
#' \item{`J_segment="dodgerblue2"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'C_segment'.}
#'
#' \item{`Mt_rRNA="yellow"`:}{ Character or integer scalar.  Fill color
#' for annotation objects of type 'Mt_rRNA'.}
#'
#' \item{`Mt_tRNA="darkgoldenrod"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'Mt_tRNA'.}
#'
#' \item{`Mt_tRNA_pseudogene="darkgoldenrod1"`:}{ Character or integer
#' scalar. Fill color for annotation objects of type 'Mt_tRNA_pseudogene'.}
#'
#' \item{`V_segment="aquamarine"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'V_segment'.}
#'
#' \item{`miRNA="cornflowerblue"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'L_segment'.}
#'
#' \item{`miRNA_pseudogene="cornsilk"`:}{ Character or integer scalar.
#' Fill color for annotation objects of type 'miRNA_pseudogene'.}
#'
#' \item{`misc_RNA="cornsilk3"`:}{ Character or integer scalar. Fill color
#' for annotation objects of type 'misc_RNA'.}
#'
#' \item{`misc_RNA_pseudogene="cornsilk4"`:}{ Character or integer scalar.
#' Fill color for annotation objects of type 'misc_RNA_pseudogene'.}
#'
#' \item{`protein_coding="#FFD58A"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'protein_coding'.}
#'
#' \item{`pseudogene="brown1"`:}{ Character or integer scalar.  Fill color
#' for annotation objects of type 'pseudogene'.}
#'
#' \item{`rRNA="darkolivegreen1"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'rRNA'.}
#'
#' \item{`rRNA_pseudogene="darkolivegreen"`:}{ Character or integer
#' scalar. Fill color for annotation objects of type 'rRNA_pseudogene'.}
#'
#' \item{`retrotransposed="blueviolet"`:}{ Character or integer scalar.
#' Fill color for annotation objects of type 'retrotransposed'.}
#'
#' \item{`scRNA="gold4"`:}{ Character or integer scalar. Fill color for
#' annotation objects of type 'scRNA'.}
#'
#' \item{`scRNA_pseudogene="darkorange2"`:}{ Character or integer scalar.
#' Fill color for annotation objects of type 'scRNA_pseudogene'.}
#'
#' \item{`snRNA="coral"`:}{ Character or integer scalar. Fill color for
#' annotation objects of type 'snRNA'.}
#'
#' \item{`snRNA_pseudogene="coral3"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'snRNA_pseudogene'.}
#'
#' \item{`snoRNA="cyan"`:}{ Character or integer scalar. Fill color for
#' annotation objects of type 'snoRNA'.}
#'
#' \item{`snoRNA_pseudogene="cyan2"`:}{ Character or integer scalar. Fill
#' color for annotation objects of type 'snoRNA_pseudogene'.}
#'
#' \item{`tRNA_pseudogene="antiquewhite3"`:}{ Character or integer scalar.
#' Fill color for annotation objects of type 'tRNA_pseudogene'.}
#'
#' \item{`utr3="#FFD58A"`:}{ Character or integer scalar.  Fill color for
#' annotation objects of type 'utr3'.}
#'
#' \item{`utr5="#FFD58A"`:}{ Character or integer scalar.  Fill color for
#' annotation objects of type 'utr5'.}
#'
#' \item{`verbose=FALSE`:}{ Logical scalar. Report data loading events
#' from Biomart or retrieval from cache.}
#'
#' }
#'
#' **_Inherited from class [`GeneRegionTrack`][GeneRegionTrack-class]:_**
#'
#' \describe{
#'
#' \item{`arrowHeadWidth=10`:}{ Numeric scalar. The width of the arrow
#' head in pixels if `shape` is `fixedArrow`.}
#'
#' \item{`arrowHeadMaxWidth=20`:}{ Numeric scalar. The maximum width of
#' the arrow head in pixels if `shape` is `arrow`.}
#'
#' \item{`col=NULL`:}{ Character or integer scalar. The border color for
#' all track items. Defaults to using the same color as in `fill`, also
#' taking into account different track `features`.}
#'
#' \item{`collapseTranscripts=FALSE`:}{ Logical or character scalar. Can
#' be one of `gene`, `longest`, `shortest` or `meta`. Merge
#' all transcripts of the same gene into one single gene model. In the case of
#' `gene` (or `TRUE`), this will only keep the start location of the
#' first exon and the end location of the last exon from all transcripts of the
#' gene. For `shortest` and `longest`, only the longest or shortest
#' transcript model is retained.  For `meta`, a meta-transcript containing
#' the union of all exons is formed (essentially identical to the operation
#' `reduce(geneModel)`).}
#'
#' \item{`exonAnnotation=NULL`:}{ Character scalar. Add annotation
#' information to the individual exon models. This can be a value in
#' `symbol`, `gene`, `transcript`, `exon` or
#' `feature`. Defaults to `exon`. Only works if `showExonId` is
#' not `FALSE`.}
#'
#' \item{`fill="orange"`:}{ Character or integer scalar. The fill color
#' for untyped items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`min.distance=0`:}{ Numeric scalar. The minimum pixel distance
#' before collapsing range items, only if `collapse==TRUE`. See
#' [`collapsing`] for details. Note that a value larger than 0 may
#' lead to UTR regions being merged to CDS regions, which in most cases is not
#' particularly useful.}
#'
#' \item{`shape=c("smallArrow", "box")`:}{ Character scalar.  The shape in
#' which to display the track items. Currently only `box`, `arrow`,
#' `ellipse`, and `smallArrow` are implemented.}
#'
#' \item{`showExonId=NULL`:}{ Logical scalar. Control whether to plot the
#' individual exon identifiers.}
#'
#' \item{`thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "3UTR",
#' "5UTR", "miRNA", "lincRNA", "three_prime_UTR", "five_prime_UTR")`:}{
#' Character vector. A listing of feature types that should be drawn with thin
#' boxes. Typically those are non-coding elements.}
#'
#' \item{`transcriptAnnotation=NULL`:}{ Character scalar.  Add annotation
#' information as transcript labels. This can be a value in `symbol`,
#' `gene`, `transcript`, `exon` or `feature`. Defaults to
#' `symbol`.  Only works if `showId` is not `FALSE`.}
#'
#' }
#'
#' **_Inherited from class [`AnnotationTrack`][AnnotationTrack-class]:_**
#'
#' \describe{
#'
#' \item{`cex.group=0.6`:}{ Numeric scalar. The font expansion factor for
#' the group-level annotation.}
#'
#' \item{`cex=1`:}{ Numeric scalar. The font expansion factor for item
#' identifiers.}
#'
#' \item{`col.line="darkgray"`:}{ Character scalar. The color used for
#' connecting lines between grouped items. Defaults to a light gray, but if set
#' to `NULL` the same color as for the first item in the group is used.}
#'
#' \item{`featureAnnotation=NULL`:}{ Character scalar. Add annotation
#' information to the individual track elements. This can be a value in
#' `id`, `group` or `feature`.  Defaults to `id`. Only
#' works if `showFeatureId` is not `FALSE`.}
#'
#' \item{`fontfamily.group="sans"`:}{ Character scalar. The font family
#' for the group-level annotation.}
#'
#' \item{`fontcolor.group="#808080"`:}{ Character or integer scalar. The
#' font color for the group-level annotation.}
#'
#' \item{`fontcolor.item="white"`:}{ Character or integer scalar. The font
#' color for item identifiers.}
#'
#' \item{`fontface.group=2`:}{ Numeric scalar. The font face for the
#' group-level annotation.}
#'
#' \item{`fontsize.group=12`:}{ Numeric scalar. The font size for the
#' group-level annotation.}
#'
#' \item{`groupAnnotation=NULL`:}{ Character scalar. Add annotation
#' information as group labels. This can be a value in `id`, `group`
#' or `feature`. Defaults to `group`. Only works if `showId` is
#' not `FALSE`.}
#'
#' \item{`just.group="left"`:}{ Character scalar. the justification of
#' group labels. Either `left`, `right`, `above` or
#' `below`.}
#'
#' \item{`lex=1`:}{ Numeric scalar. The line expansion factor for all
#' track items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for item
#' identifiers.}
#'
#' \item{`lty="solid"`:}{ Character or integer scalar. The line type for
#' all track items. This is also used to connect grouped items. See
#' [`grouping`] for details.}
#'
#' \item{`lwd=1`:}{ Integer scalar. The line width for all track items.
#' This is also used to connect grouped items. See [`grouping`] for
#' details.}
#'
#' \item{`mergeGroups=FALSE`:}{ Logical scalar. Merge fully overlapping
#' groups if `collapse==TRUE`.}
#'
#' \item{`min.height=3`:}{ Numeric scalar. The minimum range height in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details. For feathered
#' bars indicating the strandedness of grouped items this also controls the
#' height of the arrow feathers.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`rotation=0`:}{ Numeric scalar. The degree of text rotation for
#' item identifiers.}
#'
#' \item{`rotation.group=0`:}{ Numeric scalar. The degree of text rotation
#' for group labels.}
#'
#' \item{`rotation.item=0`:}{ Numeric scalar. The degree of text rotation
#' for item identifiers.}
#'
#' \item{`showFeatureId=FALSE`:}{ Logical scalar. Control whether to plot
#' the individual track item identifiers.}
#'
#' \item{`showId=FALSE`:}{ Logical scalar. Control whether to annotate
#' individual groups.}
#'
#' \item{`showOverplotting=FALSE`:}{ Logical scalar. Use a color gradient
#' to show the amount of overplotting for collapsed items. This implies that
#' `collapse==TRUE`}
#'
#' \item{`size=1`:}{ Numeric scalar. The relative size of the track. Can
#' be overridden in the [`plotTracks`] function.}
#'
#' }
#'
#' **_Inherited from class [`StackedTrack`][StackedTrack-class]:_**
#'
#' \describe{
#'
#' \item{`stackHeight=0.75`:}{ Numeric between 0 and 1.  Controls the
#' vertical size and spacing between stacked elements. The number defines the
#' proportion of the total available space for the stack that is used to draw
#' the glyphs.  E.g., a value of 0.5 means that half of the available vertical
#' drawing space (for each stacking line) is used for the glyphs, and thus one
#' quarter of the available space each is used for spacing above and below the
#' glyph. Defaults to 0.75.}
#'
#' \item{`reverseStacking=FALSE`:}{ Logical flag. Reverse the y-ordering
#' of stacked items. I.e., features that are plotted on the bottom-most stacks
#' will be moved to the top-most stack and vice versa.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.title="lightgray"`:}{ Integer or character scalar.
#' The background color for the title panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.border.title="white"`:}{ Integer or character scalar. The
#' border color for the title panels.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`collapse=TRUE`:}{ Boolean controlling whether to collapse the
#' content of the track to accommodate the minimum current device resolution.
#' See [`collapsing`] for details.}
#'
#' \item{`fontcolor="black"`:}{ Integer or character scalar.  The font
#' color for all text, unless a more specific definition exists.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontface=1`:}{ Integer or character scalar. The font face for
#' all text, unless a more specific definition exists.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`fontfamily="sans"`:}{ Integer or character scalar.  The font
#' family for all text, unless a more specific definition exists.}
#'
#' \item{`fontsize=12`:}{ Numeric scalar. The font size for all text,
#' unless a more specific definition exists.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`grid=FALSE`:}{ Boolean, switching on/off the plotting of a
#' grid.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`showTitle=TRUE`:}{ Boolean controlling whether to plot a title
#' panel. Although this can be set individually for each track, in multi-track
#' plots as created by [`plotTracks`] there will still be an empty
#' placeholder in case any of the other tracks include a title.  The same holds
#' true for axes. Note that the title panel background color could be set
#' to transparent in order to completely hide it.}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' \item{AlignmentsTrack:}{
#'
#' \describe{
#'
#' \item{`alpha.reads=0.5`:}{ Numeric scalar between 0 and 1. The
#' transparency of the individual read icons. Can be used to indicate
#' overlapping regions in read pairs. Only on supported devices.}
#'
#' \item{`alpha.mismatch=1`:}{ Numeric scalar between 0 and 1. The
#' transparency of the mismatch base information.}
#'
#' \item{`cex=0.7`:}{ Numeric Scalar. The global character expansion
#' factor.}
#'
#' \item{`cex.mismatch=NULL`:}{ Numeric Scalar. The character expansion
#' factor for the mismatch base letters.}
#'
#' \item{`col.coverage=NULL`:}{ Integer or character scalar.  The line
#' color for the coverage profile.}
#'
#' \item{`col.gap="#808080"`:}{ Integer or character scalar.  The color of
#' the line that is bridging the gap regions in gapped alignments.}
#'
#' \item{`col.mates="#E0E0E0"`:}{ Integer or character scalar.  The color
#' of the line that is connecting two paired reads.}
#'
#' \item{`col.deletion="#000000"`:}{ Integer or character scalar. The
#' color of the line that is bridging the deleted regions in alignments.}
#'
#' \item{`col.insertion="#984EA3"`:}{ Integer or character scalar. The
#' color of the line that highlighting insertions in alignments.}
#'
#' \item{`col.mismatch="#808080"`:}{ Integer or character scalar. The box
#' color around mismatch bases.}
#'
#' \item{`col.reads=NULL`:}{ Integer or character scalar.  The box color
#' around reads.}
#'
#' \item{`col.sashimi=NULL`:}{ Integer or character scalar.  The line
#' color for sashimi plots.}
#'
#' \item{`col="#808080"`:}{ Integer or character scalar. The default color
#' of all line elements.}
#'
#' \item{`collapse=FALSE`:}{ Logical scalar. Do not perform any collapsing
#' of overlapping elements. Currently not supported.}
#'
#' \item{`coverageHeight=0.1`:}{ Numeric scalar. The height of the
#' coverage region of the track. Can either be a value between 0 and 1 in which
#' case it is taken as a relative height, or a positive value greater 1 in
#' which case it is interpreted as pixels.}
#'
#' \item{`fill.coverage=NULL`:}{ Integer or character scalar.  The fill
#' color for the coverage profile.}
#'
#' \item{`fill.reads=NULL`:}{ Integer or character scalar.  The fill color
#' for the read icons.}
#'
#' \item{`fill="#BABABA"`:}{ Integer or character scalar.  The default
#' fill color of all plot elements.}
#'
#' \item{`fontface.mismatch=2`:}{ Integer scalar. The font face for
#' mismatch bases.}
#'
#' \item{`lty.coverage=NULL`:}{ Integer or character scalar.  The line
#' type of the coverage profile.}
#'
#' \item{`lty.gap=NULL`:}{ Integer or character scalar. The type of the
#' line that is bridging the gap regions in gapped alignments.}
#'
#' \item{`lty.mates=NULL`:}{ Integer or character scalar.  The type of the
#' line that is connecting two paired reads.}
#'
#' \item{`lty.deletion=NULL`:}{ Integer or character scalar.  The type of
#' the line that is bridging the deleted regions in alignments.}
#'
#' \item{`lty.insertion=NULL`:}{ Integer or character scalar.  The type of
#' the line that highlighting insertions in alignments.}
#'
#' \item{`lty.mismatch=NULL`:}{ Integer or character scalar.  The box line
#' type around mismatch bases.}
#'
#' \item{`lty.reads=NULL`:}{ Integer or character scalar.  The box line
#' type around mismatch reads.}
#'
#' \item{`lty=1`:}{ Integer or character scalar. The default type of all
#' line elements.}
#'
#' \item{`lwd.coverage=NULL`:}{ Integer or character scalar.  The line
#' width of the coverage profile.}
#'
#' \item{`lwd.gap=NULL`:}{ Integer scalar. The width of the line that is
#' bridging the gap regions in gapped alignments.}
#'
#' \item{`lwd.mates=NULL`:}{ Integer scalar. The width of the line that is
#' connecting two paired reads.}
#'
#' \item{`lwd.deletion=NULL`:}{ Integer scalar. The width of the line that
#' is bridging the deleted regions in alignments.}
#'
#' \item{`lwd.insertion=NULL`:}{ Integer scalar. The width of the line
#' that highlighting insertions in alignments.}
#'
#' \item{`lwd.mismatch=NULL`:}{ Integer scalar. The box line width around
#' mismatch bases.}
#'
#' \item{`lwd.reads=NULL`:}{ Integer scalar. The box line width around
#' reads.}
#'
#' \item{`lwd.sashimiMax=10`:}{ Integer scalar. The maximal width of the
#' line in sashimi plots.}
#'
#' \item{`lwd=1`:}{ Integer scalar. The default width of all line
#' elements.}
#'
#' \item{`max.height=10`:}{ Integer scalar. The maximum height of an
#' individual read in pixels. Can be used in combination with `min.height`
#' to control the read and stacking appearance.}
#'
#' \item{`min.height=5`:}{ Integer scalar. The minimum height of an
#' individual read in pixels. Can be used in combination with `max.height`
#' to control the read and stacking appearance.}
#'
#' \item{`minCoverageHeight=50`:}{ Integer scalar. The minimum height of
#' the coverage section. Useful in combination with a relative setting of
#' `coverageHeight`.}
#'
#' \item{`minSashimiHeight=50`:}{ Integer scalar. The minimum height of
#' the sashimi section. Useful in combination with a relative setting of
#' `sashimiHeight`.}
#'
#' \item{`noLetters=FALSE`:}{ Logical scalar. Always plot colored boxes
#' for mismatch bases regardless of the available space.}
#'
#' \item{`sashimiFilter=NULL`:}{ GRanges object. Only junctions which
#' overlap equally with `sashimiFilter` GRanges are shown. Default
#' `NULL`, no filtering.}
#'
#' \item{`sashimiFilterTolerance=0`:}{ Integer scalar. Only used in
#' combination with `sashimiFilter`. It allows to include junctions whose
#' starts/ends are within specified distance from `sashimiFilter` GRanges.
#' This is useful for cases where the aligner did not place the junction reads
#' precisely. Default `0L` , no tolerance.}
#'
#' \item{`sashimiHeight=0.1`:}{ Integer scalar. The height of the sashimi
#' part of the track. Can either be a value between 0 and 1 in which case it is
#' taken as a relative height, or a positive value greater 1 in which case it
#' is interpreted as pixels.}
#'
#' \item{`sashimiScore=1`:}{ Integer scalar. The minimum number of reads
#' supporting the junction.}
#'
#' \item{`sashimiStrand="*"`:}{ Integer scalar. Only reads which have the
#' specified strand are considered to count the junctions.}
#'
#' \item{`sashimiTransformation=NULL`:}{ Function. Applied to the junction
#' score vector prior to plotting. The function should accept exactly one input
#' argument and its return value needs to be a numeric vector of identical
#' length as the input data.}
#'
#' \item{`showIndels=FALSE`:}{ Logical scalar. Consider insertions and
#' deletions in coverage and pile-up. Default is `FALSE`. If set to
#' `TRUE` the deletions defined in CIGAR string are not considered in
#' coverage plot. The deletions are displayed as bridging lines in pile-up
#' track. Insertions are shown as vertical bars.}
#'
#' \item{`showMismatches=TRUE`:}{ Logical scalar. Add mismatch
#' information, either as individual base letters or using color coded bars.
#' This implies that the reference sequence has been provided, either to the
#' class constructor or as part of the track list.}
#'
#' \item{`size=NULL`:}{ Numeric scalar. The size of the track.  Defaults
#' to automatic sizing.}
#'
#' \item{`transformation=NULL`:}{ Function. Applied to the coverage vector
#' prior to plotting. The function should accept exactly one input argument and
#' its return value needs to be a numeric Rle of identical length as the input
#' data.}
#'
#' \item{`type=c("coverage", "pileup")`:}{ Character vector.  The type of
#' information to plot. For `coverage` a coverage plot, potentially
#' augmented by base mismatch information, for `sashimi` a sashimi plot,
#' showing the junctions, and for `pileup` the pileups of the individual
#' reads. These three can be combined.}
#'
#' }
#'
#' **_Inherited from class [`StackedTrack`][StackedTrack-class]:_**
#'
#' \describe{
#'
#' \item{`stackHeight=0.75`:}{ Numeric between 0 and 1.  Controls the
#' vertical size and spacing between stacked elements. The number defines the
#' proportion of the total available space for the stack that is used to draw
#' the glyphs.  E.g., a value of 0.5 means that half of the available vertical
#' drawing space (for each stacking line) is used for the glyphs, and thus one
#' quarter of the available space each is used for spacing above and below the
#' glyph. Defaults to 0.75.}
#'
#' \item{`reverseStacking=FALSE`:}{ Logical flag. Reverse the y-ordering
#' of stacked items. I.e., features that are plotted on the bottom-most stacks
#' will be moved to the top-most stack and vice versa.}
#'
#' }
#'
#' **_Inherited from class [`GdObject`][GdObject-class]:_**
#'
#' \describe{
#'
#' \item{`alpha=1`:}{ Numeric scalar. The transparency for all track
#' items.}
#'
#' \item{`alpha.title=NULL`:}{ Numeric scalar. The transparency for the
#' title panel.}
#'
#' \item{`background.panel="transparent"`:}{ Integer or character scalar.
#' The background color of the content panel.}
#'
#' \item{`background.title="lightgray"`:}{ Integer or character scalar.
#' The background color for the title panel.}
#'
#' \item{`background.legend="transparent"`:}{ Integer or character scalar.
#' The background color for the legend.}
#'
#' \item{`cex.axis=NULL`:}{ Numeric scalar. The expansion factor for the
#' axis annotation. Defaults to `NULL`, in which case it is automatically
#' determined based on the available space.}
#'
#' \item{`cex.title=NULL`:}{ Numeric scalar. The expansion factor for the
#' title panel. This effects the fontsize of both the title and the axis, if
#' any. Defaults to `NULL`, which means that the text size is
#' automatically adjusted to the available space.}
#'
#' \item{`col.axis="white"`:}{ Integer or character scalar.  The font and
#' line color for the y axis, if any.}
#'
#' \item{`col.border.title="white"`:}{ Integer or character scalar. The
#' border color for the title panels.}
#'
#' \item{`col.frame="lightgray"`:}{ Integer or character scalar. The line
#' color used for the panel frame, if `frame==TRUE`}
#'
#' \item{`col.grid="#808080"`:}{ Integer or character scalar.  Default
#' line color for grid lines, both when `type=="g"` in
#' [`DataTrack`][DataTrack-class]s and when display parameter `grid==TRUE`.}
#'
#' \item{`col.line=NULL`:}{ Integer or character scalar.  Default colors
#' for plot lines. Usually the same as the global `col` parameter.}
#'
#' \item{`col.symbol=NULL`:}{ Integer or character scalar.  Default colors
#' for plot symbols. Usually the same as the global `col` parameter.}
#'
#' \item{`col.title="white"`:}{ Integer or character scalar.  The border
#' color for the title panels}
#'
#' \item{`fontcolor="black"`:}{ Integer or character scalar.  The font
#' color for all text, unless a more specific definition exists.}
#'
#' \item{`fontface.title=2`:}{ Integer or character scalar.  The font face
#' for the title panels.}
#'
#' \item{`fontface=1`:}{ Integer or character scalar. The font face for
#' all text, unless a more specific definition exists.}
#'
#' \item{`fontfamily.title="sans"`:}{ Integer or character scalar. The
#' font family for the title panels.}
#'
#' \item{`fontfamily="sans"`:}{ Integer or character scalar.  The font
#' family for all text, unless a more specific definition exists.}
#'
#' \item{`fontsize=12`:}{ Numeric scalar. The font size for all text,
#' unless a more specific definition exists.}
#'
#' \item{`frame=FALSE`:}{ Boolean. Draw a frame around the track when
#' plotting.}
#'
#' \item{`grid=FALSE`:}{ Boolean, switching on/off the plotting of a
#' grid.}
#'
#' \item{`h=-1`:}{ Integer scalar. Parameter controlling the number of
#' horizontal grid lines, see [lattice::panel.grid] for details.}
#'
#' \item{`lineheight=1`:}{ Numeric scalar. The font line height for all
#' text, unless a more specific definition exists.}
#'
#' \item{`lty.grid="solid"`:}{ Integer or character scalar.  Default line
#' type for grid lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s
#' and when display parameter `grid==TRUE`.}
#'
#' \item{`lwd.border.title=1`:}{ Integer scalar. The border width for the
#' title panels.}
#'
#' \item{`lwd.title=1`:}{ Integer scalar. The border width for the title
#' panels}
#'
#' \item{`lwd.grid=1`:}{ Numeric scalar. Default line width for grid
#' lines, both when `type=="g"` in [`DataTrack`][DataTrack-class]s and when
#' display parameter `grid==TRUE`.}
#'
#' \item{`min.distance=1`:}{ Numeric scalar. The minimum pixel distance
#' before collapsing range items, only if `collapse==TRUE`. See
#' [`collapsing`] for details.}
#'
#' \item{`min.width=1`:}{ Numeric scalar. The minimum range width in
#' pixels to display. All ranges are expanded to this size in order to avoid
#' rendering issues. See [`collapsing`] for details.}
#'
#' \item{`reverseStrand=FALSE`:}{ Logical scalar. Set up the plotting
#' coordinates in 3' -> 5' direction if `TRUE`.  This will effectively
#' mirror the plot on the vertical axis.}
#'
#' \item{`rotation.title=90`:}{ The rotation angle for the text in the
#' title panel. Even though this can be adjusted, the automatic resizing of the
#' title panel will currently not work, so use at own risk.}
#'
#' \item{`rotation=0`:}{ The rotation angle for all text unless a more
#' specific definition exists.}
#'
#' \item{`showAxis=TRUE`:}{ Boolean controlling whether to plot a y axis
#' (only applies to track types where axes are implemented).}
#'
#' \item{`showTitle=TRUE`:}{ Boolean controlling whether to plot a title
#' panel. Although this can be set individually for each track, in multi-track
#' plots as created by [`plotTracks`] there will still be an empty
#' placeholder in case any of the other tracks include a title.  The same holds
#' true for axes. Note that the title panel background color could be set
#' to transparent in order to completely hide it.}
#'
#' \item{`v=-1`:}{ Integer scalar. Parameter controlling the number of
#' vertical grid lines, see [lattice::panel.grid] for details.}
#'
#' }
#'
#' }
#'
#' }
#'
#' @name settings
#'
#' @author Florian Hahne
#'
#' @seealso
#' [`AnnotationTrack`][AnnotationTrack-class]
#'
#' [`DataTrack`][DataTrack-class]
#'
#' [`DisplayPars`][DisplayPars-class]
#'
#' [`GdObject`][GdObject-class]
#'
#' [`availableDisplayPars`]
#'
#' [`collapsing`]
#'
#' [`grouping`]
#'
#' [latticeExtra::panel.horizonplot]
#'
#' [lattice::panel.bwplot]
#'
#' [lattice::panel.grid]
#'
#' [lattice::panel.loess]
#'
#' [lattice::panel.xyplot]
#'
#' [`plotTracks`]
#'
#' @examples
#'
#' ## Which scheme is used?
#' getOption("Gviz.scheme")
#'
#' ## Change default settings for GeneRegionTrack
#' scheme <- getScheme()
#' scheme$GeneRegionTrack$fill <- "salmon"
#' scheme$GeneRegionTrack$col <- NULL
#' scheme$GeneRegionTrack$transcriptAnnotation <- "transcript"
#'
#' ## replace default scheme with myScheme
#' addScheme(scheme, "myScheme")
#' options(Gviz.scheme = "myScheme")
#' getOption("Gviz.scheme")
#'
#' data(geneModels)
#' grtrack <- GeneRegionTrack(geneModels,
#'     genome = "hg19", chromosome = "chr7", name = "Gene Model"
#' )
#' plotTracks(grtrack)
NULL
