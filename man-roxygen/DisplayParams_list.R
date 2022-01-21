#' @details
#' ## Display Parameters
#'
#' The following display parameters are set for objects of class `GdObject`
#' upon instantiation, unless one or more of them have already been set by
#' one of the optional sub-class initializers, which always get precedence
#' over these global defaults. See `settings` for details on setting
#' graphical parameters for tracks.
#'
#' * `alpha=1` Numeric scalar. The transparency for all track items.
#' * `alpha.title=NULL` Numeric scalar. The transparency for the title panel.
#' * `background.legend="transparent"` Integer or character scalar.
#'    The background colour for the legend.
#' * `background.panel="transparent"` Integer or character scalar.
#' The background colour of the content panel.
#' * `background.title="lightgray"` Integer or character scalar.
#' The background colour for the title panel.
#' * `cex=1` Numeric scalar. The overall font expansion factor for all text
#' and glyphs, unless a more specific definition exists.
#' * `cex.axis=NULL` Numeric scalar. The expansion factor for the axis
#' annotation. Defaults to `NULL`, in which case it is automatically determined
#'  based on the available space.
#' * `cex.title=NULL` Numeric scalar. The expansion factor for the title
#' panel. This effects the font size of both the title and the axis, if any.
#'  Defaults to NULL, which means that the text size is automatically adjusted
#'   to the available space.
#' * `col="#0080FF"` Integer or character scalar. Default line colour setting
#' for all plotting elements, unless there is a more specific control defined
#' elsewhere.
#' * `col.axis="white"` Integer or character scalar. The font and line colour
#' for the y axis, if any.
#' * `col.border.title="white"` Integer or character scalar. The border
#' colour for the title panels.
#' * `col.frame="lightgray"` Integer or character scalar. The line colour
#' used for the panel frame, if `frame==TRUE`
#' * `col.grid="#808080"` Integer or character scalar. Default line colour
#' for grid lines, both when `type=="g"` in `DataTrack`s and when display
#' parameter `grid==TRUE`.
#' * `col.line=NULL` Integer or character scalar. Default colours for plot
#'  lines. Usually the same as the global col parameter.
#' * `col.symbol=NULL` Integer or character scalar. Default colours for plot
#' symbols. Usually the same as the global col parameter.
#' * `col.title="white"` (Aliases `fontcolour.title`) Integer or character
#' scalar. The border colour for the title panels
#' * `collapse=TRUE` Boolean controlling whether to collapse the content of
#' the track to accommodate the minimum current device resolution.
#' See collapsing for details.
#' * `fill="lightgray"` Integer or character scalar. Default fill colour
#' setting for all plotting elements, unless there is a more specific control
#' defined elsewhere.
#' * `fontcolour="black"` Integer or character scalar. The font colour for
#' all text, unless a more specific definition exists.
#' * `fontface=1` Integer or character scalar. The font face for all text,
#' unless a more specific definition exists.
#' * `fontface.title=2` Integer or character scalar. The font face for the
#' title panels.
#' * `fontfamily="sans"` Integer or character scalar. The font family for all
#' text, unless a more specific definition exists.
#' * `fontfamily.title="sans"` Integer or character scalar. The font family
#' for the title panels.
#' * `fontsize=12` Numeric scalar. The font size for all text, unless a more
#'  specific definition exists.
#' * `frame=FALSE` Boolean. Draw a frame around the track when plotting.
#' * `grid=FALSE` Boolean, switching on/off the plotting of a grid.
#' * `h=-1` Integer scalar. Parameter controlling the number of horizontal
#' grid lines, see panel.grid for details.
#' * `lineheight=1` Numeric scalar. The font line height for all text, unless
#'  a more specific definition exists.
#' * `lty="solid"` Numeric scalar. Default line type setting for all plotting
#' elements, unless there is a more specific control defined elsewhere.
#' * `lty.grid="solid"` Integer or character scalar. Default line type for
#'  grid lines, both when `type=="g"` in `DataTrack`s and when display parameter
#'  `grid==TRUE`.
#' * `lwd=1` Numeric scalar. Default line width setting for all plotting
#'  elements, unless there is a more specific control defined elsewhere.
#' * `lwd.border.title=1` Integer scalar. The border width for the title
#' panels.
#' * `lwd.grid=1` Numeric scalar. Default line width for grid lines, both
#' when `type=="g"` in `DataTrack`s and when display parameter `grid==TRUE`.
#' * `lwd.title=1` Integer scalar. The border width for the title panels
#' * `min.distance=1` Numeric scalar. The minimum pixel distance before
#' collapsing range items, only if `collapse==TRUE`. See collapsing for details.
#' * `min.height=3` Numeric scalar. The minimum range height in pixels to
#' display. All ranges are expanded to this size in order to avoid rendering
#'  issues. See collapsing for details.
#' * `min.width=1` Numeric scalar. The minimum range width in pixels to
#' display. All ranges are expanded to this size in order to avoid rendering
#' issues. See collapsing for details.
#' * `reverseStrand=FALSE` Logical scalar. Set up the plotting coordinates
#'  in 3' -> 5' direction if `TRUE.` This will effectively mirror the plot
#'   on the vertical axis.
#' * `rotation=0` The rotation angle for all text unless a more specific
#' definition exists.
#' * `rotation.title=90` (Aliases rotation.title) The rotation angle for
#' the text in the title panel. Even though this can be adjusted, the automatic
#'  resizing of the title panel will currently not work, so use at own risk.
#' * `showAxis=TRUE` Boolean controlling whether to plot a y axis (only
#' applies to track types where axes are implemented).
#' * `showTitle=TRUE` Boolean controlling whether to plot a title panel.
#'  Although this can be set individually for each track, in multi-track plots
#'  as created by `plotTracks` there will still be an empty place holder in case
#'  any of the other tracks include a title. The same holds true for axes. Note
#'  that the the title panel background colour could be set to transparent in
#'  order to completely hide it.
#' * `size=1` Numeric scalar. The relative size of the track. Can be
#' overridden in the `plotTracks` function.
#' * `v=-1` Integer scalar. Parameter controlling the number of vertical
#' grid lines, see panel.grid for details.
#' * `...` additional display parameters are allowed. Those typically
#' take the value of a valid R colour descriptors. The parameter names will
#' later be matched to optional track item types as defined in the 'feature'
#' range attribute, and all tracks of the matched types are coloured
#' accordingly. See the documentation of the `GeneRegionTrack` and
#' `AnnotationTrack` classes as well as grouping for details.
