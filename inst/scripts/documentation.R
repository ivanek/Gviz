## smart character coercion
toChar <- function(x)
{
    if(length(x)>1)
        return(sprintf("c(\"%s\")", paste(x, collapse="\", \"")))
    if(is.null(x)) return("NULL")
    if(is.function(x)) return(substitute(x))
    if(is.numeric(x)) return(x)
    if(is.logical(x)) return(as.character(x))
    else return(sprintf("\"%s\"", as.character(x)))
}

## Find the best line break point in a text
.findBestBreak <- function(x, chars=70)
{
    if(nchar(x)<=chars)
        return(x)
    xs <- strsplit(x, " ")[[1]]
    ind <- which.min(abs(cumsum(nchar(xs)+1)-chars))
    xres <- paste(xs[1:ind], collapse=" ", sep="")
    if(length(xs)>=(ind+1))
        xres <- c(xres, paste(xs[(ind+1):length(xs)], collapse=" ", sep=""))
    return(xres)
}

## Emacs-style code formatting
indent <- function(x, level=0, block=TRUE, chars=70)
{
    space <- "  "
    level <- rep(level, length(x))[1:length(x)]
    indent <- sapply(level, function(y) paste(rep(space, y), collapse=""))
    indent2 <- as.vector(sapply(indent, function(y) paste(y, rep(space, as.integer(!block)), collapse="", sep="")))
    xc <- mapply(function(y, z) paste(y, gsub(" {2,}", " ",  gsub(" *\\}", "}", gsub("\\{ *", "{", gsub("\n+", "",z)))),
                                      collapse="", sep=""), indent, x)
    res <- mapply(function(y,z){
        tmp <- .findBestBreak(y)
        xres <- tmp[1]
        while(length(tmp)==2 & nchar(tmp[2])>68)
        {
            tmp <- .findBestBreak(paste(z, tmp[2], collapse="", sep=""), chars=chars)
            xres <- c(xres, tmp[1])
        }
        if(length(tmp)==2)
            xres <- c(xres, paste(z, tmp[2], collapse="", sep=""))
        return(paste(xres, collapse="\n", sep=""))
    }, xc, indent2)
    return(res)
}


.tag <- function(x){
  tmp <- sub("\\\\", "", attr(x, "Rd_tag"))
  if(!length(tmp)) return(NA) else return(tmp)
}


.tags <- function(x) sapply(x, .tag)


.tagValue <- function(x)
{
  if(is.na(.tag(x))) return(x)
  if(.tag(x) == "TEXT") return(as.character(x))
  if(length(x)==1 && length(.tag(x[[1]]))) return(x[[1]])
  if(length(x)==1 && !length(.tag(x[[1]]))) return(as.character(x))
  if(length(x)>1 && length(.tags(x)[!is.na(.tags(x))])==length(x)) return(x)
  if(length(x)==2 && all(is.na(.tags(x))))
  {
    attr(x[[1]], "Rd_tag") <- "_sectionContent"
    attr(x[[2]], "Rd_tag") <- "_sectionContent"
    return(x)
  }
}


.traverseRd <- function(x, tag, output=NULL)
{
  tag <- tolower(tag)
  if(!is.na(.tag(x)) && tag == tolower(.tag(x)))
  {
      output <- c(output, .tagValue(x))
  }
  for(child in x)
  {
    thisTag <- .tag(child)
    if(!is.na(thisTag))
    {
      if(tag == tolower(thisTag))
      {
	output <- c(output, .tagValue(child))
      }
      output <- .traverseRd(.tagValue(child), tag, output)
    }
  }
  return(output)
}


## create a documentation skeleton for the display parameters
displayParsDoc <- function(class, details)
{
    parents <- names(getClassDef(class)@contains)
    pars <- sapply(c(class, parents), function(x) as.list(getClassDef(x)@prototype@dp), simplify=FALSE)
    text <- c("\\section{Display Parameters}{", if(length(pars[[1]]))
          {
              pars[[1]] <- pars[[1]][order(names(pars[[1]]))]
              det <- details[[class]][names(pars[[1]])]
              if(is.null(det))
              {
                  warning("No details available for class '", class, "'")
                  det <- "FIXME: PLEASE ADD PARAMETER DESCRIPTION." } else {
                      det[is.na(det)] <- "FIXME: PLEASE ADD PARAMETER DESCRIPTION."}
              if(any(is.na(det)))
                  warning("There are details missing for class '", class, "'. Please update the documentation")
              c(indent(paste("The following display parameters are set for objects of class \\code{",
                             class, "} upon instantiation, unless one or more of them have already been set ",
                             "by one of the optional sub-class initializers, which always get precedence over ",
                             "these global defaults. See \\code{\\link{settings}} for details on ",
                             "setting graphical parameters for tracks.\n\n  \\describe{\n", sep=""), level=1),
                indent(sprintf("\\item{}{\\code{%s=%s}: %s}\n\n", names(pars[[1]]), sapply(pars[[1]], toChar), det),
                       level=2, block=FALSE), indent("}", level=1))
          } else indent(paste("No formal display parameters are defined for objects of class \\code{",
                              class, "}.\n", sep=""), level=1))
    pp <- pars[-1]
    pp <- pp[sapply(pp, length)>0]
    done <- names(pars[[1]])
    if(!is.null(pp) && length(pp)>0)
    {
        text <- c(text, indent(c(paste("Additional display parameters are being inherited from the respective parent ",
                                       "classes. Note that not all of them may have an effect on the plotting of  \\code{",
                                       class, "} objects.", sep=""),"\\describe{"), level=1:2))
        for(i in names(pp))
        {
            leftovers <- setdiff(names(pp[[i]]), done)
            done <- union(names(pp[[i]]), done)
            if(length(leftovers))
            {
                pp[[i]] <- pp[[i]][leftovers][order(names(pp[[i]][leftovers]))]
                det <- details[[i]][names(pp[[i]])]
                if(is.null(det))
                {
                    warning("No details available for class '", i, "'")
                    det <- "FIXME: PLEASE ADD PARAMETER DESCRIPTION." } else {
                        det[is.na(det)] <- "FIXME: PLEASE ADD PARAMETER DESCRIPTION."}
                if(any(is.na(det)))
                    warning("There are details missing for class '", i, "'. Please update the documentation")
                text <- c(text, indent(c(sprintf("\\item{}{\\code{\\linkS4class{%s}}:", i), "\\describe{"),
                                       level=2:3),
                          indent(sprintf("\\item{}{\\code{%s=%s}: %s}\n", names(pp[[i]]), sapply(pp[[i]], toChar), det),
                                 level=4, block=FALSE),
                          indent(rep("}", 2), level=3:2))
            }
        }
        text <- c(text, indent("}", 1))
    }
    text <- c(text, indent("}", 0))
    return(paste(text, collapse=" \n\n"))
}

## Parse though an Rd file, find 'section' in there and replace its content by 'content'.
## If there are parse errors in the Rd file either before or after injection the file
## will not be altered.
injectContent <- function(file, content, section)
{
    require(tools)
    tmp <- tryCatch(parse_Rd(file), warning=function(x) stop("Error parsing rd file:\n", x))
    tags <- sapply(tmp, attr, "Rd_tag")
    ind1 <- grep(section, tags, ignore.case=TRUE)
    ind <- if(!length(ind1))
    {
        sind <- grep("\\\\section", tags, ignore.case=TRUE)
        sind[grep(section, sapply(tmp[sind], function(x) as.character(x[[1]])), ignore.case=TRUE)]
    }else ind1
    if(length(ind)>1)
        stop("Section '", section, "'is not unique in file '", file, "'")
    if(!length(ind)) tmp <- c(tmp, list(content)) else tmp[[ind]] <- content
    class(tmp) <- "Rd"
    file.copy(file, file.path(dirname(file), paste("~", basename(file), sep="")), overwrite=TRUE)
    writeLines(paste(as.character(as(tmp, "Rd")), collapse=""), file)
    trash <- tryCatch(parse_Rd(file), warning=function(x) {
        file.copy(file.path(dirname(file), paste("~", basename(file), sep="")), file, overwrite=TRUE)
        unlink(file.path(dirname(file), paste("~", basename(file), sep="")))
        warning("Injected code is syntatically incorrect. File '", file, "' has not been changed.\n",
                "Message:\n", x)})
    return(invisible(content))
}

## Create display parameters desction for the settings man page that list all availabe parameters for all classes
allDisplayParsDoc <- function(details)
{
    text <- indent(c("\\section{Display Parameters}{", "\\describe{"), level=0:1)
    for(cl in c("GenomeAxisTrack", "DataTrack", "IdeogramTrack", 
         "AnnotationTrack", "GeneRegionTrack", "BiomartGeneRegionTrack", "AlignedReadTrack"))
    {
        parents <- names(getClassDef(cl)@contains)
        pars <- sapply(c(cl, parents), function(x) as.list(getClassDef(x)@prototype@dp), simplify=FALSE)
        text <- c(text, indent(c(sprintf("\\item{%s}{:", cl), ifelse(length(pars[[1]]), "\\describe{", "{")), level=2:3))
        done <- NULL
        for(p in names(pars))
        {
            todo <- setdiff(names(pars[[p]]), done)
            done <- union(done, names(pars[[p]]))
            if(length(todo))
            {
                det <- details[[p]][todo]
                if(is.null(det))
                {
                    warning("No details available for class '", p, "'")
                    det <- "FIXME: PLEASE ADD PARAMETER DESCRIPTION." }else {
                        det[is.na(det)] <- "FIXME: PLEASE ADD PARAMETER DESCRIPTION."}
                if(any(is.na(det)))
                    warning("There are details missing for class '", p, "'. Please update the documentation")
                if(p!=cl)
                {
                    text <- c(text, indent(c("}", sprintf("\\bold{\\emph{Inherited from class %s:}}", p), "\\describe{"), level=3))
                }
                text <- c(text, indent(sprintf("\\item{}{\\code{%s=%s}: %s}", todo, sapply(pars[[p]][todo], toChar), det), level=4, block=FALSE))
            }
                
        }
        text <- c(text, indent(rep("}", 2), level=3:2))
    }
    text <- c(text, indent(rep("}", 2), level=1:0))
    return(paste(text, collapse=" \n\n", sep=""))
}


updateRdFile <- function(class, outdir)
{
    file <- file.path(outdir, paste(class, "class.Rd", sep="-"))
    content <- displayParsDoc(class, details=details)
    injectContent(file, content, "Display Parameters")
}

updateSettingsFile <- function(outdir)
{
    file <- file.path(outdir, "settings.Rd")
    content <- allDisplayParsDoc(details=details)
    injectContent(file, content, "Display Parameters")
}

updateLinks <- function(outdir, toUpdate)
{
    if(missing(toUpdate))
        toUpdate <-  dir(outdir, pattern="^[^~].*")
    toUpdate <- file.path(outdir, basename(toUpdate))
    res <- NULL
    for(f in toUpdate)
    {
        tmp <- suppressWarnings(parse_Rd(f))
        classes <- unique(.traverseRd(tmp, "linkS4class"))
        functions <- setdiff(unique(.traverseRd(tmp, "link")), classes)
        if(!length(functions) && !length(classes))
            res <- c(res, injectContent(f, "", "seealso"))
        else {
            text <- c("\\seealso{",
                      if(!is.null(classes)) indent(paste("\\code{\\linkS4class{", sort(classes), "}}", sep=""), 1) else "",
                      if(!is.null(functions)) indent(paste("\\code{\\link{", sort(functions), "}}", sep=""), 1) else "",
                      "}")
            res <- c(res, injectContent(f, paste(text, collapse="\n\n"), "seealso"))
        }
    }
    return(res)
}

details <- list(

                IdeogramTrack=c(fill="Character scalar. The fill color used for the highlighting of the currently displayed genomic region.",
                                col="Character scalar. The border color used for the highlighting of the currently displayed genomic region.",
                                lwd="Numeric scalar. The line width used for the highlighting of the currently displayed genomic region.",
                                lty="Character or integer scalar. The line type used for the highlighting of the currently displayed genomic region.",
                                fontcolor="Character scalar. The font color for the chromosome name text.",
                                fontface="Character scalar. The font face for the chromosome name text.",
                                fontfamily="Character scalar. The font family for the chromosome name text.",
                                cex="Numeric scalar. The overall font expansion factor for the chromosome name text.",
                                size="Numeric scalar. The relative size of the track. Defaults to automatic size setting. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                showId="Logical scalar. Indicate the chromosome name next to the ideogram.",	
                                bevel="Numeric scalar, between 0 and 1. The level of smoothness for the two ends of the ideogram.",
                                showTitle="Logical scalar. Plot a title panel. Defaults to omit the title panel.",
                                background.title="Character scalar. The background color for the title panel. Defaults to omit the background.",
                                fontsize="Numeric scalar. The font size for the chromosome name text."),
                
                
                DataTrack=c(jitter.x="Logical scalar. Toggle on jittering on the x axis in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            jitter.y="Logical scalar. Toggle off jittering on the y axis in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            factor="Numeric scalar. Factor to control amount of jittering in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            amount="Numeric scalar. Amount of jittering in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            span="Numeric scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            degree="Numeric scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            family="Character scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            evaluation="Numeric scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            baseline="Numeric scalar. Y-axis position of an optional baseline. This parameter has a special meaning for mountain-type plots, see the 'Details' section in \\code{\\linkS4class{DataTrack}} for more information.",
                            col.baseline="Character scalar. Color for the optional baseline, defaults to the setting of \\code{col}.",
                            pch="Integer scalar. The type of glyph used for plotting symbols.",
                            lwd.baseline="Numeric scalar. Line width of the optional baseline, defaults to the setting of \\code{lwd}.",
                            lty.baseline="Character or numeric scalar. Line type of the optional baseline, defaults to the setting of \\code{lty}.",
                            col="Character vector. The base colors to use for all plot types. Unless \\code{groups} are specified, only the first color in the vector is usually taken.", 
                            col.mountain="Character scalar. Line color in mountain-type plots, defaults to the setting of \\code{col}.",
                            lwd.mountain="Numeric scalar. Line width in mountain-type plots, defaults to the setting of \\code{lwd}.",
                            lty.mountain="Character or numeric scalar. Line type in mountain-type plots, defaults to the setting of \\code{lty}.",
                            fill.mountain="Character vector of length 2. Fill color in mountain-type plots.",
                            fill.histogram="Character scalar. Fill color in histogram-type plots, defaults to the setting of \\code{fill}.",
                            col.histogram="Character scalar. Line color in histogram-type plots.",
                            stackedBars="Logical scalar. When there are several data groups, draw the histogram-type plots as stacked barplots or grouped side by side.",
                            box.ratio="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            box.width="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            varwidth="Logical scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            notch="Logical scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            notch.frac="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            levels.fos="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            stats="Function. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            coef="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            do.out="Logical scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function. By default the size will be set automatically based on the selected plotting type.",
                            type="Character vector. The plot type, one or several in \\code{c(\"p\",\"l\", \"b\", \"a\", \"s\", \"g\", \"r\", \"S\", \"smooth\", \"histogram\", \"mountain\", \"h\", \"boxplot\", \"gradient\", \"heatmap\")}. See 'Details' section in \\code{\\linkS4class{DataTrack}} for more information on the individual plotting types.",
                            cex="Numeric scalar. The default pixel size for plotting symbols.",
                            ncolor="Integer scalar. The number of colors for the 'gradient' plotting type",
                            gradient="Character vector. The base colors for the 'gradient' plotting type.",
                            collapse="Logical scalar. Collapse overlapping ranges and aggregate the underlying data.",
                            min.distance="Numeric scalar. The mimimum distance in pixel below which to collapse ranges.",
                            window="Numeric or character scalar. Aggregate the rows values of the data matrix to \\code{window} equally sized slices on the data range using the method defined in \\code{aggregation}. If negative, apply a running window of size \\code{windowSize} using the same aggregation method. Alternatively, the special value \\code{auto} causes the function to determine the optimal window size to avoid overplotting.",
                            windowSize="Numeric scalar. The size of the running window when the value of \\code{window} is negative.",
                            separator="Numeric scalar. Number of pixels used to separate individual samples in heatmap-type plots.",
                            transformation="Function. Applied to the data  matrix prior to plotting or when calling the \\code{score} method. The function should accept exactly one input argument and its return value needs to be a numeric vector which can be coerced back into a data matrix of identical dimensionality as the input data.",
                            groups="Vector coercable to a factor. Optional sample grouping. See 'Details' section in \\code{\\linkS4class{DataTrack}} for further information.",
                            aggregation="Function or character scalar. Used to aggregate values in windows or for collapsing overlapping items. The function has to accept a numeric vector as a single input parameter and has to return a numeric scalar with the aggregated value. Alternatively, one of the predefined options \\code{mean}, \\code{median} \\code{sum}, \\code{min},  \\code{max} or \\code{extreme} can be supplied as a character scalar. Defaults to \\code{mean}.",
                            aggregateGroups="Logical scalar. Aggregate the values within a sample group using the aggregation funnction specified in the \\code{aggregate} parameter.",
                            ylim="Numeric vector of length 2. The range of the y-axis scale.",
                            h="Integer scalar. Parameter controlling the number of vertical grid lines, see \\code{\\link{panel.grid}} for details.",
                            v="Integer scalar. Parameter controlling the number of vertical grid lines, see \\code{\\link{panel.grid}} for details.",
                            col="Character or integer scalar. The color used for all line and symbol elements, unless there is a more specific control defined elsewhere.",
                            lwd="Integer scalar. The line width for all line elements, unless there is a more specific control defined elsewhere.",
                            lty="Character or integer scalar. The type for all line elements, unless there is a more specific control defined elsewhere.",
                            fill="Character scalar. The fill color for area elements, unless there is a more specific control defined elsewhere.",
                            alpha="Numeric scalar between 0 and 1. The opacity of the plotting elements, if supported by the device.",
                            lwd.grid="Integer scalar. The line width for grid elements. Defaults to the setting of \\code{lwd}.",
                            col.grid="Integer scalar. The line color for grid elements.",
                            lty.grid="Integer scalar. The line type for grid elements. Defaults to the setting of \\code{lty}.",
                            col.line="Character or integer scalar. The color used for line elements. Defaults to the setting of \\code{col}.",
                            col.symbol="Character or integer scalar. The color used for symbol elements. Defaults to the setting of \\code{col}.",
                            na.rm="Boolean controlling whether to discard all NA values when plotting or to keep empty spaces for NAs",
                            legend="Boolean triggering the addition of a legend to the track to indicate groups. This only has an effect if at least two groups are presen.",
                            cex.legend="Numeric scalar. The size factor for the legend text.",
                            fontsize.legend="Numeric scalar. The pixel size for the legend text.",
                            fontface.legend="Integer or character scalar. The font face for the legend text.",
                            fontfamily.legend="Integer or character scalar. The font family for the legend text.",
                            lineheight.legend="Numeric scalar. The line height for the legend text.",
                            fontcolor.legend="Integer or character scalar. The font color for the legend text."),

                StackedTrack=c(reverseStacking="Logical flag. Reverse the y-ordering of stacked items. I.e., features that are plotted on the bottom-most stacks will be moved to the top-most stack and vice versa."),
                
                GdObject=c(fontsize="Numeric scalar. The font size for all text.",
                           fontface="Integer or character scalar. The font face for all text.",
                           fontcolor="Integer or character scalar. The font color for all text.",
                           fontfamily="Integer or character scalar. The font family for all text.",
                           lineheight="Numeric scalar. The font line height for all text.",
                           cex="Numeric scalar. The overall font expansion factor for all text.",
                           col="Integer or character scalar. Default line color setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           fill="Integer or character scalar. Default fill color setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           lwd="Numeric scalar. Default line width setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           lty="Numeric scalar. Default line type setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           col.line="Integer or character scalar. Default colors for plot lines. Usually the same as the global \\code{col} parameter.",
                           col.symbol="Integer or character scalar. Default colors for plot symbols. Usually the same as the global \\code{col} parameter.",
                           col.grid="Integer or character scalar. Default line color for grid lines, both when \\code{type==\"g\"} in \\code{\\link{DataTrack}}s and when display parameter \\code{grid==TRUE}.",
                           lwd.grid="Numeric scalar. Default line width for grid lines, both when \\code{type==\"g\"} in \\code{\\link{DataTrack}}s and when display parameter \\code{grid==TRUE}.",
                           lty.grid="Integer or character scalar. Default line type for grid lines, both when \\code{type==\"g\"} in \\code{\\link{DataTrack}}s and when display parameter \\code{grid==TRUE}.",
                           v="Integer scalar. Parameter controlling the number of vertical grid lines, see \\code{\\link{panel.grid}} for details.",
                           h="Integer scalar. Parameter controlling the number of horizontal grid lines, see \\code{\\link{panel.grid}} for details.",
                           alpha="Numeric scalar. The transparency for all track items.",
                           background.title="Integer or character scalar. The background color for the title panels.",
                           col.title="Integer or character scalar. The font color for the title panels.",
                           col.frame="Integer or character scalar. The line color used for the panel frame, if \\code{frame==TRUE}",
                           cex.title="Numeric scalar. The expansion factor for the title panel. This effects the fontsize of both the title and the axis, if any. Defaults to \\code{NULL}, which means that the text size is automatically adjusted to the available space.",
                           fontfamily.title="Integer or character scalar. The font family for the title panels.",
                           fontface.title="Integer or character scalar. The font face for the title panels.",
                           col.axis="Integer or character scalar. The font and line color for the y axis, if any.",
                           cex.axis="Numeric scalar. The expansion factor for the axis  annotation. Defaults to \\code{NULL}, in which case it is computed based on the available space.",
                           background.panel="Integer or character scalar. The background color of the content panel.",
                           showTitle="Boolean controlling whether to plot a title panel. Although this can be set individually for each track, in multi-track plots as created by \\code{\\link{plotTracks}} there will still be an empty placeholder in case any of the other tracks include a title. The same holds true for axes. Note that the the title panel background color could be set to transparent in order to completely hide it.",
                           showAxis="Boolean controlling whether to plot a y axis (only applies to track types where axes are implemented).",
                           grid="Boolean, switching on/off the plotting of a grid.",
                           collapse="Boolean controlling wether to collapse the content of the track to accomodate the minimum current device resolution. See \\code{\\link{collapsing}} for details.",
                           min.width="Numeric scalar. The minimum range width in pixels to display. All ranges are expanded to this size in order to avoid rendering issues. See \\code{\\link{collapsing}} for details.",
                           min.height="Numeric scalar. The minimum range height in pixels to display. All ranges are expanded to this size in order to avoid rendering issues.  See \\code{\\link{collapsing}} for details.",
                           min.distance="Numeric scalar. The minimum pixel distance before collapsing range items, only if \\code{collapse==TRUE}. See \\code{\\link{collapsing}} for details.",
                           frame="Boolean. Draw a frame around the track when plotting.",
                           size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function.",
                           "..."="additional display parameters are allowed. Those typically take the value of a valid R color descriptors. The parameter names will later be matched to optional track item types as defined in the 'feature' range attribute, and all tracks of the matched types are colored accordingly. See the documentation of the \\code{\\link{GeneRegionTrack}} and \\code{\\link{AnnotationTrack}} classes as well as \\code{\\link{grouping}} for details."),

                
                GenomeAxisTrack=c(col="Character scalar. The color for the axis lines and tickmarks.",
                                  col.range="Character scalar. The border color for highlighted regions on the axis.",
                                  fill.range="Character scalar. The fill color for highlighted regions on the axis.",
                                  col.id="Character scalar. The text color for the optional range annotation.",
                                  cex.id="Numeric scalar. The text size for the optional range annotation.",
                                  showId="Logical scalar. Show the optional range highlighting annotation.",
                                  fontcolor="Character scalar. The font color for the axis annotation text.",
                                  fontsize="Numeric scalar. Font size for the axis annotation text in points.",
                                  lwd="Numeric scalar. The line width for the axis elementes.",
                                  cex="Numeric scalar. The overall font expansion factor for the axis annotation text.",
                                  size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function. Defaults to the ideal size based on the other track settings.",
                                  showTitle="Logical scalar. Plot a title panel. Defaults to omit the title panel.",
                                  background.title="Character scalar. The background color for the title panel. Defaults to omit the background.",
                                  add53="Logical scalar. Add 5' to 3' direction indicators.",
                                  add35="Logical scalar. Add 3' to 5' direction indicators.",
                                  exponent="Numeric scalar. The exponent for the axis coordinates, e.g., 3 means mb, 6 means gb, etc. The default is to automatically determine the optimal exponent.",
                                  labelPos="Character vector, one in \"alternating\", \"revAlternating\", \"above\" or \"below\". The vertical positioning of the axis labels. If \\code{scale} is not \\code{NULL}, the possible values are \"above\", \"below\" and \"beside\".",
                                  littleTicks="Logical scalar. Add more fine-grained tick marks.",
                                  distFromAxis="Numeric scalar. Control the distance of the axis annotation from the tick marks.",
                                  fontface="Character scalar. The font face for the axis annotation text.",
                                  fontfamily="Character scalar. The font family for the axis annotation text.",
                                  scale="Numeric scalar. If not \\code{NULL} a small scale is drawn instead of the full axis, if the value is between 0 and 1 it is interpreted as a fraction of the current plotting region, otherwise as an absolute length value in genomic coordinates."),
                                    
                AnnotationTrack=c(fill="Character or integer scalar. The fill color for untyped items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  col="Character or integer scalar. The border color for all track items.",
                                  col.line="Character scalar. The color used for connecting lines between grouped items. Defaults to a light gray, but if set to \\code{NULL} the same color as for the first item in the group is used.",
                                  lty="Character or integer scalar. The line type for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lwd="Integer scalar. The line width for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lex="Numeric scalar. The line expansion factor for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  fontface="Integer scalar. The font face for item identifiers.",
                                  fontcolor="Character or integer scalar. The font color for item identifiers.",
                                  fontsize="Numeric scalar. The font size for item identifiers.",
                                  fontfamily="Character scalar. The font family for item identifiers.",
                                  lineheight="Numeric scalar. The font line height for item identifiers.",
                                  cex="Numeric scalar. The font expansion factor for item identifiers.",
                                  rotation="Numeric scalar. The degree of text rotation for item identifiers.",
                                  size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                  showFeatureId="Logical scalar. Control whether to plot the individual track item identifiers.",
                                  showId="Logical scalar. Control whether to annotate individual groups.",
                                  cex.group="Numeric scalar. The font expansion factor for the group-level annotation.",
                                  fontface.group="Numeric scalar. The font face for the group-level annotation.",
                                  fontfamily.group="Character scalar. The font family for the group-level annotation.",
                                  fontsize.group="Numeric scalar. The font size for the group-level annotation.",
                                  fontcolor.group="Character or integer scalar. The font color for the group-level annotation.",
                                  shape="Character scalar. The shape in which to display the track items. Currently only \\code{box}, \\code{arrow}, \\code{ellipse}, and \\code{smallArrow} are implemented.",
                                  showOverplotting="Logical scalar. Use a color gradient to show the amount of overplotting for collapsed items. This implies that \\code{collapse==TRUE}",
                                  min.width="Numeric scalar. The minimum range width in pixels to display. All ranges are expanded to this size in order to avoid rendering issues. See \\code{\\link{collapsing}} for details.",
                                  min.height="Numeric scalar. The minimum range height in pixels to display. All ranges are expanded to this size in order to avoid rendering issues.  See \\code{\\link{collapsing}} for details. For feathered bars indicating the strandedness of grouped items this also controls the height of the arrow feathers.",
                                  alpha="Numeric scalar between 0 and 1. The opacity of the plotting elements, if supported by the device.",
                                  mergeGroups="Logical scalar. Merge fully overlapping groups if \\code{collapse==TRUE}."),
                
                DetailsAnnotationTrack=c(details.size="Numeric scalar. The fraction of vertical space of the track used for the details section.",
                                         details.minWidth="Numeric scalar. The minium width in pixels for a details panel, if less space is available no details are plotted.",
                                         detailsConnector.col="Character or integer scalar. Color of the line connecting the \\code{AnnotstionTrack} item with its details panel.",
                                         detailsConnector.lty="Character or integer scalar. Type of connecting line.", 
                                         detailsConnector.lwd="Integer scalar. Line width of the connector.",
                                         detailsConnector.pch="Integer scalar. Type of the connector's ends.",        
                                         detailsConnector.cex="Numeric scalar. Relative size of the connector's end points.",
                                         detailsBorder.lty="Character or integer scalar. Line type of the border around each details panel.",
                                         detailsBorder.lwd="Integer scalar. Line width of the border.",
                                         detailsBorder.col="Character or integer scalar. Line color of the border.",
                                         detailsBorder.fill="Character or integer scalar. Background color of the border.",
                                         details.ratio="Numeric scalar. By default, the plotting method tries to fill all available space of the details panel tiles. Depending on the dimensions of your plot and the number of tiles this may lead to fairly stretched plots. Restricting the ration of width over height can help to fine tune for somewhat more sane graphics in these cases. Essentially this adds some white space in between individual tiles to force the desired ratio. Together with the \\code{size} and \\code{details.size} arguments, which control the vertical extension of the whole track and of the details section, this allows for some fairly generic resizing of the tiles.",
                                         detailsFunArgs="List.Additional arguments that get passed on the the details plotting function.",
                                         groupDetails="Logial scalar. Plot details for feature groups rather than for individual features."),
                
                GeneRegionTrack=c(min.distance="Numeric scalar. The minimum pixel distance before collapsing range items, only if \\code{collapse==TRUE}. See \\code{\\link{collapsing}} for details. Note that a value larger than 0 may lead to UTR regions being merged to CDS regions, which in most cases is not particularly useful.",
                                  fill="Character or integer scalar. The fill color for untyped items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  col="Character or integer scalar. The border color for all track items. Defaults to using the same color as in \\code{fill}, also taking into account different track \\code{features}.",
                                  lty="Character or integer scalar. The line type for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lwd="Integer scalar. The line width for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lex="Numeric scalar. The line expansion factor for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  fontface="Integer scalar. The font face for item identifiers.",
                                  fontcolor="Character or integer scalar. The font color for item identifiers.",
                                  fontsize="Numeric scalar. The font size for item identifiers.",
                                  fontfamily="Character scalar. The font family for item identifiers.",
                                  lineheight="Numeric scalar. The font line height for item identifiers.",
                                  cex="Numeric scalar. The font expansion factor for item identifiers.",
                                  rotation="Numeric scalar. The degree of text rotation for item identifiers.",
                                  size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                  showExonId="Logical scalar. Control whether to plot the individual exon identifiers.",
                                  showId="Logical scalar. Control whether to annotate individual groups.",
                                  cex.group="Numeric scalar. The font expansion factor for the group-level annotation.",
                                  fontface.group="Numeric scalar. The font face for the group-level annotation.",
                                  fontfamily.group="Character scalar. The font family for the group-level annotation.",
                                  fontsize.group="Numeric scalar. The font size for the group-level annotation.",
                                  fontcolor.group="Character or integer scalar. The font color for the group-level annotation.",
                                  shape="Character scalar. The shape in which to display the track items. Currently only \\code{box}, \\code{arrow}, \\code{ellipse}, and \\code{smallArrow} are implemented.",
                                  showOverplotting="Logical scalar. Use a color gradient to show the amount of overplotting for collapsed items. This implies that \\code{collapse==TRUE}",
                                  min.width="Numeric scalar. The minimum range width in pixels to display. All ranges are expanded to this size in order to avoid rendering issues. See \\code{\\link{collapsing}} for details.",
                                  alpha="Numeric scalar between 0 and 1. The opacity of the plotting elements, if supported by the device.",
                                  geneSymbols="Logical scalar. Use human-readable gene symbols or gene IDs for the transcript annotation.",
                                  collapseTranscripts="Logical scalar. Merge all transcripts of the same gene into one single gene model. Essentially, this will only keep the start location of the first exon and the end location of the last exon from all transcripts of a gene.",
                                  thinBoxFeature="Character vector. A listing of feature types that should be drawn with thin boxes. Typically those are non-coding elements."),
                                  BiomartGeneRegionTrack=c("C_segment"="Character or integer scalar. Fill color for annotation objects of type 'C_segment'.",
                                  "D_segment"="Character or integer scalar. Fill color for annotation objects of type 'C_segment'.",
                                  "J_segment"="Character or integer scalar. Fill color for annotation objects of type 'C_segment'.",
                                  "miRNA"="Character or integer scalar. Fill color for annotation objects of type 'L_segment'.",
                                  "miRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'miRNA_pseudogene'.",
                                  "misc_RNA"="Character or integer scalar. Fill color for annotation objects of type 'misc_RNA'.",
                                  "misc_RNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'misc_RNA_pseudogene'.",
                                  "Mt_rRNA"="Character or integer scalar. Fill color for annotation objects of type 'Mt_rRNA'.",
                                  "Mt_tRNA"="Character or integer scalar. Fill color for annotation objects of type 'Mt_tRNA'.",
                                  "Mt_tRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'Mt_tRNA_pseudogene'.",
                                  "protein_coding"="Character or integer scalar. Fill color for annotation objects of type 'protein_coding'.",
                                  "pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'pseudogene'.",
                                  "retrotransposed"="Character or integer scalar. Fill color for annotation objects of type 'retrotransposed'.",
                                  "rRNA"="Character or integer scalar. Fill color for annotation objects of type 'rRNA'.",
                                  "rRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'rRNA_pseudogene'.",
                                  "scRNA"="Character or integer scalar. Fill color for annotation objects of type 'scRNA'.",
                                  "scRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'scRNA_pseudogene'.",
                                  "snoRNA"="Character or integer scalar. Fill color for annotation objects of type 'snoRNA'.",
                                  "snoRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'snoRNA_pseudogene'.",
                                  "snRNA"="Character or integer scalar. Fill color for annotation objects of type 'snRNA'.",
                                  "snRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'snRNA_pseudogene'.",
                                  "tRNA_pseudogene"="Character or integer scalar. Fill color for annotation objects of type 'tRNA_pseudogene'.",
                                  "V_segment"="Character or integer scalar. Fill color for annotation objects of type 'V_segment'."),

                AlignedReadTrack=c(fill="the fill color for the coverage indicator.",
                                   type="the plot type, one or several in \\code{c(\"p\",\"l\", \"b\", \"a\", \"s\", \"g\", \"r\", \"S\", \"smooth\", \"histogram\", \"mountain\", \"h\", \"boxplot\", \"gradient\", \"heatmap\")}. See
	the 'Details' section in \\code{\\linkS4class{DataTrack}} for more information on the individual plotting types.",
                                   size="the relative size of the track. Defaults to size selection based on the underlying data. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                   detail="the amount of detail to plot the data. Either \\code{coverage} to show the coverage only, or \\code{reads} to show individual reads. For large data sets the latter can be very inefficient. Please note that \\code{reads} is only available when the object has been created with option \\code{coverageOnly=FALSE}.",
                                   collapse="collapse overlapping ranges and aggregate the underlying data."),

                SequenceTrack=c(size="Numeric scalar. The size of the track item. Defaults to auto-detect the size based on the other parameter settings.",
                                fontcolor="Character vector. The colors used for the 5 possible nucleotides (G, A, T, C, N). Defaults to use colors as defined in the \\code{biovizBase} package.",
                                fontsize="Numeric scalar. Controls the size of the sequence and thus also the level of plotable details.",
                                fontface="Numeric scalar. The face of the font.",
                                lwd="Numeric scalar. The width of the line when no indiviual letters can be plotted due to size limitations.",
                                col="Character scalar. The color of the line when no indiviual letters can be plotted due to size limitations.",
                                min.width="Numeric scalar. The minimum width of the colored boxes that are drawn when no indiviual letters can be plotted due to size limitations.",
                                showTitle="Logical scalar. Do not show a title panel by default.",
                                background.title="Character scalar. Make the title panel transparent by default.",
                                noLetters="Logical scalar. Always plot colored boxes (or a line) regardles of the available space",
                                add53="Logical scalar. Add a direction indicator",
                                complement="Logical scalar. Plot the sequence complement."))


updateDocumentation <- function(outdir="~/Rpacks/Gviz/man")
{
    library(IRanges)
    library(rtracklayer)
    library(GenomicRanges)
    library(lattice)
    library(biomaRt)
    library(RColorBrewer)
    library(Biobase)
    library(grid)
    library(AnnotationDbi)
    ##source(file.path(dirname(outdir), "inst/scripts/sourcePackage.R"))
    dps <- sapply(c("GdObject", "GenomeAxisTrack", "RangeTrack", "NumericTrack", "DataTrack", "IdeogramTrack", "StackedTrack",
                    "AnnotationTrack", "GeneRegionTrack", "BiomartGeneRegionTrack", "AlignedReadTrack"), updateRdFile, outdir)
    settings <- updateSettingsFile(outdir)              
    links <- updateLinks(outdir)
}
