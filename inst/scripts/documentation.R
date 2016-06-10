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

## Create display parameters section for the settings man page that list all availabe parameters for all classes
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

                IdeogramTrack=c(

                                background.title="Character scalar. The background color for the title panel. Defaults to omit the background.",
                                bevel="Numeric scalar, between 0 and 1. The level of smoothness for the two ends of the ideogram.",
                                cex.bands="Numeric scalar. The  font expansion factor for the chromosome band identifier text.",
                                cex="Numeric scalar. The overall font expansion factor for the chromosome name text.",
                                col="Character scalar. The border color used for the highlighting of the currently displayed genomic region.",
                                fill="Character scalar. The fill color used for the highlighting of the currently displayed genomic region.",
                                fontcolor="Character scalar. The font color for the chromosome name text.",
                                fontface="Character scalar. The font face for the chromosome name text.",
                                fontfamily="Character scalar. The font family for the chromosome name text.",
                                fontsize="Numeric scalar. The font size for the chromosome name text.",
                                lty="Character or integer scalar. The line type used for the highlighting of the currently displayed genomic region.",
                                lwd="Numeric scalar. The line width used for the highlighting of the currently displayed genomic region.",
                                outline="Logical scalar. Add borders to the individual chromosome staining bands.",
                                showBandId="Logical scalar. Show the identifier for the chromosome bands if there is space for it.",
                                showId="Logical scalar. Indicate the chromosome name next to the ideogram.",
                                showTitle="Logical scalar. Plot a title panel. Defaults to omit the title panel.",
                                size="Numeric scalar. The relative size of the track. Defaults to automatic size setting. Can also be overridden in the \\code{\\link{plotTracks}} function."

                                ),


                DataTrack=c(

                            aggregateGroups="Logical scalar. Aggregate the values within a sample group using the aggregation funnction specified in the \\code{aggregation} parameter.",
                            aggregation="Function or character scalar. Used to aggregate values in windows or for collapsing overlapping items. The function has to accept a numeric vector as a single input parameter and has to return a numeric scalar with the aggregated value. Alternatively, one of the predefined options \\code{mean}, \\code{median} \\code{sum}, \\code{min}, \\code{max} or \\code{extreme} can be supplied as a character scalar. Defaults to \\code{mean}.",
                            alpha="Numeric scalar between 0 and 1. The opacity of the plotting elements, if supported by the device.",
                            alpha.confint="Numeric scalar. The tranasparency for the confidence intervalls in confint-type plots.",
                            amount="Numeric scalar. Amount of jittering in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            baseline="Numeric scalar. Y-axis position of an optional baseline. This parameter has a special meaning for mountain-type and polygon-type plots, see the 'Details' section in \\code{\\linkS4class{DataTrack}} for more information.",
                            box.legend="Logical scalar. Draw a box around a legend",
                            box.ratio="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            box.width="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            cex.legend="Numeric scalar. The size factor for the legend text.",
                            cex.sampleNames="Numeric scalar. The size factor for the sample names text in heatmap or horizon plots. Defaults to an automatic setting.",
                            cex="Numeric scalar. The default pixel size for plotting symbols.",
                            coef="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            col.baseline="Character scalar. Color for the optional baseline, defaults to the setting of \\code{col}.",
                            col.confint="Character vector. Border colors for the confidence intervals for confint-type plots.",
                            col.grid="Integer scalar. The line color for grid elements.",
                            col.histogram="Character scalar. Line color in histogram-type plots.",
                            col.horizon="The line color for the segments in the \\code{horizon}-type plot. See \\code{\\link{horizonplot}} for details.",
                            col.line="Character or integer scalar. The color used for line elements. Defaults to the setting of \\code{col}.",
                            col.mountain="Character scalar. Line color in mountain-type and polygon-type plots, defaults to the setting of \\code{col}.",
                            col.sampleNames="Character or integer scalar. The color used for the sample names in heatmap plots.",
                            col.symbol="Character or integer scalar. The color used for symbol elements. Defaults to the setting of \\code{col}.",
                            col="Character or integer vector. The color used for all line and symbol elements, unless there is a more specific control defined elsewhere. Unless \\code{groups} are specified, only the first color in the vector is usually regarded.",
                            collapse="Logical scalar. Collapse overlapping ranges and aggregate the underlying data.",
                            degree="Numeric scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            do.out="Logical scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            evaluation="Numeric scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            factor="Numeric scalar. Factor to control amount of jittering in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            family="Character scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            col.confint="Character vector. Fill colors for the confidence intervals for confint-type plots.",
                            fill.histogram="Character scalar. Fill color in histogram-type plots, defaults to the setting of \\code{fill}.",
                            fill.horizon="The fill colors for the segments in the \\code{horizon}-type plot. This should be a vector of length six, where the first three entries are the colors for positive changes, and the latter three entries are the colors for negative changes. Defaults to a red-blue color scheme. See \\code{\\link{horizonplot}} for details.",
                            fill.mountain="Character vector of length 2. Fill color in mountain-type and polygon-type plots.",
                            fill="Character scalar. The fill color for area elements, unless there is a more specific control defined elsewhere.",
                            fontcolor.legend="Integer or character scalar. The font color for the legend text.",
                            fontface.legend="Integer or character scalar. The font face for the legend text.",
                            fontfamily.legend="Integer or character scalar. The font family for the legend text.",
                            fontsize.legend="Numeric scalar. The pixel size for the legend text.",
                            gradient="Character vector. The base colors for the 'gradient' plotting type or the 'heatmap' type with a single group. When plotting heatmaps with more than one group, the 'col' parameter can be used to control the group color scheme, however the gradient will always be from white to 'col' and thus does not offer as much flexibility as this 'gradient' parameter.",
                            grid="Logical vector. Draw a line grid under the track content.",
                            groups="Vector coercable to a factor. Optional sample grouping. See 'Details' section in \\code{\\linkS4class{DataTrack}} for further information.",
                            h="Integer scalar. Parameter controlling the number of vertical grid lines, see \\code{\\link{panel.grid}} for details.",
                            horizon.origin="The baseline relative to which changes are indicated on the \\code{horizon}-type plot. See \\code{\\link{horizonplot}} for details.",
                            horizon.scale="The scale for each of the segments in the \\code{horizon}-type plot. Defaults to 1/3 of the absolute data range. See \\code{\\link{horizonplot}} for details.",
                            jitter.x="Logical scalar. Toggle on jittering on the x axis in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            jitter.y="Logical scalar. Toggle off jittering on the y axis in xy-type plots. See \\code{\\link{panel.xyplot}} for details.",
                            legend="Boolean triggering the addition of a legend to the track to indicate groups. This only has an effect if at least two groups are presen.",
                            levels.fos="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            lineheight.legend="Numeric scalar. The line height for the legend text.",
                            lty.baseline="Character or numeric scalar. Line type of the optional baseline, defaults to the setting of \\code{lty}.",
                            lty.grid="Integer scalar. The line type for grid elements. Defaults to the setting of \\code{lty}.",
                            lty.mountain="Character or numeric scalar. Line type in mountain-type and polygon-type plots, defaults to the setting of \\code{lty}.",
                            lty="Character or integer scalar. The type for all line elements, unless there is a more specific control defined elsewhere.",
                            lwd.baseline="Numeric scalar. Line width of the optional baseline, defaults to the setting of \\code{lwd}.",
                            lwd.grid="Integer scalar. The line width for grid elements. Defaults to the setting of \\code{lwd}.",
                            lwd.mountain="Numeric scalar. Line width in mountain-type and polygon-type plots, defaults to the setting of \\code{lwd}.",
                            lwd="Integer scalar. The line width for all line elements, unless there is a more specific control defined elsewhere.",
                            min.distance="Numeric scalar. The mimimum distance in pixel below which to collapse ranges.",
                            min.width="Numeric scalar. The mimimum horizontal size of a plotted object in pixels.",
                            na.rm="Boolean controlling whether to discard all NA values when plotting or to keep empty spaces for NAs",
                            ncolor="Integer scalar. The number of colors for the 'gradient' plotting type",
                            notch.frac="Numeric scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            notch="Logical scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            pch="Integer scalar. The type of glyph used for plotting symbols.",
                            separator="Numeric scalar. Number of pixels used to separate individual samples in heatmap- and horizon-type plots.",
                            showColorBar="Boolean. Indicate the data range color mapping in the axis for 'heatmap' or 'gradient' types.",
                            showSampleNames="Boolean. Display the names of the individual samples in a heatmap or a horizon plot.",
                            size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function. By default the size will be set automatically based on the selected plotting type.",
                            span="Numeric scalar. Parameter controlling the loess calculation for smooth and mountain-type plots. See \\code{\\link{panel.loess}} for details.",
                            stackedBars="Logical scalar. When there are several data groups, draw the histogram-type plots as stacked barplots or grouped side by side.",
                            stats="Function. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            transformation="Function. Applied to the data  matrix prior to plotting or when calling the \\code{score} method. The function should accept exactly one input argument and its return value needs to be a numeric vector which can be coerced back into a data matrix of identical dimensionality as the input data.",
                            type="Character vector. The plot type, one or several in \\code{c(\"p\",\"l\", \"b\", \"a\", \"a_confint\", \"s\", \"g\", \"r\", \"S\", \"smooth\", \"histogram\", \"mountain\", \"polygon\", \"h\", \"boxplot\", \"gradient\", \"heatmap\", \"horizon\")}. See 'Details' section in \\code{\\linkS4class{DataTrack}} for more information on the individual plotting types.",
                            v="Integer scalar. Parameter controlling the number of vertical grid lines, see \\code{\\link{panel.grid}} for details.",
                            varwidth="Logical scalar. Parameter controlling the boxplot appearance. See \\code{\\link{panel.bwplot}} for details.",
                            window="Numeric or character scalar. Aggregate the rows values of the data matrix to \\code{window} equally sized slices on the data range using the method defined in \\code{aggregation}. If negative, apply a running window of size \\code{windowSize} using the same aggregation method. Alternatively, the special value \\code{auto} causes the function to determine the optimal window size to avoid overplotting.",
                            windowSize="Numeric scalar. The size of the running window when the value of \\code{window} is negative.",
                            ylim="Numeric vector of length 2. The range of the y-axis scale."

                            ),

                StackedTrack=c(

                               stackHeight="Numeric between 0 and 1. Controls the vertical size and spacing between stacked elements. The number defines the proportion of the total available space for the stack that is used to draw the glyphs. E.g., a value of 0.5 means that half of the available vertical drawing space (for each stacking line) is used for the glyphs, and thus one quarter of the available space each is used for spacing above and below the glyph. Defaults to 0.75.",
                               reverseStacking="Logical flag. Reverse the y-ordering of stacked items. I.e., features that are plotted on the bottom-most stacks will be moved to the top-most stack and vice versa."

                               ),

                GdObject=c(

                           alpha="Numeric scalar. The transparency for all track items.",
                           background.panel="Integer or character scalar. The background color of the content panel.",
                           background.title="Integer or character scalar. The background color for the title panel.",
                           cex.axis="Numeric scalar. The expansion factor for the axis  annotation. Defaults to \\code{NULL}, in which case it is automatically determined based on the available space.",
                           cex.title="Numeric scalar. The expansion factor for the title panel. This effects the fontsize of both the title and the axis, if any. Defaults to \\code{NULL}, which means that the text size is automatically adjusted to the available space.",
                           cex="Numeric scalar. The overall font expansion factor for all text and glyphs, unless a more specific definition exists.",
                           col.axis="Integer or character scalar. The font and line color for the y axis, if any.",
                           col.frame="Integer or character scalar. The line color used for the panel frame, if \\code{frame==TRUE}",
                           col.grid="Integer or character scalar. Default line color for grid lines, both when \\code{type==\"g\"} in \\code{\\link{DataTrack}}s and when display parameter \\code{grid==TRUE}.",
                           col.line="Integer or character scalar. Default colors for plot lines. Usually the same as the global \\code{col} parameter.",
                           col.symbol="Integer or character scalar. Default colors for plot symbols. Usually the same as the global \\code{col} parameter.",
                           col.title="Integer or character scalar. The border color for the title panels",
                           col="Integer or character scalar. Default line color setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           collapse="Boolean controlling whether to collapse the content of the track to accomodate the minimum current device resolution. See \\code{\\link{collapsing}} for details.",
                           fill="Integer or character scalar. Default fill color setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           fontcolor.title="Integer or character scalar. The font color for the title panels.",
                           fontcolor="Integer or character scalar. The font color for all text, unless a more specific definition exists.",
                           fontface.title="Integer or character scalar. The font face for the title panels.",
                           fontface="Integer or character scalar. The font face for all text, unless a more specific definition exists.",
                           fontfamily.title="Integer or character scalar. The font family for the title panels.",
                           fontfamily="Integer or character scalar. The font family for all text, unless a more specific definition exists.",
                           fontsize="Numeric scalar. The font size for all text, unless a more specific definition exists.",
                           frame="Boolean. Draw a frame around the track when plotting.",
                           grid="Boolean, switching on/off the plotting of a grid.",
                           h="Integer scalar. Parameter controlling the number of horizontal grid lines, see \\code{\\link{panel.grid}} for details.",
                           lineheight="Numeric scalar. The font line height for all text, unless a more specific definition exists.",
                           lty.grid="Integer or character scalar. Default line type for grid lines, both when \\code{type==\"g\"} in \\code{\\link{DataTrack}}s and when display parameter \\code{grid==TRUE}.",
                           lty="Numeric scalar. Default line type setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           lwd.grid="Numeric scalar. Default line width for grid lines, both when \\code{type==\"g\"} in \\code{\\link{DataTrack}}s and when display parameter \\code{grid==TRUE}.",
                           lwd.title="Integer scalar. The border width for the title panels",
                           lwd="Numeric scalar. Default line width setting for all plotting elements, unless there is a more specific control defined elsewhere.",
                           min.distance="Numeric scalar. The minimum pixel distance before collapsing range items, only if \\code{collapse==TRUE}. See \\code{\\link{collapsing}} for details.",
                           min.height="Numeric scalar. The minimum range height in pixels to display. All ranges are expanded to this size in order to avoid rendering issues.  See \\code{\\link{collapsing}} for details.",
                           min.width="Numeric scalar. The minimum range width in pixels to display. All ranges are expanded to this size in order to avoid rendering issues. See \\code{\\link{collapsing}} for details.",
                           reverseStrand="Logical scalar. Set up the plotting coordinates in 3' -> 5' direction if \\code{TRUE}. This will effectively mirror the plot on the vertical axis.",
                           rotation.title="The rotation angle for the text in the title panel. Even though this can be adjusted, the automatic resizing of the title panel will currently not work, so use at own risk.",
                           rotation="The rotation angle for all text unless a more specific definiton exists",
                           showAxis="Boolean controlling whether to plot a y axis (only applies to track types where axes are implemented).",
                           showTitle="Boolean controlling whether to plot a title panel. Although this can be set individually for each track, in multi-track plots as created by \\code{\\link{plotTracks}} there will still be an empty placeholder in case any of the other tracks include a title. The same holds true for axes. Note that the the title panel background color could be set to transparent in order to completely hide it.",
                           size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function.",
                           v="Integer scalar. Parameter controlling the number of vertical grid lines, see \\code{\\link{panel.grid}} for details.",
                           "..."="additional display parameters are allowed. Those typically take the value of a valid R color descriptors. The parameter names will later be matched to optional track item types as defined in the 'feature' range attribute, and all tracks of the matched types are colored accordingly. See the documentation of the \\code{\\link{GeneRegionTrack}} and \\code{\\link{AnnotationTrack}} classes as well as \\code{\\link{grouping}} for details."

                           ),

                GenomeAxisTrack=c(

                                  add35="Logical scalar. Add 3' to 5' direction indicators.",
                                  add53="Logical scalar. Add 5' to 3' direction indicators.",
                                  background.title="Character scalar. The background color for the title panel. Defaults to omit the background.",
                                  cex.id="Numeric scalar. The text size for the optional range annotation.",
                                  cex="Numeric scalar. The overall font expansion factor for the axis annotation text.",
                                  col.id="Character scalar. The text color for the optional range annotation.",
                                  col.range="Character scalar. The border color for highlighted regions on the axis.",
                                  col="Character scalar. The color for the axis lines and tickmarks.",
                                  distFromAxis="Numeric scalar. Control the distance of the axis annotation from the tick marks.",
                                  exponent="Numeric scalar. The exponent for the axis coordinates, e.g., 3 means mb, 6 means gb, etc. The default is to automatically determine the optimal exponent.",
                                  fill.range="Character scalar. The fill color for highlighted regions on the axis.",
                                  fontcolor="Character scalar. The font color for the axis annotation text.",
                                  fontface="Character scalar. The font face for the axis annotation text.",
                                  fontfamily="Character scalar. The font family for the axis annotation text.",
                                  fontsize="Numeric scalar. Font size for the axis annotation text in points.",
                                  labelPos="Character vector, one in \"alternating\", \"revAlternating\", \"above\" or \"below\". The vertical positioning of the axis labels. If \\code{scale} is not \\code{NULL}, the possible values are \"above\", \"below\" and \"beside\".",
                                  littleTicks="Logical scalar. Add more fine-grained tick marks.",
                                  lwd="Numeric scalar. The line width for the axis elementes.",
                                  scale="Numeric scalar. If not \\code{NULL} a small scale is drawn instead of the full axis, if the value is between 0 and 1 it is interpreted as a fraction of the current plotting region, otherwise as an absolute length value in genomic coordinates.",
                                  showId="Logical scalar. Show the optional range highlighting annotation.",
                                  showTitle="Logical scalar. Plot a title panel. Defaults to omit the title panel.",
                                  size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function. Defaults to the ideal size based on the other track settings."

                                  ),

                AnnotationTrack=c(

                                  alpha="Numeric scalar between 0 and 1. The opacity of the plotting elements, if supported by the device.",
                                  arrowHeadMaxWidth="Numeric scalar. The maximum width of the arrow head in pixels if \\code{shape} is \\code{arrow}.",
                                  arrowHeadWidth="Numeric scalar. The width of the arrow head in pixels if \\code{shape} is \\code{fixedArrow}.",
                                  cex.group="Numeric scalar. The font expansion factor for the group-level annotation.",
                                  cex="Numeric scalar. The font expansion factor for item identifiers.",
                                  col.line="Character scalar. The color used for connecting lines between grouped items. Defaults to a light gray, but if set to \\code{NULL} the same color as for the first item in the group is used.",
                                  col="Character or integer scalar. The border color for all track items.",
                                  featureAnnotation="Character scalar. Add annotation information to the individual track elements. This can be a value in \\code{id}, \\code{group} or \\code{feature}. Defaults to \\code{id}. Only works if \\code{showFeatureId} is not \\code{FALSE}.",
                                  fontcolor.group="Character or integer scalar. The font color for the group-level annotation.",
                                  fontcolor="Character or integer scalar. The font color for item identifiers.",
                                  fontface.group="Numeric scalar. The font face for the group-level annotation.",
                                  fontface="Integer scalar. The font face for item identifiers.",
                                  fontfamily.group="Character scalar. The font family for the group-level annotation.",
                                  fontfamily="Character scalar. The font family for item identifiers.",
                                  fontsize.group="Numeric scalar. The font size for the group-level annotation.",
                                  fontsize="Numeric scalar. The font size for item identifiers.",
                                  groupAnnotation="Character scalar. Add annotation information as group labels. This can be a value in \\code{id}, \\code{group} or \\code{feature}. Defaults to \\code{group}. Only works if \\code{showId} is not \\code{FALSE}.",
                                  just.group="Character scalar. the justification of group labels. Either \\code{left}, \\code{right}, \\code{above} or \\code{below}.",
                                  lex="Numeric scalar. The line expansion factor for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lineheight="Numeric scalar. The font line height for item identifiers.",
                                  lty="Character or integer scalar. The line type for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lwd="Integer scalar. The line width for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  mergeGroups="Logical scalar. Merge fully overlapping groups if \\code{collapse==TRUE}.",
                                  min.height="Numeric scalar. The minimum range height in pixels to display. All ranges are expanded to this size in order to avoid rendering issues.  See \\code{\\link{collapsing}} for details. For feathered bars indicating the strandedness of grouped items this also controls the height of the arrow feathers.",
                                  min.width="Numeric scalar. The minimum range width in pixels to display. All ranges are expanded to this size in order to avoid rendering issues. See \\code{\\link{collapsing}} for details.",
                                  rotation.group="Numeric scalar. The degree of text rotation for group labels.",
                                  rotation.item="Numeric scalar. The degree of text rotation for item identifiers.",
                                  shape="Character scalar. The shape in which to display the track items. Currently only \\code{box}, \\code{arrow}, \\code{fixedArrow}, \\code{ellipse}, and \\code{smallArrow} are implemented.",
                                  showFeatureId="Logical scalar. Control whether to plot the individual track item identifiers.",
                                  showId="Logical scalar. Control whether to annotate individual groups.",
                                  showOverplotting="Logical scalar. Use a color gradient to show the amount of overplotting for collapsed items. This implies that \\code{collapse==TRUE}",
                                  size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                  fill="Character or integer scalar. The fill color for untyped items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details."

                                  ),

                DetailsAnnotationTrack=c(

                                         details.minWidth="Numeric scalar. The minium width in pixels for a details panel, if less space is available no details are plotted.",
                                         details.ratio="Numeric scalar. By default, the plotting method tries to fill all available space of the details panel tiles. Depending on the dimensions of your plot and the number of tiles this may lead to fairly stretched plots. Restricting the ration of width over height can help to fine tune for somewhat more sane graphics in these cases. Essentially this adds some white space in between individual tiles to force the desired ratio. Together with the \\code{size} and \\code{details.size} arguments, which control the vertical extension of the whole track and of the details section, this allows for some fairly generic resizing of the tiles.",
                                         detailsBorder.col="Character or integer scalar. Line color of the border.",
                                         detailsBorder.fill="Character or integer scalar. Background color of the border.",
                                         detailsBorder.lty="Character or integer scalar. Line type of the border around each details panel.",
                                         detailsBorder.lwd="Integer scalar. Line width of the border.",
                                         detailsConnector.cex="Numeric scalar. Relative size of the connector's end points.",
                                         detailsConnector.col="Character or integer scalar. Color of the line connecting the \\code{AnnotstionTrack} item with its details panel.",
                                         detailsConnector.lty="Character or integer scalar. Type of connecting line.",
                                         detailsConnector.lwd="Integer scalar. Line width of the connector.",
                                         detailsConnector.pch="Integer scalar. Type of the connector's ends.",
                                         detailsFunArgs="List.Additional arguments that get passed on the the details plotting function.",
                                         details.size="Numeric scalar. The fraction of vertical space of the track used for the details section.",
                                         groupDetails="Logial scalar. Plot details for feature groups rather than for individual features."

                                         ),

                GeneRegionTrack=c(

                                  alpha="Numeric scalar between 0 and 1. The opacity of the plotting elements, if supported by the device.",
                                  arrowHeadMaxWidth="Numeric scalar. The maximum width of the arrow head in pixels if \\code{shape} is \\code{arrow}.",
                                  arrowHeadWidth="Numeric scalar. The width of the arrow head in pixels if \\code{shape} is \\code{fixedArrow}.",
                                  cex.group="Numeric scalar. The font expansion factor for the group-level annotation.",
                                  cex="Numeric scalar. The font expansion factor for item identifiers.",
                                  col="Character or integer scalar. The border color for all track items. Defaults to using the same color as in \\code{fill}, also taking into account different track \\code{features}.",
                                  collapseTranscripts="Logical or character scalar. Can be one in \\code{gene}, \\code{longest}, \\code{shortest} or \\code{meta}. Merge all transcripts of the same gene into one single gene model. In the case of \\code{gene} (or \\code{TRUE}), this will only keep the start location of the first exon and the end location of the last exon from all transcripts of the gene. For \\code{shortest} and \\code{longest}, only the longest or shortest transcript model is retained. For \\code{meta}, a meta-transcript containing the union of all exons is formed (essentially identical to the operation \\code{reduce(geneModel)}).",
                                  exonAnnotation="Character scalar. Add annotation information to the individual exon models. This can be a value in \\code{symbol}, \\code{gene}, \\code{transcript}, \\code{exon} or \\code{feature}. Defaults to \\code{exon}. Only works if \\code{showExonId} is not \\code{FALSE}.",
                                  fill="Character or integer scalar. The fill color for untyped items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  fontcolor.group="Character or integer scalar. The font color for the group-level annotation.",
                                  fontcolor="Character or integer scalar. The font color for item identifiers.",
                                  fontface.group="Numeric scalar. The font face for the group-level annotation.",
                                  fontface="Integer scalar. The font face for item identifiers.",
                                  fontfamily.group="Character scalar. The font family for the group-level annotation.",
                                  fontfamily="Character scalar. The font family for item identifiers.",
                                  fontsize.group="Numeric scalar. The font size for the group-level annotation.",
                                  fontsize="Numeric scalar. The font size for item identifiers.",
                                  geneSymbols="Logical scalar. Use human-readable gene symbols or gene IDs for the transcript annotation.",
                                  lex="Numeric scalar. The line expansion factor for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lineheight="Numeric scalar. The font line height for item identifiers.",
                                  lty="Character or integer scalar. The line type for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  lwd="Integer scalar. The line width for all track items. This is also used to connect grouped items. See \\code{\\link{grouping}} for details.",
                                  min.distance="Numeric scalar. The minimum pixel distance before collapsing range items, only if \\code{collapse==TRUE}. See \\code{\\link{collapsing}} for details. Note that a value larger than 0 may lead to UTR regions being merged to CDS regions, which in most cases is not particularly useful.",
                                  min.width="Numeric scalar. The minimum range width in pixels to display. All ranges are expanded to this size in order to avoid rendering issues. See \\code{\\link{collapsing}} for details.",
                                  rotation="Numeric scalar. The degree of text rotation for item identifiers.",
                                  shape="Character scalar. The shape in which to display the track items. Currently only \\code{box}, \\code{arrow}, \\code{ellipse}, and \\code{smallArrow} are implemented.",
                                  showExonId="Logical scalar. Control whether to plot the individual exon identifiers.",
                                  showId="Logical scalar. Control whether to annotate individual groups.",
                                  showOverplotting="Logical scalar. Use a color gradient to show the amount of overplotting for collapsed items. This implies that \\code{collapse==TRUE}",
                                  size="Numeric scalar. The relative size of the track. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                  thinBoxFeature="Character vector. A listing of feature types that should be drawn with thin boxes. Typically those are non-coding elements.",
                                  transcriptAnnotation="Character scalar. Add annotation information as transcript labels. This can be a value in \\code{symbol}, \\code{gene}, \\code{transcript}, \\code{exon} or \\code{feature}. Defaults to \\code{symbol}. Only works if \\code{showId} is not \\code{FALSE}."

                                  ),


                BiomartGeneRegionTrack=c(

                                         "C_segment"="Character or integer scalar. Fill color for annotation objects of type 'C_segment'.",
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
                                         "V_segment"="Character or integer scalar. Fill color for annotation objects of type 'V_segment'.",
                                         "verbose"="Logical scalar. Report data loading events from Bioamart or retrieval from cache."
                                         ),

                AlignedReadTrack=c(

                                   collapse="collapse overlapping ranges and aggregate the underlying data.",
                                   detail="the amount of detail to plot the data. Either \\code{coverage} to show the coverage only, or \\code{reads} to show individual reads. For large data sets the latter can be very inefficient. Please note that \\code{reads} is only available when the object has been created with option \\code{coverageOnly=FALSE}.",
                                   fill="the fill color for the coverage indicator.",
                                   size="the relative size of the track. Defaults to size selection based on the underlying data. Can be overridden in the \\code{\\link{plotTracks}} function.",
                                   type="the plot type, one or several in \\code{c(\"p\",\"l\", \"b\", \"a\", \"s\", \"g\", \"r\", \"S\", \"smooth\", \"histogram\", \"mountain\", \"polygon\", \"h\", \"boxplot\", \"gradient\", \"heatmap\", \"horizon\")}. See the 'Details' section in \\code{\\linkS4class{DataTrack}} for more information on the individual plotting types."

                                   ),

                SequenceTrack=c(

                                add53="Logical scalar. Add a direction indicator",
                                background.title="Character scalar. Make the title panel background transparent by default.",
                                cex="The character expansion factor for the size of the sequence letters. Together with \\code{fontsize} this determines the final font size and thus the level of plotable details.",
                                col="Character scalar. The color of the line when no indiviual letters can be plotted due to size limitations.",
                                complement="Logical scalar. Plot the sequence complement.",
                                fontcolor="Character vector. The colors used for the 5 possible nucleotides (G, A, T, C, N). Defaults to use colors as defined in the \\code{biovizBase} package.",
                                fontface="Numeric scalar. The face of the font.",
                                fontsize="Numeric scalar. Controls the size of the sequence letters and thus also the level of plotable details.",
                                lwd="Numeric scalar. The width of the line when no indiviual letters can be plotted due to size limitations.",
                                min.width="Numeric scalar. The minimum width in pixels of the colored boxes that are drawn when no indiviual letters can be plotted due to size limitations. If the horizontal space that a single base occupies is smaller than this value, only a horizontal line is drawn to indicate the presence of a sequence.",
                                noLetters="Logical scalar. Always plot colored boxes (or a line) regardles of the available space.",
                                rotation="Numeric scalar. The rotation angle for each individual letter in the sequence.",
                                showTitle="Logical scalar. Do not show a title panel by default.",
                                size="Numeric scalar. The size of the track item. Defaults to auto-detect the size based on the other parameter settings."

                                ),

                HighlightTrack=c(

                                 col="Integer or character scalar. The boder color for the highlighting regions.",
                                 fill="Integer or character scalar. The fill color for the highlighting regions.",
                                 inBackground="Logical scalar. Place the box in the background or in the foreground."

                                 ),

                AlignmentsTrack=c(

                                  alpha.mismatch="Numeric scalar between 0 and 1. The transparency of the mismatch base information.",
                                  alpha.reads="Numeric scalar between 0 and 1. The transparency of the individual read icons. Can be used to indicate overlapping regions in read pairs. Only on supported devices.",
                                  cex.mismatch="Numeric Scalar. The character expansion factor for the mismatch base letters.",
                                  col.coverage="Integer or character scalar. The line color for the coverage profile.",
                                  col.gap="Integer or character scalar. The color of the line that is bridging the gap regions in gapped alignments.",
                                  col.mates="Integer or character scalar. The color of the line that is connecting two paired reads.",
                                  col.mismatch="Integer or character scalar. The box color around mismatch bases.",
                                  col.reads="Integer or character scalar. The box color around reads.",
                                  col.sashimi="Integer or character scalar. The line color for sashimi plots.",
                                  col="Integer or character scalar. The default color of all line elements.",
                                  collapse="Logical scalar. Do not perform any collapsing of overlapping elements. Currently not supported.",
                                  coverageHeight="Numeric scalar. The height of the coverage region of the track. Can either be a value between 0 and 1 in which case it is taken as a relative height, or a positive value greater 1 in which case it is interpreted as pixels.",
                                  fill.coverage="Integer or character scalar. The fill color for the coverage profile.",
                                  fill.reads="Integer or character scalar. The fill color for the read icons.",
                                  fill="Integer or character scalar. The default fill color of all plot elements.",
                                  fontface.mismatch="Integer scalar. The font face for mismatch bases.",
                                  lty.coverage="Integer or character scalar. The line type of the coverage profile.",
                                  lty.gap="Integer or character scalar. The type of the line that is bridging the gap regions in gapped alignments.",
                                  lty.mates="Integer or character scalar. The type of the line that is connecting two paired reads.",
                                  lty.mismatch="Integer or character scalar. The box line type around mismatch bases.",
                                  lty.reads="Integer or character scalar. The box line type around mismatch reads.",
                                  lty="Integer or character scalar. The default type of all line elements.",
                                  lwd.coverage="Integer or character scalar. The line width of the coverage profile.",
                                  lwd.gap="Integer scalar. The width of the line that is bridging the gap regions in gapped alignments.",
                                  lwd.mates="Integer scalar. The width of the line that is connecting two paired reads.",
                                  lwd.mismatch="Integer scalar. The box line width around mismatch bases.",
                                  lwd.reads="Integer scalar. The box line width around reads.",
                                  lwd.sashimiMax="Integer scalar. The maximal width of the line in sashimi plots.",
                                  lwd="Integer scalar. The default width of all line elements.",
                                  noLetters="Logical scalar. Always plot colored boxes for mismatch bases regardles of the available space.",
                                  max.height="Integer scalar. The maximum height of an individual read in pixels. Can be used in combination with \\code{min.height} to control the read and stacking appearance. ",
                                  min.height="Integer scalar. The minimum height of an individual read in pixels. Can be used in combination with \\code{max.height} to control the read and stacking appearance.",
                                  minCoverageHeight="Integer scalar. The minimum height of the coverage section. Uselful in combination with a relative setting of \\code{coverageHeight}.",
                                  minSashimi="Integer scalar. The minimum height of the sashimi section. Uselful in combination with a relative setting of \\code{sashimiHeight}.",
                                  showMismatches="Logical scalar. Add mismatch information, either as individual base letters or using color coded bars. This implies that the reference sequence has been provided, either to the class constructor or as part of the track list.",
                                  sashimiFilter="GRanges object. Only junctions which overlap equally with \\code{sashimiFilter} GRanges are shown. Default \\code{NULL}, no filtering.",
                                  sashimiFilterTolerance="Integer scalar. Only used in combination with \\code{sashimiFilter}. It allows to include junctions whose starts/ends are within specified distance from \\code{sashimiFilter} GRanges. This is useful for cases where the aligner did not place the junction reads precisely. Default \\code{0L} , no tolerance.",
                                  sashimiHeight="Integer scalar. The height of the sashimi part of the track. Can either be a value between 0 and 1 in which case it is taken as a relative height, or a positive value greater 1 in which case it is interpreted as pixels.",
                                  sashimiScore="Integer scalar. The minimum number of reads supporting the junction.",
                                  sashimiStrand="Integer scalar. Only reads which have the specified strand are considered to count the junctions.",
                                  size="Numeric scalar. The size of the track. Defaults to automatic sizing.",
                                  transformation: "Function. Applied to the coverage vector prior to plotting. The function should accept exactly one input argument and its return value needs to be a numeric Rle of identical length as the input data."
                                  type="Character vactor. The type of information to plot. For \\code{coverage} a coverage plot, potentially augmented by base mismatch information, for \\code{sashimi} a sashimi plot, showing the juctions, and for \\code{pileup} the pileups of the individual reads. Theese three can be combined."

                                  )
                )


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
