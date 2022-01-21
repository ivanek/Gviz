#' @include GdObject-class.R
NULL

## GenomeAxisTrack Class -----------------------------------------------------

#' GenomeAxisTrack class and methods
#'
#' A class representing a customizable genomic axis.
#'
#'
#' A \code{GenomeAxisTrack} can be customized using the familiar display
#' parameters. By providing a \code{GRanges} or \code{IRanges} object to the
#' constructor, ranges on the axis can be further highlighted.
#'
#' With the \code{scale} display parameter, a small scale indicator can be
#' shown instead of the entire genomic axis. The scale can either be provided
#' as a fraction of the plotting region (it will be rounded to the nearest
#' human readable absolute value) or as an absolute value and is always
#' displayed in bp, kb, mb or gb units. Note that most display parameters for
#' the \code{GenomeAxisTrack} are ignored when a scale is used instead of the
#' full axis. In particular, only the parameters \code{exponent}, \code{alpha},
#' \code{lwd}, \code{col}, \code{cex}, \code{distFromAxis} and \code{labelPos}
#' are used.
#'
#' @template GenomeAxisTrack-class_param
#'
#' @name GenomeAxisTrack-class
#'
#' @return The return value of the constructor function is a new object of class
#' \code{GenomeAxisTrack}.
#'
#' Objects can be created using the constructor function \code{GenomeAxisTrack}.
#'
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @examples
#' ## Construct object
#' axTrack <- GenomeAxisTrack(
#'     name = "Axis",
#'     range <- IRanges(start = c(100, 300, 800), end = c(150, 400, 1000))
#' )
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(axTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#'
#' ## Plotting
#' plotTracks(axTrack, from = 0, to = 1100)
#'
#' ## Track names
#' names(axTrack)
#' names(axTrack) <- "foo"
#'
#' ## Subsetting and splitting
#' subTrack <- subset(axTrack, from = 0, to = 500)
#' length(subTrack)
#' subTrack[1]
#' split(axTrack, c(1, 1, 2))
#'
#' ## Accessors
#' start(axTrack)
#' end(axTrack)
#' width(axTrack)
#'
#' strand(axTrack)
#'
#' range(axTrack)
#' ranges(axTrack)
#'
#' ## Annotation
#' values(axTrack)
#'
#' ## Grouping
#' group(axTrack)
#'
#' ## HTML image map
#' coords(axTrack)
#' tags(axTrack)
#' axTrack <- plotTracks(axTrack)$foo
#' coords(axTrack)
#' tags(axTrack)
#'
#' ## adding an axis to another track
#' data(cyp2b10)
#' grTrack <- GeneRegionTrack(
#'     start = 26682683, end = 26711643,
#'     rstart = cyp2b10$start, rends = cyp2b10$end, chromosome = 7, genome = "mm9",
#'     transcript = cyp2b10$transcript, gene = cyp2b10$gene, symbol = cyp2b10$symbol,
#'     name = "Cyp2b10", strand = cyp2b10$strand
#' )
#'
#' plotTracks(list(grTrack, GenomeAxisTrack()))
#' plotTracks(list(grTrack, GenomeAxisTrack(scale = 0.1)))
#' plotTracks(list(grTrack, GenomeAxisTrack(scale = 5000)))
#' plotTracks(list(grTrack, GenomeAxisTrack(scale = 0.5, labelPos = "below")))
#' @exportClass GenomeAxisTrack
setClass("GenomeAxisTrack",
    contains = "GdObject",
    representation = representation(range = "GRanges"),
    prototype(
        range = GRanges(),
        name = "GenomeAxisTrack",
        dp = DisplayPars(
            add35 = FALSE,
            add53 = FALSE,
            background.title = "transparent",
            cex.id = 0.7,
            cex = 0.8,
            col.border.title = "transparent",
            lwd.border.title = 1,
            col.id = "white",
            col.range = "cornsilk4",
            distFromAxis = 1,
            exponent = NULL,
            fill.range = "cornsilk3",
            fontcolor = "#808080",
            fontsize = 10,
            labelPos = "alternating",
            littleTicks = FALSE,
            lwd = 2,
            scale = NULL,
            showId = FALSE,
            showTitle = FALSE,
            ticksAt = NULL,
            size = NULL,
            col = "darkgray"
        )
    )
)
## Initialize ----------------------------------------------------------------

## Only pass on the stuff to the GdObject initializer

#' @describeIn GenomeAxisTrack-class Intialize.
#' @export
setMethod("initialize", "GenomeAxisTrack", function(.Object, range, ids, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "GenomeAxisTrack")
    if (missing(range) || is.null(range)) {
        range <- GRanges()
    }
    if (is(range, "IRanges")) {
        range <- GRanges(ranges = range, seqnames = "dummy", id = ids)
    }
    .Object@range <- range
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

#' @describeIn GenomeAxisTrack-class Constructor
#' @export
GenomeAxisTrack <- function(range = NULL, name = "Axis", id, ...) {
    if (missing(id)) {
        id <- names(range)
    }
    new("GenomeAxisTrack", name = name, range = range, id = id, ...)
}

## General accessors ---------------------------------------------------------

#' @describeIn GenomeAxisTrack-class return the genomic coordinates for the
#' track along with all additional annotation information as an object of
#' class `GRanges`.
setMethod("ranges", "GenomeAxisTrack", function(x) x@range)
setReplaceMethod("ranges", "GenomeAxisTrack", function(x, value) {
    x@range <- value
    return(x)
})

#' @describeIn GenomeAxisTrack-class return the genomic coordinates for the
#'  track as an object of class `IRanges`.
#'  @export
setMethod("range", "GenomeAxisTrack", function(x) ranges(x@range))

#' @describeIn GenomeAxisTrack-class return the start coordinates of the track
#' items.
#' @export
setMethod("start", "GenomeAxisTrack", function(x) if (length(x)) start(range(x)) else NULL)

#' @describeIn GenomeAxisTrack-class replace the start coordinates of the track
#' items.
#' @export
setReplaceMethod("start", "GenomeAxisTrack", function(x, value) {
    start(x@range) <- value
    return(x)
})

#' @describeIn GenomeAxisTrack-class return the end coordinates of the track
#' items.
#' @export
setMethod("end", "GenomeAxisTrack", function(x) if (length(x)) end(range(x)) else NULL)

#' @describeIn GenomeAxisTrack-class replace the end coordinates of the track
#' items.
#' @export
setReplaceMethod("end", "GenomeAxisTrack", function(x, value) {
    end(x@range) <- value
    return(x)
})

#' @describeIn GenomeAxisTrack-class return the with of the track items in
#' genomic coordinates.
#' @export
setMethod("width", "GenomeAxisTrack", function(x) if (length(x)) as.integer(width(range(x))) else NULL)

#' @describeIn GenomeAxisTrack-class return the number of items stored in the ranges slot.
#' @export
setMethod("length", "GenomeAxisTrack", function(x) length(ranges(x)))

#' @describeIn GenomeAxisTrack-class return all additional annotation information
#' except for the genomic coordinates for the track items.
#' @export
setMethod("values", "GenomeAxisTrack", function(x) as.data.frame(values(ranges(x))))

#' @describeIn GenomeAxisTrack-class return a vector of strand specifiers for
#' all track items, in the form '+' for the Watson strand, '-' for the Crick
#' strand or '*' for either of the two.
#' @export
setMethod("strand", "GenomeAxisTrack", function(x) as.character(strand(ranges(x))))

## Annotation Accessors ------------------------------------------------------

## Collapse  -----------------------------------------------------------------

#' @describeIn GenomeAxisTrack-class preprocess the track before plotting.
#' This will collapse overlapping track items based on the available resolution
#' and increase the width and height of all track objects to a minimum value
#' to avoid rendering issues. See collapsing for details.
#' @keywords internal
setMethod(
    "collapseTrack", signature(GdObject = "GenomeAxisTrack"),
    function(GdObject, min.width = 1, min.distance = 0, collapse = TRUE, diff = .pxResolution(coord = "x"), xrange) {

        ## Collapse overlapping ranges (less than minXDist space between them) including the associated attributes using
        ## "|" as separator. For both "strand" and "feature" we take the first available entry, which is not optimal but
        ## seems to be the sanest thing to do here...
        if (collapse) {
            GdObject <- GdObject[order(range(GdObject))]
            r <- ranges(GdObject)
            minXDist <- min.distance * diff
            r <- reduce(r, min.gapwidth = minXDist)
        }
        r <- .resize(r, min.width, diff)
        GdObject@range <- r
        return(GdObject)
    }
)

## Subset --------------------------------------------------------------------

#' @describeIn GenomeAxisTrack-class subset the items in the `GenomeAxisTrack` object.
#' This is essentially similar to subsetting of the `GRanges` object in the
#' `range` slot. For most applications, the subset method may be more appropriate.
#' @export
setMethod("[", signature(x = "GenomeAxisTrack"), function(x, i, j, ..., drop = TRUE) {
    x <- .deepCopyPars(x)
    x@range <- x@range[i, , drop = drop]
    return(x)
})

#' @describeIn GenomeAxisTrack-class plot subset all the contained tracks in an
#' `GenomeAxisTrack` by coordinates and sort if necessary.
#' @export
setMethod("subset", signature(x = "GenomeAxisTrack"), function(x, from = NULL, to = NULL, sort = FALSE, ...) {
    if (!length(x)) {
        return(x)
    }
    ranges <- .defaultRange(x, from = from, to = to)
    lsel <- end(x) < ranges["from"]
    rsel <- start(x) > ranges["to"]
    x <- x[!(lsel | rsel), ]
    if (sort) {
        x <- x[order(range(x)), ]
    }
    return(x)
})

## DrawGD --------------------------------------------------------------------

.expLabel <- function(GdObject, tckText, prune = FALSE) {
    tck <- tckText
    exponent <- if (is.null(.dpOrDefault(GdObject, "exponent"))) {
        exp <- 0
        while (all(abs(tck)[abs(tck) > 0] / 10^exp >= 1)) {
            exp <- exp + 3
        }
        exp - 3
    } else {
        max(0, .dpOrDefault(GdObject, "exponent"))
    }
    if (exponent > 0) {
        tckText <- tckText / (10^exponent)
    }
    if (prune) {
        tmp <- as.character(tckText)
        count <- max(nchar(gsub("*.\\.", "", tmp)))
        while (count > 1 && !any(duplicated(round(tckText, count)))) {
            count <- count - 1
        }
        tckText <- round(tckText, count + 1)
    }
    ## alternative might be to show "bp" also for small integers
    ## "0" = sprintf("%s bp", tckText),
    return(switch(as.character(exponent),
        "0" = sprintf("%i", as.integer(tckText)),
        "3" = sprintf("%s kb", tckText),
        "6" = sprintf("%s mb", tckText),
        "9" = sprintf("%s gb", tckText),
        lapply(tckText, function(x) bquote(paste(.(x), " ", 10^.(exponent))))
    ))
}

#' @describeIn GenomeAxisTrack-class plot the object to a graphics device.
#' The return value of this method is the input object, potentially updated
#' during the plotting operation. Internally, there are two modes in which the
#' method can be called. Either in 'prepare' mode, in which case no plotting is
#' done but the object is preprocessed based on the available space, or in
#' 'plotting' mode, in which case the actual graphical output is created.
#' Since subsetting of the object can be potentially costly, this can be
#' switched off in case subsetting has already been performed before or
#' is not necessary.
#'
#' @export
setMethod("drawGD", signature("GenomeAxisTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, subset = TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if ((is.logical(debug) && debug) || debug == "prepare") {
        browser()
    }
    ## Nothing to do if the coordinate width is 0, so we can quit right away
    imageMap(GdObject) <- NULL
    if ((maxBase - minBase) == 0) {
        return(invisible(GdObject))
    }
    ## We start by setting up the drawing canvas
    if (subset) {
        GdObject <- subset(GdObject, from = minBase, to = maxBase)
    }
    xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
    ## pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(-1, 1), extension=0))
    pushViewport(dataViewport(xscale = xscale, yscale = c(-1, 1), extension = 0))
    ## Create a useful data range for the axis
    pres <- .pxResolution()
    curVp <- vpLocation()
    cex <- .dpOrDefault(GdObject, "cex", 0.8)
    lwd <- .dpOrDefault(GdObject, "lwd", 1)
    fontface <- .dpOrDefault(GdObject, "fontface", 1)
    add53 <- .dpOrDefault(GdObject, "add53", FALSE)
    add35 <- .dpOrDefault(GdObject, "add35", FALSE)
    lcex <- cex * 0.75
    textYOff <- pres["y"] * 3
    textXOff <- pres["x"] * 2
    endMargin <- if (add53 || add35) {
        abs(as.numeric(convertWidth(stringWidth("5'"), "native")) * lcex) + (textXOff * 2)
    } else {
        pres["x"] * 5
    }
    axRange <- c(minBase + endMargin, maxBase - endMargin)
    ## We want fixed vertical sizes for axis tracks to avoid akward stretching effects.
    color <- .dpOrDefault(GdObject, "col", "darkgray")[1]
    littleTicks <- .dpOrDefault(GdObject, "littleTicks", FALSE)
    dfact <- max(1, .dpOrDefault(GdObject, "distFromAxis", 1))
    labelPos <- .dpOrDefault(GdObject, "labelPos", "alternating")
    lwdAdd <- (lwd - 1) / 2
    tickHeight <- (ifelse(littleTicks, 2, 1) * 3 * dfact + lwdAdd) * pres["y"]
    ids <- as.character(values(GdObject)$id)
    showIds <- .dpOrDefault(GdObject, "showId", FALSE) && !is.null(ids) && !all(ids == "")
    rcex <- .dpOrDefault(GdObject, "cex.id", 0.7)
    rcol <- .dpOrDefault(GdObject, "col.id", "white")
    sep <- (if (length(GdObject)) {
        if (showIds) {
            max(1.5, ((max(as.numeric(convertHeight(stringHeight(ids), "native")) * rcex) + textYOff) / pres["y"]) / 2)
        } else {
            1.5
        }
    } else {
        1
    }) + lwdAdd
    pyOff <- pres["y"] * sep
    ## In prepare mode we just want to figure out the optimal size
    if (prepare) {
        nsp <- if (is.null(.dpOrDefault(GdObject, "scale"))) {
            (sum(tickHeight, pyOff * 2, textYOff * 2 + (as.numeric(convertHeight(stringHeight("1"), "native")) / 2) * cex) * 2 * 1.3) / pres["y"]
        } else {
            labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
            if (labelPos %in% c("above", "below")) {
                (sum(tickHeight, pyOff * 2 + (as.numeric(convertHeight(stringHeight("1"), "native")) / 2) * cex) * 2) / pres["y"]
            } else {
                (sum(tickHeight, pyOff * 2 + (as.numeric(convertHeight(stringHeight("1"), "native")) / 2) * cex)) / pres["y"]
            }
        }
        displayPars(GdObject) <- list("neededVerticalSpace" = nsp)
        popViewport(1)
        return(invisible(GdObject))
    }
    if ((is.logical(debug) && debug) || debug == "draw") {
        browser()
    }
    ## Plot range if there is any
    alpha <- .dpOrDefault(GdObject, "alpha", 1)

    ## in "scale" mode we just plot a simple scale and return ...
    scaleLen <- .dpOrDefault(GdObject, "scale")
    if (!is.null(scaleLen)) {
        len <- (maxBase - minBase + 1)
        if (scaleLen > len) {
            warning(sprintf("scale (%d) cannot be larger than plotted region %d - setting to ~5%\n", scaleLen, len))
            scaleLen <- 0.05
        }
        xoff <- len * 0.03 + minBase
        labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
        if (scaleLen <= 1 && scaleLen > 0) { # calculate and round the scale
            scaleLen <- len * scaleLen
            ex <- .dpOrDefault(GdObject, "exponent", floor(log10(scaleLen)))
            v <- round(scaleLen, -ex)
            if (v == 0) v <- scaleLen
        } else { # if the scale is an absolute value don't round
            ex <- .dpOrDefault(GdObject, "exponent", floor(log10(scaleLen)))
            v <- scaleLen
        }

        ## work out exponent/unit
        label <- .expLabel(GdObject, v)
        grid.lines(x = c(xoff, v + xoff), y = c(0, 0), default.units = "native", gp = gpar(col = color, lwd = lwd, alpha = alpha))
        grid.segments(
            x0 = c(xoff, v + xoff), y0 = c(0 - tickHeight, 0 - tickHeight),
            x1 = c(xoff, v + xoff), y1 = c(tickHeight, tickHeight, tickHeight),
            default.units = "native", gp = gpar(col = color, lwd = lwd, alpha = alpha)
        )
        z <- len * 0.01
        if (labelPos == "below") {
            grid.text(
                label = if (is.character(label)) label else label[[1]], x = xoff + v / 2, y = 0 - (tickHeight / 1.5 * dfact),
                just = c("center", "top"), gp = gpar(alpha = alpha, col = color, cex = cex, fontface = fontface), default.units = "native"
            )
        } else if (labelPos == "above") {
            grid.text(
                label = if (is.character(label)) label else label[[1]], x = xoff + v / 2, y = tickHeight / 1.5 * dfact, just = c("center", "bottom"),
                gp = gpar(alpha = alpha, col = color, cex = cex, fontface = fontface), default.units = "native"
            )
        } else {
            grid.text(
                label = if (is.character(label)) label else label[[1]], x = v + xoff + z, y = 0, just = c("left", "center"),
                gp = gpar(alpha = alpha, col = color, cex = cex, fontface = fontface), default.units = "native"
            )
        }
        popViewport(1)
        return(invisible(GdObject))
    }

    GdObject <- GdObject[end(GdObject) > axRange[1] & start(GdObject) < axRange[2]]
    if (length(GdObject)) {
        rfill <- .dpOrDefault(GdObject, "fill.range", "cornsilk3")
        rcolor <- .dpOrDefault(GdObject, "col.range", "cornsilk4")[1]
        diff <- .pxResolution(coord = "x")
        GdObject <- collapseTrack(GdObject, diff = diff, xrange = c(minBase, maxBase))
        start(GdObject) <- pmax(axRange[1], start(GdObject))
        end(GdObject) <- pmin(axRange[2], end(GdObject))
        coords <- cbind(start(GdObject), -0.1, end(GdObject), 0.1)
        grid.rect(
            x = start(GdObject), y = -pyOff, width = width(GdObject), height = pyOff * 2,
            default.units = "native", just = c("left", "bottom"), gp = gpar(col = rcolor, fill = rfill, alpha = alpha)
        )
        vals <- values(GdObject)
        if (showIds) {
            grid.text(ids,
                x = start(GdObject) + width(GdObject) / 2, y = 0,
                gp = gpar(col = rcol, cex = rcex, fontface = fontface),
                default.units = "native", just = c("center", "center")
            )
        }
        ## Calculate the coordinates for the image map
        map <- as.matrix(.getImageMap(coords))
        if (is.null(ids) || length(ids) == 0) {
            ids <- as.character(seq_len(nrow(map)))
        }
        rownames(map) <- make.unique(as.character(ids))
        tags <- lapply(
            list(title = ids, start = as.character(start(GdObject)), end = as.character(end(GdObject))),
            function(x) {
                names(x) <- rownames(map)
                x
            }
        )
        imageMap(GdObject) <- ImageMap(coords = map, tags = tags)
    }
    ## width<1, we can return here, no need for tick marks
    if (abs(diff(axRange)) < 1) {
        popViewport()
        return(invisible(GdObject))
    }
    ## We want two parallel lines with little hooks on the ends
    pyHook <- pres["y"] * (sep + 2 + lwdAdd)
    pxOff <- pres["x"] * 5
    grid.segments(
        x0 = rep(axRange[1], 2), y0 = c(-1, 1) * pyOff, x1 = rep(axRange[2], 2), y1 = c(-1, 1) * pyOff,
        default.units = "native", gp = gpar(col = color, alpha = alpha, lwd = lwd)
    )
    grid.segments(
        x0 = c(axRange[2] - pxOff, axRange[1]), y0 = c(pyHook, -pyOff),
        x1 = c(axRange[2], axRange[1] + pxOff), y1 = c(pyOff, -pyHook),
        default.units = "native", gp = gpar(col = color, alpha = alpha, lwd = lwd)
    )
    ## Here we plot the top level ticks
    tck <- .dpOrDefault(GdObject, "ticksAt", .ticks(axRange))
    tck <- tck[tck < axRange[2] - pxOff * 2 & tck > axRange[1] + pxOff * 2]
    y0t <- rep(c(1, -1) * pyOff, length(tck))[seq_along(tck)]
    y1t <- y0t + rep(c(tickHeight, -tickHeight), length(tck))[seq_along(tck)]
    labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
    y0t <- switch(labelPos,
        "alternating" = y0t,
        "revAlternating" = -y0t,
        "above" = abs(y0t),
        "below" = -abs(y0t),
        "beside" = y0t
    )
    y1t <- switch(labelPos,
        "alternating" = y1t,
        "revAlternating" = -y1t,
        "above" = abs(y1t),
        "below" = -abs(y1t),
        "beside" = y1t
    )
    ## ttck <- if(min(diff(tck))==1) tck+0.5 else tck # to align labels with ticks, shift by 0.5
    ttck <- tck + 0.5 # to align labels with ticks, shift by 0.5
    grid.segments(x0 = ttck, x1 = ttck, y0 = y0t, y1 = y1t, default.units = "native", gp = gpar(col = color, alpha = alpha, lwd = lwd, lineend = "square"))
    ## The top level tick labels
    label <- .expLabel(GdObject, tck)
    ylabs <- y1t + (ifelse(y1t > 0, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"), "native")) / 2) * cex))
    if (is.character(label)) {
        grid.text(
            label = label, x = ttck, y = ylabs, just = c("centre", "centre"),
            gp = gpar(cex = cex, fontface = fontface), default.units = "native"
        )
    } else {
        for (i in seq_along(label)) {
            grid.text(
                label = label[[i]], x = ttck[i], y = ylabs[i], just = c("centre", "centre"),
                gp = gpar(cex = cex, fontface = fontface), default.units = "native"
            )
        }
    }
    ## The second level ticks and labels if necessary
    if (.dpOrDefault(GdObject, "littleTicks", FALSE) && length(tck) > 1) {
        avSpace <- min(diff(tck))
        spaceFac <- 1.8
        spaceNeeded <- min(as.numeric(convertWidth(stringWidth(if (is.character(label)) label else "000000000"), "native")) / 2) * lcex * spaceFac
        nTcks <- (avSpace %/% spaceNeeded)
        if (nTcks %% 2 == 0) {
            nTcks <- nTcks - 1
        }
        btck <- tck
        if (!(minBase %in% btck)) {
            btck <- c(minBase, btck)
        }
        if (!(maxBase %in% btck)) {
            btck <- c(btck, maxBase)
        }
        y0lt <- y1lt <- ltck <- NULL
        for (i in seq_len(length(btck) - 1))
        {
            toFill <- btck[i:(i + 1)]
            ttck <- if (i == 1) rev(toFill[2] - (avSpace / nTcks) * seq_len(nTcks - 1)) else toFill[1] + (avSpace / nTcks) * seq_len(nTcks - 1)
            ltck <- c(ltck, ttck)
            ord <- if (i == 1) {
                if (y0t[1] > 0) c(1, -1) else c(-1, 1)
            } else if (y0t[i - 1] < 0) c(1, -1) else c(-1, 1)
            y0 <- rep(ord * pyOff, length(ttck))[seq_along(ttck)]
            y1 <- y0 + rep(ord * tickHeight / 2, length(ttck))[seq_along(ttck)]
            y0lt <- c(y0lt, switch(labelPos,
                "alternating" = y0,
                "revAlternating" = y0,
                "above" = abs(y0),
                "below" = -abs(y0)
            ))
            y1lt <- c(y1lt, switch(labelPos,
                "alternating" = y1,
                "revAlternating" = y1,
                "above" = abs(y1),
                "below" = -abs(y1)
            ))
        }
        endPadding <- pres["x"] * 15
        sel <- ltck > min(tck, axRange + endPadding) & ltck < max(tck, axRange - endPadding)
        if (length(ltck[sel]) && min(diff(tck)) > nTcks) {
            grid.segments(
                x0 = ltck[sel], x1 = ltck[sel], y0 = y0lt[sel], y1 = y1lt[sel], default.units = "native",
                gp = gpar(col = color, alpha = alpha, lwd = lwd, lineend = "square")
            )
            llabel <- .expLabel(GdObject, ltck[sel], prune = TRUE)
            ytlabs <- y1lt + (ifelse(y1lt > 0, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"), "native")) / 2) * lcex))
            if (is.character(label)) {
                grid.text(
                    label = llabel, x = ltck[sel], y = ytlabs[sel], just = c("centre", "centre"),
                    gp = gpar(cex = lcex, fontface = fontface), default.units = "native"
                )
            } else {
                for (i in seq_along(llabel)) {
                    grid.text(
                        label = llabel[[i]], x = ltck[sel][i], y = ytlabs[sel][i], just = c("centre", "centre"),
                        gp = gpar(cex = lcex, fontface = fontface), default.units = "native"
                    )
                }
            }
        }
    }
    ## The direction indicators
    rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
    p5 <- expression("5'")
    p3 <- expression("3'")
    if (add53) {
        grid.text(
            label = p5, x = min(axRange) - textXOff, y = pyOff,
            just = c(ifelse(rev, "left", "right"), "bottom"), gp = gpar(cex = cex * .75, fontface = fontface),
            default.units = "native"
        )
        grid.text(
            label = p3, x = max(axRange) + textXOff, y = pyOff,
            just = c(ifelse(rev, "right", "left"), "bottom"), gp = gpar(cex = cex * .75, fontface = fontface),
            default.units = "native"
        )
    }
    if (add35) {
        grid.text(
            label = p3, x = axRange[1] - textXOff, y = -pyOff,
            just = c(ifelse(rev, "left", "right"), "top"), gp = gpar(cex = cex * .75, fontface = fontface),
            default.units = "native"
        )
        grid.text(
            label = p5, x = axRange[2] + textXOff, y = -pyOff,
            just = c(ifelse(rev, "right", "left"), "top"), gp = gpar(cex = cex * 0.75, fontface = fontface),
            default.units = "native"
        )
    }
    popViewport()
    return(invisible(GdObject))
})


## Show ----------------------------------------------------------------------

#' @describeIn GenomeAxisTrack-class Show method.
#' @export
setMethod(
    "show", signature(object = "GenomeAxisTrack"),
    function(object) {
        cat(sprintf("Genome axis '%s'\n", names(object)))
        if (.dpOrDefault(object, "add53", FALSE)) {
            cat("5->3 label is set\n")
        }
        if (.dpOrDefault(object, "add35", FALSE)) {
            cat("3->5 label is set\n")
        }
        if (.dpOrDefault(object, "littleTicks", FALSE)) {
            cat("littleTicks label is set\n")
        }
        if (length(object)) {
            cat("There are annotated axis regions:\n")
            print(ranges(object))
        }
    }
)
