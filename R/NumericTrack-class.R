#' @include RangeTrack-class.R
NULL

## NumericTrack-class --------------------------------------------------------

#' NumericTrack class and methods
#'
#' The virtual parent class for all track items in the Gviz package designed to
#' contain numeric data. This class merely exists for dispatching purpose.
#'
#' @template GdObject-class_slot
#' @template RangeTrack-class_slot
#'
#' @template NumericTrack_param
#'
#' @name NumericTrack-class
#'
#' @return A virtual class: No objects may be created from it.
#'
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @exportClass NumericTrack
setClass("NumericTrack",
    representation = representation("VIRTUAL"),
    prototype = prototype(
        name = "NumericTrack",
        dp = DisplayPars()
    ),
    contains = "RangeTrack"
)


## NumericTrack Methods drawAxis ---------------------------------------------

#' @describeIn NumericTrack-class add a y-axis to the title panel of a track.
#' @importFrom grDevices extendrange colorRampPalette
#' @export
setMethod("drawAxis", signature(GdObject = "NumericTrack"), function(GdObject, from, to, ...) {
    type <- match.arg(.dpOrDefault(GdObject, "type", "p"), .PLOT_TYPES, several.ok = TRUE)
    yvals <- values(GdObject)
    ylim <- .dpOrDefault(GdObject, "ylim", if (!is.null(yvals) && length(yvals)) {
        range(yvals, na.rm = TRUE, finite = TRUE)
    } else {
        c(-1, 1)
    })
    if (diff(ylim) == 0) {
        ylim <- ylim + c(-1, 1)
    }
    hSpaceAvail <- vpLocation()$isize["width"] / 6
    yscale <- extendrange(r = ylim, f = 0.05)
    col <- .dpOrDefault(GdObject, "col.axis", "white")
    acex <- .dpOrDefault(GdObject, "cex.axis")
    acol <- .dpOrDefault(GdObject, "col.axis", "white")
    at <- .dpOrDefault(GdObject, "yTicksAt", pretty(yscale))
    at <- at[which(at >= sort(ylim)[1] & at <= sort(ylim)[2])]
    if (!length(at)) {
        at <- sort(ylim)[c(1, 2)]
    }
    if (is.null(acex)) {
        vSpaceNeeded <- max(as.numeric(convertWidth(stringHeight(at), "inches"))) * length(at) * 1.5
        hSpaceNeeded <- max(as.numeric(convertWidth(stringWidth(at), "inches")))
        vSpaceAvail <- abs(diff(range(at))) / abs(diff(yscale)) * vpLocation()$isize["height"]
        acex <- max(0.6, min(vSpaceAvail / vSpaceNeeded, hSpaceAvail / hSpaceNeeded))
    }
    nlevs <- max(1, nlevels(factor(.dpOrDefault(GdObject, "groups"))))
    if (any(type %in% c("heatmap", "horizon")) && .dpOrDefault(GdObject, "showSampleNames", FALSE)) {
        groups <- .dpOrDefault(GdObject, "groups")
        sn <- if (is.null(groups)) rownames(values(GdObject)) else rev(unlist(split(rownames(values(GdObject)), factor(groups))))
        cex.sn <- .dpOrDefault(GdObject, "cex.sampleNames", acex)
        col.cn <- .dpOrDefault(GdObject, "col.sampleNames", "white")
        wd <- max(as.numeric(convertWidth(stringWidth(sn) + unit(10, "points"), "npc"))) * cex.sn
        samNames <- viewport(x = 1, width = wd, just = 1, yscale = c(-0.05, 1.05))
        pushViewport(samNames)
        nr <- nrow(values(GdObject))
        if (nr > 1) {
            yy <- head(seq(0.05, 0.95, len = nr + 1), -1)
            yy <- yy + diff(yy)[[1]] / 2
        } else {
            yy <- 0.5
        }
        grid.text(x = rep(0.5, nr), y = yy, label = rev(sn), just = 0.5, gp = gpar(cex = cex.sn, col = col.cn))
        popViewport(1)
        samAxis <- viewport(x = 1 - wd, width = 1 - wd, just = 1)
        pushViewport(samAxis)
        on.exit(popViewport(1))
    }
    ## if any of the types are gradient or heatmap we want the gradient scale
    if (any(type %in% c("gradient", "heatmap")) && .dpOrDefault(GdObject, "showColorBar", TRUE)) {
        ## viewport to hold the color strip
        shift <- ifelse(all(type %in% c("gradient", "heatmap")), 1, 0)
        pcols <- .getPlottingFeatures(GdObject)
        ncolor <- .dpOrDefault(GdObject, "ncolor", 100)
        vpAxisCont <- viewport(x = unit(1, "npc") - unit(2 - shift, "points"), width = unit(1, "npc") - unit(2 - shift, "points"), just = 1)
        pushViewport(vpAxisCont)
        for (i in seq_len(nlevs)) {
            ## create color palette
            cr <- c("white", pcols$col[i])
            if (nlevs < 2) {
                cr <- .dpOrDefault(GdObject, "gradient", cr)
            }
            palette <- colorRampPalette(cr)(ncolor + 5)[-seq_len(5)]
            pshift <- ifelse(i == nlevs, 1 - shift, 0)
            vpTitleAxis <- viewport(
                x = unit(1, "npc") - unit(4 * (i - 1), "points"), width = unit(4 + pshift, "points"),
                yscale = yscale, just = 1
            )
            pushViewport(vpTitleAxis)
            ## draw a rectangle for each color
            if (all(type %in% c("gradient", "heatmap"))) {
                if (i == nlevs) {
                    suppressWarnings(grid.yaxis(gp = gpar(col = acol, cex = acex), at = at))
                }
                grid.rect(
                    y = unit(seq(ylim[1], ylim[2], length.out = ncolor + 1), "native")[-(ncolor + 1)], x = unit(0, "npc") - unit(1, "points"),
                    width = 1, height = 1 / ncolor, gp = gpar(fill = palette, lty = 0), just = c("left", "bottom")
                )
            } else {
                grid.rect(
                    y = unit(seq(ylim[1], ylim[2], length.out = ncolor + 1), "native")[-(ncolor + 1)], x = 0,
                    width = 1, height = 1 / ncolor, gp = gpar(fill = palette, lty = 0), just = c("left", "bottom")
                )
                if (i == nlevs) {
                    suppressWarnings(grid.yaxis(gp = gpar(col = acol, cex = acex), at = at))
                    grid.lines(x = c(0, 0), y = ylim, gp = gpar(col = acol), default.units = "native")
                }
            }
            popViewport(1)
        }
        popViewport(1)
    } else {
        vpTitleAxis <- viewport(x = 0.95, width = 0.2, yscale = yscale, just = 0)
        pushViewport(vpTitleAxis)
        suppressWarnings(grid.yaxis(gp = gpar(col = acol, cex = acex), at = at))
        grid.lines(x = c(0, 0), y = ylim, gp = gpar(col = acol), default.units = "native")
        popViewport(1)
    }
})

## NumericTrack Methods drawGrid ---------------------------------------------

#' @describeIn NumericTrack-class superpose a grid on top of a track.
setMethod("drawGrid", signature(GdObject = "NumericTrack"), function(GdObject, from, to) {
    if (.dpOrDefault(GdObject, "grid", FALSE)) {
        vals <- score(GdObject)
        ylim <- .dpOrDefault(GdObject, "ylim", range(vals, na.rm = TRUE, finite = TRUE))
        if (diff(ylim)) {
            pushViewport(dataViewport(xData = c(from, to), yData = ylim, extension = c(0, 0.1), clip = TRUE))
            panel.grid(
                h = .dpOrDefault(GdObject, "h", -1), v = .dpOrDefault(GdObject, "v", -1),
                col = .dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty = .dpOrDefault(GdObject, "lty.grid", 1),
                lwd = .dpOrDefault(GdObject, "lwd.grid", 1)
            )
            popViewport(1)
        }
    }
})
