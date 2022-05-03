## Plot a list of GdObjects as individual tracks similar to the display on the UCSC genome browser
## Arguments:
##    o trackList: a list of GdObjects
##    o from, to: the plotting range, will be figured out automatically from the tracks if missing
##    o sized: a vector of relative vertical sizes, or NULL to auto-detect
##    o panel.only: don't draw track titles, useful to embed in a lattice-like function, this also implies add=TRUE
##    o extend.right, extend.left: extend the coordinates in 'from' and 'too'
##    o title.width: the expansion factor for the width of the title track
## Value: the function is called for its side-effect of drawing on the graphics device


#' The main plotting function for one or several Gviz tracks.
#'
#' `plotTracks` is the main interface when plotting single track objects,
#' or lists of tracks linked together across the same genomic coordinates.
#' Essentially, the resulting plots are very similar to the graphical output of
#' the UCSC Genome Browser, except for all of the interactivity.
#'
#' Gviz tracks are plotted in a vertically stacked layout. Each track
#' panel is split up into a title section containing the track name, as well as
#' an optional axis for tracks containing numeric data, and a data section
#' showing the actual data along genomic coordinates. In that sense, the output
#' is very similar to the UCSC Genome Browser.
#'
#' The layout of the individual tracks is highly customizable though so called
#' "display parameters". See \code{\link{settings}} for details.
#'
#' While plotting a track, the software automatically computes HTML image map
#' coordinates based on the current graphics device. These coordinates as well
#' as the associated annotation information can later be used to embed images
#' of the plots in semi-interactive HTML pages. See
#' \code{\linkS4class{ImageMap}} for details.
#'
#' @name plotTracks
#'
#' @param trackList A list of Gviz track objects, all inheriting from
#' class \code{\linkS4class{GdObject}}. The tracks will all be drawn to the
#' same genomic coordinates, either as defined by the `from` and `to`
#' arguments if supplied, or by the maximum range across all individual items
#' in the list.
#' @param from,to Character scalar, giving the range of genomic coordinates to
#' draw the tracks in. Note that `to` cannot be larger than `from`.
#' If `NULL`, the plotting ranges are derived from the individual tracks.
#' See `extend.left` and `extend.right` below for the definition of
#' the final plotting ranges.
#' @param \dots Additional arguments which are all interpreted as display
#' parameters to tweak the appearance of the plot. These parameters are global,
#' meaning that they will be used for all tracks in the list where they
#' actually make sense, and they override the track-internal settings. See
#' [settings] for details on display parameters.
#' @param sizes A numeric vector of relative vertical sizes for the individual
#' tracks of length equal to the number of tracks in `trackList`, or
#' `NULL` to auto-detect the most appropriate vertical size proportions.
#' @param panel.only Logical flag, causing the tracks to be plotted as
#' lattice-like panel functions without resetting the plotting canvas and
#' omitting the title pane. This allows to embed tracks into a trellis layout.
#' Usually the function is called for a single track only when
#' `panel.only==TRUE`.
#' @param extend.right,extend.left Numeric scalar, extend the plotting range to
#' the right or to the left by a fixed number of bases. The final plotting
#' range is defined as `from-extend.left` to `to+extend.right`.
#' @param title.width A expansion factor for the width of the title panels.
#' This can be used to make more space, e.g. to accommodate for more detailed
#' data axes. The default is to use as much space as needed to fit all the
#' annotation text.
#' @param add Logical flag, add the plot to an existing plotting canvas without
#' re-initialising.
#' @param main Character scalar, the plots main header.
#' @param cex.main,fontface.main,col.main The fontface, color and expansion
#' factor settings for the main header.
#' @param margin The margin width to add to the plot in pixels.
#' @param innerMargin The inner margin width to add to the plot in pixels.
#' @param chromosome Set the chromosome for all the tracks in the track list.
#'
#' @return A list of Gviz tracks, each one augmented by the computed image map
#' coordinates in the `imageMap` slot, along with the additional `ImageMap`
#' object `titles` containing information about the title panels.
#'
#' @author Florian Hahne
#'
#' @seealso
#' \code{\linkS4class{GdObject}}
#'
#' \code{\linkS4class{ImageMap}}
#'
#' \code{\linkS4class{ImageMap}}
#'
#' \code{\linkS4class{StackedTrack}}
#'
#' \code{\link{settings}}
#'
#' @examples
#' ## Create some tracks to plot
#' st <- c(2000000, 2070000, 2100000, 2160000)
#' ed <- c(2050000, 2130000, 2150000, 2170000)
#' str <- c("-", "+", "-", "-")
#' gr <- c("Group1", "Group2", "Group1", "Group3")
#' annTrack <- AnnotationTrack(
#'     start = st, end = ed, strand = str, chromosome = 7,
#'     genome = "hg19", feature = "test", group = gr,
#'     id = paste("annTrack item", 1:4),
#'     name = "annotation track foo",
#'     stacking = "squish"
#' )
#'
#' ax <- GenomeAxisTrack()
#'
#' dt <- DataTrack(
#'     start = seq(min(st), max(ed), len = 10), width = 18000,
#'     data = matrix(runif(40), nrow = 4), genome = "hg19", chromosome = 7,
#'     type = "histogram", name = "data track bar"
#' )
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(annTrack) <- list(fontfamily = font, fontfamily.title = font)
#'     displayPars(ax) <- list(fontfamily = font, fontfamily.title = font)
#'     displayPars(dt) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#'
#' ## Now plot the tracks
#' res <- plotTracks(list(ax, annTrack, dt))
#'
#' ## Plot only a subrange
#' res <- plotTracks(list(ax, annTrack, dt), from = 2080000, to = 2156000)
#'
#' ## Extend plotting ranges
#' res <- plotTracks(list(ax, annTrack, dt), extend.left = 200000, extend.right = 200000)
#'
#' ## Add a header
#' res <- plotTracks(list(ax, annTrack, dt),
#'     main = "A GenomGraphs plot",
#'     col.main = "darkgray"
#' )
#'
#' ## Change vertical size and title width
#' res <- plotTracks(list(ax, annTrack, dt), sizes = c(1, 1, 5))
#'
#' names(annTrack) <- "foo"
#' res <- plotTracks(list(ax, annTrack), title.width = 0.6)
#'
#' ## Adding and lattice like plots
#' library(grid)
#' grid.newpage()
#' pushViewport(viewport(height = 0.5, y = 1, just = "top"))
#' grid.rect()
#' plotTracks(annTrack, add = TRUE)
#' popViewport(1)
#' pushViewport(viewport(height = 0.5, y = 0, just = "bottom"))
#' grid.rect()
#' plotTracks(dt, add = TRUE)
#' popViewport(1)
#' \dontrun{
#' library(lattice)
#' myPanel <- function(x, ...) {
#'     plotTracks(annTrack,
#'         panel.only = TRUE,
#'         from = min(x), to = max(x), shape = "box"
#'     )
#' }
#' a <- seq(1900000, 2250000, len = 40)
#' xyplot(b ~ a | c, data.frame(a = a, b = 1, c = cut(a, 4)),
#'     panel = myPanel,
#'     scales = list(x = "free")
#' )
#' }
#'
#' @export
#' @importFrom lattice current.panel.limits panel.abline panel.grid panel.lines panel.points panel.polygon panel.segments panel.xyplot panel.text trellis.par.get
#'
plotTracks <- function(trackList, from = NULL, to = NULL, ..., sizes = NULL, panel.only = FALSE, extend.right = 0,
                       extend.left = 0, title.width = NULL, add = FALSE, main, cex.main = 2, fontface.main = 2,
                       col.main = "black", margin = 6, chromosome = NULL, innerMargin = 3) {
    ## If we have to open a new device for this but do not run through the whole function because of errors we want to
    ## clean up in the end
    done <- FALSE
    cdev <- dev.cur()
    on.exit(if (cdev == 1 && !done) dev.off())
    ## We only need a new plot for regular calls to the function. Both add==TRUE and panel.only=TRUE will add to an existing grid plot
    if (!panel.only && !add) {
        grid.newpage()
    }
    if (!is.list(trackList)) {
        trackList <- list(trackList)
    }
    ## All arguments in ... are considered to be additional display parameters and need to be attached to each item in the track list
    dps <- list(...)
    trackList <- lapply(trackList, function(x) {
        displayPars(x, recursive = TRUE) <- dps
        return(x)
    })

    ## OverlayTracks and HighlightTracks can be discarded if they are empty
    trackList <- trackList[!vapply(trackList, function(x) (is(x, "HighlightTrack") || is(x, "OverlayTrack")) && length(x) < 1, FUN.VALUE = logical(1L))]
    isHt <- which(vapply(trackList, is, "HighlightTrack", FUN.VALUE = logical(1L)))
    isOt <- which(vapply(trackList, is, "OverlayTrack", FUN.VALUE = logical(1L)))
    ## A mix between forward and reverse strand tracks should trigger an alarm
    strds <- unique(.whichStrand(trackList))
    if (!is.null(strds) && length(strds) > 1) {
        warning("Plotting a mixture of forward strand and reverse strand tracks.\n Are you sure this is correct?")
    }
    ## We first run very general housekeeping tasks on the tracks for which we don't really need to know anything about device
    ## size, resolution or plotting ranges.
    ## Chromosomes should all be the same for all tracks, if not we will force them to be set to the first one that can be detected.
    ## If plotting ranges are supplied we can speed up a lot of the downstream operations by subsetting first.
    ## We may want to use alpha blending on those devices that support it, but also fall back to non-transparent colors without causing
    ## warnings.
    hasAlpha <- .supportsAlpha()
    chrms <- unique(unlist(lapply(trackList, .recChromosome)))
    if (is.null(chromosome)) {
        chrms <- if (!is.null(chrms)) chrms[gsub("^chr", "", chrms) != "NA"] else chrms
        chromosome <- head(chrms, 1)
        if (length(chromosome) == 0) {
            chromosome <- "chrNA"
        }
        if (!is.null(chrms) && length(unique(chrms)) != 1) {
            warning("The track chromosomes in 'trackList' differ. Setting all tracks to chromosome '", chromosome, "'", sep = "")
        }
    }
    if (!is.null(from) || !(is.null(to))) {
        trackList <- lapply(trackList, function(x) {
            chromosome(x) <- chromosome
            subset(x, from = from, to = to, chromosome = chromosome, sort = FALSE, stacks = FALSE, use.defaults = FALSE)
        })
    }
    trackList <- lapply(trackList, consolidateTrack,
        chromosome = chromosome, any(.needsAxis(trackList)), any(.needsTitle(trackList)),
        title.width, alpha = hasAlpha, ...
    )

    ## Now we figure out the plotting ranges. If no ranges are given as function arguments we take the absolute min/max of all tracks.
    ranges <- .defaultRange(trackList, from = from, to = to, extend.left = extend.left, extend.right = extend.right, annotation = TRUE)
    ## Now we can subset all the objects in the list to the current boundaries and compute the initial stacking
    trackList <- lapply(trackList, subset, from = ranges["from"], to = ranges["to"], chromosome = chromosome)
    trackList <- lapply(trackList, setStacks, recomputeRanges = FALSE)
    ## Highlight tracks are just a way to add a common highlighting region to several tracks, but other than that we can treat the containing
    ## tracks a normal track objects, and thus unlist them. We only want to record their indexes in the expanded list for later.
    htList <- list()
    expandedTrackList <- if (length(isHt)) {
        j <- 1
        tlTemp <- list()
        for (i in seq_along(trackList)) {
            if (!i %in% isHt) {
                tlTemp <- c(tlTemp, trackList[[i]])
                j <- j + 1
            } else {
                tlTemp <- c(tlTemp, trackList[[i]]@trackList)
                htList[[as.character(i)]] <- list(
                    indexes = j:(j + length(trackList[[i]]@trackList) - 1),
                    track = trackList[[i]]
                )
                j <- j + length(trackList[[i]]@trackList)
            }
        }
        tlTemp
    } else {
        trackList
    }
    ## If there is a AlignmentsTrack and also a SequenceTrack we can tell the former to use the latter, unless already provided
    isAt <- vapply(expandedTrackList, is, "AlignmentsTrack", FUN.VALUE = logical(1L))
    isSt <- vapply(expandedTrackList, is, "SequenceTrack", FUN.VALUE = logical(1L))
    for (ai in which(isAt)) {
        if (is.null(expandedTrackList[[ai]]@referenceSequence) && any(isSt)) {
            expandedTrackList[[ai]]@referenceSequence <- expandedTrackList[[min(which(isSt))]]
        }
    }
    ## We need to reverse the list to get a top to bottom plotting order
    expandedTrackList <- rev(expandedTrackList)
    map <- vector(mode = "list", length = length(expandedTrackList))
    titleCoords <- NULL
    names(map) <- rev(vapply(expandedTrackList, names, FUN.VALUE = character(1L)))
    ## Open a fresh page and set up the bounding box, unless add==TRUE
    if (!panel.only) {
        ## We want a margin pixel border
        ## for backward compatibility, if margin has length of 2,
        ## the first one will be used as a horizontal, second as a vertical margin
        if (length(margin) == 2) {
            margin <- rev(margin)
        }
        ## we switched to same settings as in par
        ## c(bottom, left, top, right)
        margin <- rep(as.numeric(margin), length.out = 4)
        vpWidth <- vpLocation()$size["width"]
        vpHeight <- vpLocation()$size["height"]
        vpBound <- viewport(
            x = margin[2L] / vpWidth, y = margin[1L] / vpHeight,
            width = (vpWidth - sum(margin[c(2, 4)])) / vpWidth,
            height = (vpHeight - sum(margin[c(1, 3)])) / vpHeight,
            just = c("left", "bottom")
        )
        pushViewport(vpBound)
        ## If there is a header we have to make some room for it here
        if (!missing(main) && main != "") {
            vpHeader <- viewport(width = 1, height = 0.1, y = 1, just = c("center", "top"))
            pushViewport(vpHeader)
            grid.text(main, gp = gpar(col = col.main, cex = cex.main, fontface = fontface.main))
            popViewport(1)
            vpMain <- viewport(width = 1, height = 0.9, y = 0.9, just = c("center", "top"))
        } else {
            vpMain <- viewport(width = 1, height = 1)
        }
        pushViewport(vpMain)
        ## A first guestimate of the vertical space that's needed
        spaceSetup <- .setupTextSize(expandedTrackList, sizes, title.width, spacing = innerMargin)
    } else {
        vpBound <- viewport()
        pushViewport(vpBound)
        spaceSetup <- .setupTextSize(expandedTrackList, sizes, spacing = innerMargin)
    }
    ## First iteration to set up all the dimensions by calling the drawGD methods in prepare mode, i.e.,
    ## argument prepare=TRUE. Nothing is drawn at this point, and this only exists to circumvent the
    ## chicken and egg problem of not knowing how much space we need until we draw, but also not knowing
    ## where to draw until we know the space needed.
    for (i in rev(seq_along(expandedTrackList)))
    {
        fontSettings <- .fontGp(expandedTrackList[[i]], cex = NULL)
        vpTrack <- viewport(
            x = 0, y = sum(spaceSetup$spaceNeeded[seq_len(i)]), just = c(0, 1), width = 1, height = spaceSetup$spaceNeeded[i],
            gp = fontSettings
        )
        pushViewport(vpTrack)
        vpContent <- if (!panel.only) {
            viewport(
                x = spaceSetup$title.width + spaceSetup$spacing,
                width = 1 - spaceSetup$title.width - spaceSetup$spacing * 2, just = 0
            )
        } else {
            viewport(width = 1)
        }
        pushViewport(vpContent)
        expandedTrackList[[i]] <- drawGD(expandedTrackList[[i]], minBase = ranges["from"], maxBase = ranges["to"], prepare = TRUE, subset = FALSE)
        popViewport(2)
    }
    ## Now lets recalculate the space and draw for real
    spaceSetup <- .setupTextSize(expandedTrackList, sizes, title.width, spacing = innerMargin)
    ## First the highlight box backgrounds
    htBoxes <- data.frame(stringsAsFactors = FALSE)
    for (hlite in htList) {
        if (length(ranges(hlite$track))) {
            inds <- setdiff(sort(length(expandedTrackList) - hlite$index + 1), which(vapply(expandedTrackList, is, "IdeogramTrack", FUN.VALUE = logical(1L))))
            y <- reduce(IRanges(start = inds, width = 1))
            yy <- ifelse(start(y) == 1, 0, sum(spaceSetup$spaceNeeded[seq_len(start(y)) - 1])) # check
            ht <- sum(spaceSetup$spaceNeeded[start(y):end(y)])
            htBoxes <- rbind(htBoxes, data.frame(
                y = yy, height = ht, x = start(hlite$track), width = width(hlite$track),
                col = .dpOrDefault(hlite$track, "col", "orange"),
                fill = .dpOrDefault(hlite$track, "fill", "red"),
                lwd = .dpOrDefault(hlite$track, "lwd", 1),
                lty = .dpOrDefault(hlite$track, "lty", 1),
                alpha = .dpOrDefault(hlite$track, "alpha", 1),
                inBackground = .dpOrDefault(hlite$track, "inBackground", TRUE),
                stringsAsFactors = FALSE
            ))
        }
    }
    .drawHtBoxes <- function(htBoxes, background = TRUE) {
        htBoxes <- htBoxes[htBoxes$inBackground == background, , drop = FALSE]
        rscales <- if (strds[1] == "reverse") c(from = ranges["to"], to = ranges["from"]) else ranges
        if (nrow(htBoxes)) {
            vpContent <- if (!panel.only) {
                viewport(
                    x = spaceSetup$title.width + spaceSetup$spacing, xscale = rscales,
                    width = 1 - spaceSetup$title.width - spaceSetup$spacing * 2, just = 0
                )
            } else {
                viewport(width = 1, xscale = rscales)
            }
            pushViewport(vpContent)
            grid.rect(
                x = htBoxes$x, just = c(0, 1), width = htBoxes$width, y = htBoxes$y + htBoxes$height, height = htBoxes$height,
                gp = gpar(col = htBoxes$col, fill = htBoxes$fill, lwd = htBoxes$lwd, lty = htBoxes$lty, alpha = unique(htBoxes$alpha)), default.units = "native"
            )
            popViewport(1)
        }
    }
    if (nrow(htBoxes)) {
        .drawHtBoxes(htBoxes)
    }
    ## Now the track content
    for (i in rev(seq_along(expandedTrackList))) {
        vpTrack <- viewport(x = 0, y = sum(spaceSetup$spaceNeeded[seq_len(i)]), just = c(0, 1), width = 1, height = spaceSetup$spaceNeeded[i])
        pushViewport(vpTrack)
        fill <- .dpOrDefault(expandedTrackList[[i]], "background.title", .DEFAULT_SHADED_COL)
        thisTrack <- if (is(expandedTrackList[[i]], "OverlayTrack")) {
            tmpThisTrack <- expandedTrackList[[i]]@trackList
             while (any(vapply(tmpThisTrack, is, "OverlayTrack", FUN.VALUE = logical(1L)))) {
                 tmpThisTrack <- rapply(tmpThisTrack, function(x) if (is(x, "OverlayTrack")) x@trackList else x)
             }
            tmpSpaceSetup <- .setupTextSize(list(tmpThisTrack[[1]]), sizes[1], title.width, spacing = innerMargin)
            spaceSetup$nwrap[i] <- tmpSpaceSetup$nwrap[1]
            if (spaceSetup$title.width < tmpSpaceSetup$title.width) {
                spaceSetup$title.width <- tmpSpaceSetup$title.width
            }
            tmpThisTrack[[1]]
        } else {
            expandedTrackList[[i]]
        }
        if (!panel.only) {
            fontSettings <- .fontGp(thisTrack, subtype = "title", cex = NULL)
            vpTitle <- viewport(x = 0, width = spaceSetup$title.width, just = 0, gp = fontSettings)
            pushViewport(vpTitle)
            lwd.border.title <- .dpOrDefault(thisTrack, "lwd.title", 1)
            col.border.title <- .dpOrDefault(thisTrack, "col.border.title", "transparent")
            grid.rect(gp = gpar(fill = fill, col = col.border.title, lwd = lwd.border.title))
            needAxis <- .needsAxis(thisTrack)
            drawAxis(thisTrack, ranges["from"], ranges["to"], subset = FALSE)
            tit <- spaceSetup$nwrap[i]
            ## FIXME: Do we want something smarted for the image map coordinates?
            titleCoords <- rbind(titleCoords, cbind(.getImageMap(cbind(0, 0, 1, 1)),
                title = names(thisTrack)
            ))
            if (.dpOrDefault(thisTrack, "showTitle", TRUE) && !is.null(tit) && tit != "") {
                x <- if (needAxis) 0.075 else 0.4
                just <- if (needAxis) c("center", "top") else "center"
                ## FIXME: We need to deal with this when calculating the space for the title bar
                rot <- .dpOrDefault(thisTrack, "rotation.title", 90)
                gp <- .fontGp(thisTrack, "title", cex = spaceSetup$cex[i])
                suppressWarnings(grid.text(tit, unit(x, "npc"), rot = rot, gp = gp, just = just))
            }
            popViewport(1)
        }
        ## Draw the panel background, grid lines if necessary and the panel content
        vpBackground <- if (!panel.only) {
            viewport(
                x = spaceSetup$title.width,
                width = 1 - spaceSetup$title.width, just = 0
            )
        } else {
            viewport(width = 1)
        }
        pushViewport(vpBackground)
        grid.rect(gp = gpar(col = "transparent", fill = .dpOrDefault(thisTrack, "background.panel", "transparent")))
        drawGrid(thisTrack, ranges["from"], ranges["to"])
        popViewport(1)
        fontSettings <- .fontGp(expandedTrackList[[i]], cex = NULL)
        vpContentOuter <- if (!panel.only) {
            viewport(
                x = spaceSetup$title.width, width = 1 - spaceSetup$title.width,
                just = 0, gp = fontSettings, clip = TRUE
            )
        } else {
            viewport(width = 1, gp = fontSettings, clip = TRUE)
        }
        pushViewport(vpContentOuter)
        vpContent <- if (!panel.only) {
            viewport(x = spaceSetup$spacing, width = 1 - (spaceSetup$spacing * 2), just = 0, gp = fontSettings)
        } else {
            viewport(width = 1, gp = fontSettings)
        }
        pushViewport(vpContent)
        tmp <- drawGD(expandedTrackList[[i]], minBase = ranges["from"], maxBase = ranges["to"], subset = FALSE)
        if (!is.null(tmp)) {
            map[[(length(map) + 1) - i]] <- tmp
        }
        popViewport(2)
        if (.dpOrDefault(thisTrack, "frame", FALSE)) {
            grid.rect(gp = gpar(col = .dpOrDefault(thisTrack, "col.frame", .DEFAULT_SHADED_COL), fill = "transparent"))
        }
        popViewport(1)
    }
    if (nrow(htBoxes)) {
        .drawHtBoxes(htBoxes, FALSE)
    }
    popViewport(if (panel.only) 1 else 2)
    tc <- as.character(titleCoords[, 5])
    tc[which(tc == "" | is.na(tc) | is.null(tc))] <- "NA"
    names(tc) <- tc
    if (!is.null(titleCoords)) {
        tcoord <- as.matrix(titleCoords[, seq(1, 4)])
        rownames(tcoord) <- names(tc)
        map$titles <- ImageMap(coords = tcoord, tags = list(title = tc))
    }
    done <- TRUE
    return(invisible(map))
}
