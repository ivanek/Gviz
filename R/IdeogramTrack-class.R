#' @include GdObject-class.R
NULL

## IdeogramTrack Class -------------------------------------------------------

#' IdeogramTrack class and methods
#'
#'
#' A class to represent the schematic display of a chromosome, also known as an
#' ideogram. The respective information is typically directly fetched from
#' UCSC.
#'
#'
#' Ideograms are schematic depictions of chromosomes, including chromosome band
#' information and centromere location. The relevant data for various species
#' is stored in the UCSC data base. The initializer method of the class will
#' automatically fetch the respective data for a given genome and chromosome
#' from UCSC and fill the appropriate object slots. When plotting
#' \code{IdeogramTracks}, the current genomic location is indicated on the
#' chromosome by a colored box.
#'
#' The \code{Gviz.ucscUrl} option controls which URL is being used to connect
#' to UCSC. For instance, one could switch to the European UCSC mirror by
#' calling \code{options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/"}.
#'
#' @name IdeogramTrack-class
#' @aliases IdeogramTrack-class IdeogramTrack drawGD,IdeogramTrack-method
#' end,IdeogramTrack-method end<-,IdeogramTrack-method
#' initialize,IdeogramTrack-method show,IdeogramTrack-method
#' start,IdeogramTrack-method start<-,IdeogramTrack-method
#' subset,IdeogramTrack-method width,IdeogramTrack-method
#' width<-,IdeogramTrack-method length,IdeogramTrack-method
#' [,IdeogramTrack-method [,IdeogramTrack,ANY,ANY-method
#' [,IdeogramTrack,ANY,ANY,ANY-method chromosome<-,IdeogramTrack-method
#' genome<-,IdeogramTrack-method position,IdeogramTrack-method
#' @docType class
#' @param chromosome The chromosome for which to create the ideogram. Has to be
#' a valid UCSC chromosome identifier of the form \code{chrx}, or a single
#' integer or numeric character unless
#' \code{option(ucscChromosomeNames=FALSE)}. The user has to make sure that the
#' respective chromosome is indeed defined for the the track's genome.
#' @param genome The genome on which to create the ideogram. This has to be a
#' valid UCSC genome identifier if the ideogram data is to be fetched from the
#' UCSC repository.
#' @param name Character scalar of the track's name used in the title panel
#' when plotting. Defaults to the selected chromosome.
#' @param bands A \code{data.frame} with the cytoband information for all
#' available chromosomes on the genome similar to the data that would be
#' fetched from UCSC. The table needs to contain the mandatory columns
#' \code{chrom}, \code{chromStart}, \code{chromEnd}, \code{name} and
#' \code{gieStain} with the chromosome name, cytoband start and end
#' coordinates, cytoband name and coloring information, respectively. This can
#' be used when no connection to the internet is available or when the cytoband
#' information has been cached locally to avoid the somewhat slow connection to
#' UCSC.
#' @param \dots Additional items which will all be interpreted as further
#' display parameters.
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{IdeogramTrack}.
#' @note
#'
#' When fetching ideogram data from UCSC the results are cached for faster
#' acces. See \code{\link{clearSessionCache}} on details to delete these cached
#' items.
#' @section Objects from the Class:
#'
#' Objects can be created using the constructor function \code{IdeogramTrack}.
#' @author Florian Hahne
#' @inherit GdObject-class seealso
#'
#' @examples
#' \dontshow{
#' ## Load some sample data
#' data(idTrack)
#' }
#'
#' ## Construct the object
#' \dontrun{
#' idTrack <- IdeogramTrack(chromosome = 7, genome = "mm9")
#' }
#'
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(idTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#' ## Plotting
#' plotTracks(idTrack, from = 5000000, to = 9000000)
#'
#' ## Track names
#' names(idTrack)
#' names(idTrack) <- "foo"
#' plotTracks(idTrack, from = 5000000, to = 9000000)
#'
#'
#' ## Accessors
#' chromosome(idTrack)
#' \dontrun{
#' chromosome(idTrack) <- "chrX"
#' }
#'
#' genome(idTrack)
#' \dontrun{
#' genome(id) <- "hg19"
#' }
#'
#' range(idTrack)
#' ranges(idTrack)
#'
#' ## Annotation
#' values(idTrack)
#'
#' ## coercion
#' as(idTrack, "data.frame")
#' @exportClass IdeogramTrack
setClass("IdeogramTrack",
    contains = "RangeTrack",
    representation = representation(bandTable = "data.frame"),
    prototype = prototype(
        name = "IdeogramTrack",
        bandTable = data.frame(),
        dp = DisplayPars(
            background.title = "transparent",
            bevel = 0.45,
            centromereShape = "triangle",
            cex.bands = 0.7,
            cex = 0.8,
            col = "red",
            col.border.title = "transparent",
            lwd.border.title = 1,
            fill = "#FFE3E6",
            fontface = 1,
            fontfamily = "sans",
            fontcolor = .DEFAULT_SHADED_COL,
            fontsize = 10,
            outline = FALSE,
            showBandId = FALSE,
            lty = 1,
            lwd = 1,
            showId = TRUE,
            showTitle = FALSE,
            size = NULL
        )
    )
)

## Initialize ----------------------------------------------------------------

## Grab the chromosome band and length information from UCSC and fill the ranges slot.

#' @export
setMethod("initialize", "IdeogramTrack", function(.Object, genome, chromosome, bands, name, ...) {
    ## the display parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "IdeogramTrack")
    if (missing(bands)) {
        bands <- NULL
    }
    if (is.null(bands) && (missing(genome) || missing(chromosome))) {
        return(callNextMethod(.Object = .Object, range = GRanges(), genome = NULL, chromosome = NULL, ...))
    }
    if (is.null(bands)) {
        sessionInfo <- .cacheGenomes(genome = genome)
        .Object@bandTable <- sessionInfo$bands
        bands <- sessionInfo$bands
    } else {
        .checkClass(bands, "data.frame")
        cols <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
        miss <- !cols %in% colnames(bands)
        if (any(miss)) {
            stop(sprintf(
                "The following column%s missing from the bands table: %s",
                ifelse(sum(miss) > 1, "s are", " is"), paste(cols[miss], collapse = ", ")
            ))
        }
        .Object@bandTable <- bands
    }
    chromosome <- if (is.null(chromosome)) as.character(bands[1, "chrom"]) else .chrName(chromosome)[1]
    bands <- bands[bands$chrom == chromosome, ]
    if (nrow(bands) == 0) {
        stop("Chromosome '", chromosome, "' does not exist on UCSC genome '", genome, "'")
    }
    if (is.null(name)) {
        name <- .chrName(chromosome)[1]
    }
    bnames <- as.character(bands$name)
    sel <- is.na(bnames)
    if (any(sel)) {
        bnames[sel] <- paste("band", seq_len(sum(sel)), sep = "_")
    }
    if (any(bnames == "")) {
        bnames[bnames == ""] <- sprintf("band_%i", which(bnames == ""))
    }
    ranges <- GRanges(
        seqnames = bnames, ranges = IRanges(start = bands$chromStart, end = bands$chromEnd),
        name = bnames, type = as.character(bands$gieStain)
    )
    .Object <- callNextMethod(.Object = .Object, range = ranges, genome = genome, chromosome = chromosome, name = name, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters
#' @export
IdeogramTrack <- function(chromosome = NULL, genome, name = NULL, bands = NULL, ...) {
    if (missing(genome)) stop("Need to specify genome for creating an IdeogramTrack")
    new("IdeogramTrack", chromosome = chromosome, genome = genome, name = name, bands = bands, ...)
}

## General accessors ---------------------------------------------------------

#' @export
setMethod("start", "IdeogramTrack", function(x) NULL) # if(length(x)) start(range(x)) else NULL)

#' @export
setReplaceMethod("start", "IdeogramTrack", function(x, value) {
    # start(x@range) <- value
    return(x)
})

#' @export
setMethod("end", "IdeogramTrack", function(x) NULL) # if(length(x)) end(range(x)) else NULL)

#' @export
setReplaceMethod("end", "IdeogramTrack", function(x, value) {
    # end(x@range) <- value
    return(x)
})

#' @export
setMethod("width", "IdeogramTrack", function(x) NULL) # if(length(x)) width(range(x)) else NULL)

#' @export
setReplaceMethod("width", "IdeogramTrack", function(x, value) {
    return(x)
})

#' @export
setMethod("length", "IdeogramTrack", function(x) length(ranges(x)))

#' @export
setReplaceMethod("chromosome", "IdeogramTrack", function(GdObject, value) {
    ## We have changed the class definition to include the bands for all chromosomes, but still want the old objects to work
    chromosome <- .chrName(value[1])
    if (.hasSlot(GdObject, "bandTable") && chromosome %in% as.character(GdObject@bandTable$chrom)) {
        ranges <- GdObject@bandTable[GdObject@bandTable$chrom == chromosome, ]
        bnames <- as.character(ranges$name)
        sel <- is.na(bnames)
        if (any(sel)) {
            bnames[sel] <- paste("band", seq_len(sum(sel)), sep = "_")
        }
        if (any(bnames == "")) {
            bnames[bnames == ""] <- sprintf("band_%i", which(bnames == ""))
        }
        ranges <- GRanges(
            seqnames = bnames, ranges = IRanges(start = ranges$chromStart, end = ranges$chromEnd),
            name = bnames, type = ranges$gieStain
        )
        GdObject@range <- ranges
        GdObject@chromosome <- chromosome
        return(GdObject)
    }
    message("Updating chromosome band information")
    tmp <- IdeogramTrack(genome = genome(GdObject), chromosome = .chrName(value[1]), name = names(GdObject))
    displayPars(tmp) <- displayPars(GdObject, hideInternal = FALSE)
    return(tmp)
})

#' @export
setReplaceMethod("genome", "IdeogramTrack", function(x, value) {
    if (genome(x) != value) {
        message("Updating chromosome band information")
    }
    tmp <- IdeogramTrack(genome = value[1], chromosome = chromosome(x), name = names(x))
    displayPars(tmp) <- displayPars(x, hideInternal = FALSE)
    return(tmp)
})


## Annotation Accessors ------------------------------------------------------
## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------

#' @export
setMethod("[", signature(x = "IdeogramTrack"), function(x, i, j, ..., drop = TRUE) {
    return(x)
})

## Position ------------------------------------------------------------------

#' @export
setMethod("position", signature("IdeogramTrack"), definition = function(GdObject, ...) NULL)

## DrawGD --------------------------------------------------------------------

## Helper function to compute coordinates for a rounded ideogram cap
.roundedCap <- function(bl, tr, st, vals, side = c("left", "right"), bevel = 0.4, n = 100) {
    side <- match.arg(side)
    bevel <- max(1 / n, min(bevel, 0.5))
    coords <- if (bevel <= 0) {
        cbind(c(0, 1, 1, 0), c(1.5, 1.5, -1.5, -1.5))
    } else {
        cbind(
            c(0, sin(seq(0, pi / 2, len = max(2, n * bevel))) + bevel * 4, sin(seq(pi / 2, pi, len = max(2, n * bevel))) + bevel * 4, 0) / (1 + 4 * bevel),
            (c(
                1 + 0.5 - bevel, cos(seq(0, pi / 2, len = max(2, n * bevel))) + (0.5 - bevel),
                cos(seq(pi / 2, pi, len = max(2, n * bevel))) - (0.5 - bevel), cos(pi) - (0.5 - bevel)
            ) + 1 + (0.5 - bevel)) / (2 + 2 * (0.5 - bevel))
        )
    }
    if (side == "right") {
        coords[, 1] <- coords[, 1] * abs(diff(c(bl[1], tr[1]))) + bl[1]
        coords[, 2] <- coords[, 2] * abs(diff(c(bl[2], tr[2]))) + bl[2]
    } else {
        coords[, 1] <- (1 - coords[, 1]) * abs(diff(c(bl[1], tr[1]))) + bl[1]
        coords[, 2] <- (1 - coords[, 2]) * abs(diff(c(bl[2], tr[2]))) + bl[2]
    }
    lcS <- split(as.data.frame(coords), cut(coords[, 1], st, right = TRUE, include.lowest = TRUE, labels = FALSE), drop = TRUE)
    first <- TRUE
    shift <- ifelse(side == "left", 0, 1)
    for (j in names(lcS)) {
        xx <- lcS[[j]][, 1]
        yy <- lcS[[j]][, 2]
        if (!first) {
            prev <- lcS[[as.character(as.numeric(j) - 1)]]
            xx <- c(tail(prev[, 1], 1), xx, tail(prev[, 1], 1))
            yy <- c(1 - tail(prev[, 2], 1), yy, tail(prev[, 2], 1))
        }
        grid.polygon(xx, yy, gp = gpar(col = vals[as.numeric(j) + shift, "col"], fill = vals[as.numeric(j) + shift, "col"]))
        first <- FALSE
    }
    return(coords)
}

## A more generic method to come up with colors for chromosome bands that still relies a bit on biovizBase
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
.getBioColorIdeo <- function(type) {
    type <- as.character(type)
    ocols <- getBioColor("CYTOBAND")
    cols <- c(ocols[c("gneg", "stalk", "acen")], gpos = unname(ocols["gpos100"]), gvar = unname(ocols["gpos100"]))
    gpcols <- unique(grep("gpos", type, value = TRUE))
    crmp <- colorRampPalette(c(cols["gneg"], cols["gpos"]))(100)
    posCols <- setNames(crmp[as.integer(gsub("gpos", "", gpcols))], gpcols)
    return(c(cols, posCols))
}

## The actual drawing method
#' @importFrom grDevices rgb2hsv
#' @export
setMethod("drawGD", signature("IdeogramTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if ((is.logical(debug) && debug) || debug == "prepare") {
        browser()
    }
    imageMap(GdObject) <- NULL
    if (names(GdObject)[1] != chromosome(GdObject)) {
        chrnam <- names(GdObject)[1]
    } else {
        chrnam <- paste("Chromosome", gsub("chr", "", chromosome(GdObject)))
    }
    cex <- .dpOrDefault(GdObject, "cex", 1)
    ## Nothing to do if there are no ranges in the object
    if (!length(GdObject)) {
        return(invisible(GdObject))
    }
    ## In prepare mode we just want to figure out the optimal size
    if (prepare) {
        pres <- .pxResolution()
        nsp <- if (.dpOrDefault(GdObject, "showId", TRUE)) {
            (as.numeric(convertHeight(stringHeight(chrnam), "native")) * cex) +
                as.numeric(convertHeight(unit(20, "points"), "native"))
        } else {
            as.numeric(convertHeight(unit(25, "points"), "native"))
        }
        nsp <- nsp / pres["y"]
        displayPars(GdObject) <- list("neededVerticalSpace" = nsp)
        ## Augment the ranges to fill gaps if there are any
        gaps <- setdiff(IRanges(0, max(end(range(GdObject)))), range(GdObject))
        if (length(gaps)) {
            gaps <- GRanges(seqnames(GdObject)[1], gaps, name = rep(as.character(NA), length(gaps)), type = rep("gneg", length(gaps)))
            rr <- c(ranges(GdObject), gaps)
            ranges(GdObject) <- sort(rr)
        }
        return(invisible(GdObject))
    }
    if ((is.logical(debug) && debug) || debug == "draw") {
        browser()
    }
    ## Do we need some space for the chromosome name?
    if (.dpOrDefault(GdObject, "showId", TRUE)) {
        gp <- .fontGp(GdObject)
        width <- vpLocation()$isize["width"]
        width <- as.numeric(convertWidth(stringWidth(chrnam), "inches")) * cex + 0.2
        wfac <- vpLocation()$isize["width"]
        nspace <- min(width / wfac, 0.75)
        if ((width / wfac) > 0.75) {
            cex <- cex * (0.75 / (width / wfac))
        }
        pushViewport(viewport(x = 0, width = nspace, just = 0, gp = gp))
        grid.text(chrnam, 0, default.units = "native", just = c("left", "center"))
        popViewport(1)
    } else {
        nspace <- 0
    }
    pushViewport(viewport(x = nspace, width = 1 - nspace, just = 0))
    ## A box indicating the current range on the chromosome
    len <- end(range(range(GdObject)))
    fill <- .dpOrDefault(GdObject, "fill", "#FFE3E6")
    if (!missing(minBase) && !missing(maxBase)) {
        grid.rect(minBase / len, 0.1,
            width = min(1, (maxBase - minBase) / len), height = 0.8, just = c("left", "bottom"),
            gp = gpar(col = "transparent", fill = fill)
        )
    }
    ## Color mapping for the bands taken from the biovizBase package
    cols <- .getBioColorIdeo(values(GdObject)$type)
    vals <- data.frame(values(GdObject), col = cols[as.character(values(GdObject)$type)], stringsAsFactors = FALSE)
    ## For the rounded caps we need  to figure out the overlap with existing bands for proper coloring
    bevel <- 0.02
    ol <- queryHits(findOverlaps(range(GdObject), IRanges(start = c(bevel, 1 - bevel) * len, width = 1)))
    st <- start(range(GdObject)) / len
    ed <- end(range(GdObject)) / len
    stExt <- if (length(GdObject) == 1) c(0, bevel, 1 - bevel) else c(st[seq_len(ol[1])], bevel, st[(ol[1] + 1):ol[2]], 1 - bevel)
    valsExt <- if (length(GdObject) == 1) vals[rep(1, 3), ] else rbind(vals[seq_len(ol[1]), ], vals[ol[1], ], vals[(ol[1] + 1):ol[2], ], vals[ol[2], ])
    if (ol[2] < length(st)) {
        stExt <- c(stExt, st[(ol[2] + 1):length(st)])
        valsExt <- rbind(valsExt, vals[(ol[2] + 1):length(st), ])
    }
    wd <- diff(c(stExt, 1))
    ls <- ol[1] + 1
    rs <- ol[2] + 1
    ## The centromere is treated separately
    cent <- grep("acen", valsExt$type)
    if (length(cent)) {
        bef <- ls:(min(cent) - 1)
        aft <- (max(cent) + 1):rs
    } else {
        bef <- ls:rs
        aft <- NULL
    }
    margin <- 0.3
    ## First the normal bands
    lcol <- "black"
    lwd <- 1
    lty <- 1
    lcolBands <- if (.dpOrDefault(GdObject, "outline", FALSE)[1] == TRUE) rep(lcol, nrow(valsExt)) else valsExt[, "col"]
    grid.rect(stExt[bef], margin,
        width = wd[bef], height = 1 - margin * 2, gp = gpar(col = lcolBands[bef], fill = valsExt[bef, "col"]),
        just = c("left", "bottom")
    )
    if (!is.null(aft)) {
        grid.rect(stExt[aft], margin,
            width = wd[aft], height = 1 - margin * 2, gp = gpar(col = lcolBands[aft], fill = valsExt[aft, "col"]),
            just = c("left", "bottom")
        )
    }
    ## Now the centromere, if there is any
    centromereShape <- .dpOrDefault(GdObject, "centromereShape", "triangle")
    if (length(cent) && centromereShape == "triangle") {
        grid.polygon(c(stExt[min(cent)], stExt[min(cent)] + wd[min(cent)], rep(stExt[min(cent)], 2)),
            c(margin, 0.5, (1 - margin), margin),
            gp = gpar(col = cols["acen"], fill = cols["acen"])
        )
        grid.polygon(c(stExt[max(cent)], rep(stExt[max(cent)] + wd[max(cent)], 2), stExt[max(cent)]),
            c(0.5, margin, (1 - margin), 0.5),
            gp = gpar(col = cols["acen"], fill = cols["acen"])
        )
    }
    ## Now the caps
    str <- if (length(st) == 1) c(0, 1) else st
    edr <- if (length(ed) == 1) c(1, 2) else ed
    lc <- .roundedCap(c(stExt[1], margin), c(stExt[ls], 1 - margin), str, vals, side = "left", bevel = .dpOrDefault(GdObject, "bevel", 0.45))
    rc <- .roundedCap(c(tail(stExt, 1), margin), c(1, 1 - margin), edr, vals, side = "right", bevel = .dpOrDefault(GdObject, "bevel", 0.45))
    ## Now some outlines
    grid.lines(lc[, 1], lc[, 2], gp = gpar(col = lcol, lwd = lwd, lty = lty))
    grid.lines(rc[, 1], rc[, 2], gp = gpar(col = lcol, lwd = lwd, lty = lty))
    if (length(cent)) {
        x0 <- c(rep(max(lc[, 1]), 2), rep(stExt[max(cent) + 1], 2), rep(stExt[min(cent)], 2), rep(stExt[max(cent)], 2))
        y0 <- c(rep(c(margin, (1 - margin)), 3), 0.5, 0.5)
        x1 <- c(
            rep(stExt[min(cent)], 2), rep(tail(stExt, 1), 2), rep(stExt[min(cent)] + wd[min(cent)], 2),
            rep(stExt[max(cent)] + wd[max(cent)], 2)
        )
        y1 <- c(rep(c(margin, (1 - margin)), 2), 0.5, 0.5, margin, (1 - margin))
    } else {
        x0 <- rep(max(lc[, 1]), 2)
        y0 <- c(margin, (1 - margin))
        x1 <- rep(max(stExt), 2)
        y1 <- y0
    }
    grid.segments(x0, y0, x1, y1, gp = gpar(col = lcol, lwd = lwd, lty = lty))
    ## centromere (circle)
    if (length(cent) && centromereShape == "circle") {
        sc <- as.numeric(convertHeight(unit(sum(wd[cent]), "points"), "native")) /
            as.numeric(convertWidth(unit(sum(wd[cent]), "points"), "native"))
        grid.circle(
            x = stExt[min(cent)] + sum(wd[cent]) / 2, y = 0.5, r = sc * sum(wd[cent]) / 2,
            gp = gpar(col = lcol, fill = cols["acen"])
        )
    }
    ## The outlines of the box
    if (!missing(minBase) && !missing(maxBase)) {
        col <- .dpOrDefault(GdObject, "col", "red")
        lwd <- .dpOrDefault(GdObject, "lwd", 1)
        lty <- .dpOrDefault(GdObject, "lty", "solid")
        grid.rect(minBase / len, 0.1,
            width = min(1, (maxBase - minBase) / len), height = 0.8, just = c("left", "bottom"),
            gp = gpar(col = col, fill = "transparent", lwd = lwd, lty = lty)
        )
    }
    ## Finally the band annotation if we need it
    if (.dpOrDefault(GdObject, "showBandId", FALSE)) {
        bn <- as.character(values(GdObject)$name)
        cval <- rgb2hsv(col2rgb(cols[as.character(values(GdObject)$type)]))["v", ]
        tcol <- ifelse(cval > 0.9, "black", "white")
        bwidth <- (c(st[-1], 1) - st) / 2
        cex.bands <- .dpOrDefault(GdObject, "cex.bands", 0.7)
        sspace <- as.numeric(convertUnit(unit(0.01, "inches"), "native"))
        swidth <- as.numeric(convertWidth(stringWidth(bn), "native")) * cex.bands + sspace
        sel <- swidth < bwidth
        if (any(sel)) {
            grid.text(x = (st + bwidth)[sel], y = 0.5, label = bn[sel], hjust = 0.5, gp = gpar(col = tcol[sel], cex = cex.bands))
        }
    }
    popViewport(1)
    return(invisible(GdObject))
})

## Show ----------------------------------------------------------------------

#' @export
setMethod(
    "show", signature(object = "IdeogramTrack"),
    function(object) {
        cat(sprintf(
            paste("Ideogram track '%s' for chromosome %s of the %s genome"),
            names(object),
            gsub("^chr", "", chromosome(object)),
            genome(object)
        ), "\n")
    }
)
