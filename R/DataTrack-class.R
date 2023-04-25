#' @include NumericTrack-class.R
#' @include ReferenceTrack-class.R

NULL

## DataTrack Class -----------------------------------------------------------

#' DataTrack class and methods
#'
#' A class to store numeric data values along genomic coordinates. Multiple
#' samples as well as sample groupings are supported, with the restriction of
#' equal genomic coordinates for a single observation across samples.
#'
#' Depending on the setting of the \code{type} display parameter, the data can
#' be plotted in various different forms as well as combinations thereof.
#' Supported plotting types are:
#'
#' \describe{
#'
#' \item{\code{p}:}{ simple xy-plot.}
#'
#' \item{\code{l}:}{ lines plot. In the case of multiple samples this plotting
#' type is not overly usefull since the points in the data matrix are connected
#' in column-wise order. Type \code{a} might be more appropriate in these
#' situations.}
#'
#' \item{\code{b}:}{ combination of xy-plot and lines plot.}
#'
#' \item{\code{a}:}{ lines plot of the column-wise average values.}
#'
#' \item{\code{s}:}{ sort and connect data points along the x-axis}
#'
#' \item{\code{S}:}{ sort and connect data points along the y-axis}
#'
#' \item{\code{g}:}{ add grid lines. To ensure a consitant look and feel across
#' multiple tracks, grid lines should preferentially be added by using the
#' \code{grid} display parameter.}
#'
#' \item{\code{r}:}{ add a regression line to the plot.}
#'
#' \item{\code{h}:}{ histogram-like vertical lines centered in the middle of
#' the coordinate ranges.}
#'
#' \item{\code{smooth}:}{ add a loess fit to the plot. The following display
#' parameters can be used to control the loess calculation: \code{span, degree,
#' family, evaluation}. See \code{\link{panel.loess}} for details.}
#'
#' \item{\code{histogram}:}{ plot data as a histogram, where the width of the
#' histogram bars reflects the width of the genomic ranges in the \code{range}
#' slot.}
#'
#' \item{\code{mountain}:}{ plot a smoothed version of the data relative to a
#' baseline, as defined by the \code{baseline} display parameter. The following
#' display parameters can be used to control the smoothing: \code{span, degree,
#' family, evaluation}. See \code{\link{panel.loess}} for details. The layout
#' of the plot can be further customized via the following display parameters:
#' \code{col.mountain, lwd.mountain, lty.mountain, fill.mountain}.}
#'
#' \item{\code{polygon}:}{ plot data as a polygon (similar to
#' \code{mountain}-type but without smoothing). Data are plotted relative to a
#' baseline, as defined by the \code{baseline} display parameter. The layout of
#' the plot can be further customized via the following display parameters:
#' \code{col.mountain, lwd.mountain, lty.mountain, fill.mountain}.}
#'
#' \item{\code{boxplot}:}{ plot the data as box-and-whisker plots. The layout
#' of the plot can be further customized via the following display parameters:
#' \code{box.ratio, box.width, varwidt, notch, notch.frac, levels.fos, stats,
#' coef, do.out}. See \code{\link{panel.bwplot}} for details.}
#'
#' \item{\code{gradient}:}{ collapse the data across samples and plot this
#' average value as a color-coded gradient. Essenitally this is similar to the
#' heatmap-type plot of a single sample. The layout of the plot can be further
#' customized via the display parameters \code{ncolor} and \code{gradient}
#' which control the number of gradient colors as well as the gradient base
#' colors, respectively.}
#'
#' \item{\code{heatmap}:}{ plot the color-coded values for all samples in the
#' form of a heatmap. The data for individual samples can be visually separated
#' by setting the \code{separator} display parameter. It's value is taken as
#' the amount of spacing in pixels in between two heatmap rows. The layout of
#' the plot can be further customized via the display parameters \code{ncolor}
#' and \code{gradient} which control the number of gradient colors as well as
#' the gradient base colors, respectively.}
#'
#' \item{\code{horizon}:}{ plot continuous data by cutting the y range into
#' segments and overplotting them with color representing the magnitude and
#' direction of deviation. This is particularly useful when comparing multiple
#' samples, in which case the horizon strips are stacked. See
#' \code{\link{horizonplot}} for details. Please note that the \code{origin}
#' and \code{horizonscale} arguments of the Lattice \code{horizonplot} function
#' are available as display parameters \code{horizon.origin} and
#' \code{horizon.scale}.}
#'
#' }
#'
#' For some of the above plotting-types the \code{groups} display parameter can
#' be used to indicate sample sub-groupings. Its value is supposed to be a
#' factor vector of similar length as the number of samples. In most cases, the
#' groups are shown in different plotting colors and data aggregation
#' operations are done in a stratified fashion.
#'
#' The \code{window} display parameter can be used to aggregate the data prior
#' to plotting. Its value is taken as the number of equal-sized windows along
#' the genomic coordinates of the track for which to compute average values.
#' The special value \code{auto} can be used to automatically determine a
#' reasonable number of windows which can be particularly useful when plotting
#' very large genomic regions with many data points.
#'
#' The \code{aggregation} parameter can be set to define the aggregation
#' function to be used when averaging in windows or across collapsed items. It
#' takes the form of either a function which should condense a numeric vector
#' into a single number, or one of the predefined options as character scalars
#' \code{"mean"}, \code{"median"} or \code{"sum"} for mean, median or
#' summation, respectively. Defaults to computing mean values for each sample.
#' Note that the predefined options can be much faster because they are
#' optimized to work on large numeric tables.
#'
#' @template DataTrack-class_param
#'
#' @name DataTrack-class
#'
#'
#' @return
#' The return value of the constructor function is a new object of class
#' \code{DataTrack} or \code{ReferenceDataTrack}.
#' @section Objects from the class:
#'
#' Objects can be created using the constructor function \code{DataTrack}.
#'
#' @author Florian Hahne
#' @inherit GdObject-class seealso
#'
#' @examples
#' ## Object construction:
#'
#' ## An empty object
#' DataTrack()
#'
#' ## from individual arguments
#' dat <- matrix(runif(400), nrow = 4)
#' dtTrack <- DataTrack(
#'     start = seq(1, 1000, len = 100), width = 10, data = dat,
#'     chromosome = 1, genome = "mm9", name = "random data"
#' )
#'
#' ## from GRanges
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(seq(1, 1000, len = 100),
#'     width = 10
#' ))
#' values(gr) <- t(dat)
#' dtTrack <- DataTrack(range = gr, genome = "mm9", name = "random data")
#'
#' ## from IRanges
#' dtTrack <- DataTrack(
#'     range = ranges(gr), data = dat, genome = "mm9",
#'     name = "random data", chromosome = 1
#' )
#'
#' ## from a data.frame
#' df <- as.data.frame(gr)
#' colnames(df)[1] <- "chromosome"
#' dtTrack <- DataTrack(range = df, genome = "mm9", name = "random data")
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(dtTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#' ## Plotting
#' plotTracks(dtTrack)
#'
#' ## Track names
#' names(dtTrack)
#' names(dtTrack) <- "foo"
#' plotTracks(dtTrack)
#'
#' ## Subsetting and splitting
#' subTrack <- subset(dtTrack, from = 100, to = 300)
#' length(subTrack)
#' subTrack[1:2, ]
#' subTrack[, 1:2]
#' split(dtTrack, rep(1:2, each = 50))
#'
#' ## Accessors
#' start(dtTrack)
#' end(dtTrack)
#' width(dtTrack)
#' position(dtTrack)
#' width(subTrack) <- width(subTrack) - 5
#'
#' strand(dtTrack)
#' strand(subTrack) <- "-"
#'
#' chromosome(dtTrack)
#' chromosome(subTrack) <- "chrX"
#'
#' genome(dtTrack)
#' genome(subTrack) <- "mm9"
#'
#' range(dtTrack)
#' ranges(dtTrack)
#'
#' ## Data
#' values(dtTrack)
#' score(dtTrack)
#'
#' ## coercion
#' as(dtTrack, "data.frame")
#' @importFrom grDevices boxplot.stats
#'
#' @exportClass DataTrack
setClass("DataTrack",
    contains = "NumericTrack",
    representation = representation(data = "matrix", strand = "character"),
    prototype = prototype(
        columns = c("score"),
        name = "DataTrack",
        dp = DisplayPars(
            aggregateGroups = FALSE,
            aggregation = "mean",
            missingAsZero = TRUE,
            alpha.confint = 0.3,
            amount = NULL,
            baseline = NULL,
            box.legend = FALSE,
            box.ratio = 1,
            box.width = NULL,
            grid = FALSE,
            cex.legend = 0.8,
            cex.sampleNames = NULL,
            cex = 0.7,
            coef = 1.5,
            col.baseline = NULL,
            col.confint = NA,
            col.boxplotFrame = .DEFAULT_SHADED_COL,
            col.histogram = .DEFAULT_SHADED_COL,
            col.horizon = NA,
            col.mountain = NULL,
            col.sampleNames = "white",
            col = trellis.par.get("superpose.line")[["col"]],
            collapse = FALSE,
            degree = 1,
            do.out = TRUE,
            evaluation = 50,
            factor = 0.5,
            family = "symmetric",
            fill.confint = NULL,
            fill.histogram = NULL,
            fill.horizon = c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3"),
            fill.mountain = c("#CCFFFF", "#FFCCFF"),
            fontface.legend = NULL,
            fontfamily.legend = NULL,
            fontsize.legend = NULL,
            fontcolor.legend = .DEFAULT_SHADED_COL,
            gradient = RColorBrewer::brewer.pal(9, "Blues"),
            groups = NULL,
            horizon.origin = 0,
            horizon.scale = NULL,
            jitter.x = FALSE,
            jitter.y = FALSE,
            levels.fos = NULL,
            legend = TRUE,
            lineheight.legend = NULL,
            lty.baseline = NULL,
            lty.mountain = NULL,
            lwd.baseline = NULL,
            lwd.mountain = NULL,
            min.distance = 0,
            na.rm = FALSE,
            ncolor = 100,
            notch.frac = 0.5,
            notch = FALSE,
            pch = 20,
            separator = 0,
            showColorBar = TRUE,
            showSampleNames = FALSE,
            size = NULL,
            span = 1 / 5,
            stackedBars = TRUE,
            stats = boxplot.stats,
            transformation = NULL,
            type = "p",
            varwidth = FALSE,
            window = NULL,
            windowSize = NULL,
            ylim = NULL,
            yTicksAt = NULL
        )
    )
)

## Initialize ----------------------------------------------------------------

## Only pass on the stuff to the GdObject initializer

#' @describeIn DataTrack-class  Initialize.
#' @export
setMethod("initialize", "DataTrack", function(.Object, data = matrix(), strand, ...) {
    ## the display parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "DataTrack")
    .Object@data <- data
    if (!missing(strand)) {
        .Object@strand <- unique(strand)
    }
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## The file-based version of the DataTrack class. This will mainly provide a means to dispatch to
## a special 'subset' method which should stream the necessary data from disk.

#' @describeIn DataTrack-class The file-based version of the `DataTrack-class`.
#' @exportClass ReferenceDataTrack
setClass("ReferenceDataTrack", contains = c("DataTrack", "ReferenceTrack"))

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritance has some strange features with regards to method selection

#' @describeIn DataTrack-class  Initialize.
#' @export
setMethod("initialize", "ReferenceDataTrack", function(.Object, stream, reference, mapping = list(),
                                                       args = list(), defaults = list(), ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object = .Object, reference = reference, stream = stream,
        mapping = mapping, args = args, defaults = defaults)
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## Constructor ---------------------------------------------------------------

## Constructor. The following arguments are supported:
##    o range: an object of class 'GRanges' containing the data coordinates, or an object of class data.frame
##       with the two mandatroy columns 'start' and 'end'. The coordinates (i.e., the length of the GRanges object)
##       have to match the columns of the data matrix (see below).
##    o data: a numeric matrix of data points with number of columns equal to the number of coordinates in 'range',
##      or a numeric vector of appropriate length that will be coerced in such a one-row matrix.
##    o start, end, width: numeric vectors of the item start and end coordinates, have to match the columns in 'data'
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o strand: character, display data on the plus strand ("+" or 0), the minus strand ("-" or 1) or both
##      strands ("+-" or "-+" or 2). Currently has to be unique for the whole track.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters

#' @describeIn DataTrack-class Constructor function for `DataTrack-class`
#' @export
DataTrack <- function(range = NULL, start = NULL, end = NULL, width = NULL, data, chromosome, strand, genome,
                      name = "DataTrack", importFunction, stream = FALSE, ...) {
    ## Build a GRanges object from the inputs
    wasGR <- is(range, "GRanges") || is.character(range)
    fromFile <- is.character(range)
    isStream <- FALSE
    .missingToNull(c("strand", "chromosome", "importFunction", "genome"))
    args <- list(strand = strand, chromosome = chromosome, genome = genome)
    defs <- list(strand = "*", chromosome = "chrNA")
    range <- .buildRange(
        range = range, start = start, end = end, width = width, args = args, defaults = defs,
        asIRanges = FALSE, chromosome = chromosome, genome = NA, trackType = "DataTrack",
        importFun = importFunction, stream = stream
    )
    if (is.list(range)) {
        isStream <- TRUE
        slist <- range
        range <- GRanges()
    }
    ## Some default checking
    if (length(unique(strand(range))) > 1) {
        stop("The strand has to be unique for all ranges in a DataTrack object.")
    }
    if (!missing(data) && length(range) > 0) {
        if (is.character(data)) {
            if (!wasGR) {
                stop("Columns indices for the data section are only allowed when 'range' is of class 'GRanges'")
            }
            mt <- is.na(match(data, colnames(values(range))))
            if (any(mt)) {
                warning("Unable to match data columns: ", paste(data[mt], collapse = ","))
            }
            data <- as.data.frame(values(range)[, data[!mt], drop = FALSE])
        }
        if (is.null(dim(data))) {
            dim(data) <- c(1, length(data))
        }
        if (is.matrix(data)) {
            data <- as.data.frame(t(data))
        }
    } else {
        data <- if (ncol(values(range))) as.data.frame(values(range)) else matrix(nrow = 0, ncol = 0)
    }
    data <- .prepareDtData(data, len = length(range))
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- if (length(range) > 0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    }
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    values(range) <- NULL
    if (!isStream) {
        return(new("DataTrack",
            chromosome = chromosome, strand = as.character(strand(range)), range = range,
            name = name, genome = genome, data = data, ...
        ))
    } else {
        ## A bit hackish but for some functions we may want to know which track type we need but at the
        ## same time we do not want to enforce this as an additional argument
        e <- new.env()
        e[["._trackType"]] <- "DataTrack"
        environment(slist[["stream"]]) <- e
        return(new("ReferenceDataTrack",
            chromosome = chromosome, strand = as.character(strand(range)), range = range,
            name = name, genome = genome, data = data, stream = slist[["stream"]], reference = slist[["reference"]],
            mapping = slist[["mapping"]], args = args, defaults = defs, ...
        ))
    }
}

## General accessors ---------------------------------------------------------

#' @describeIn DataTrack-class return the raw data values of the object, i.e.,
#' the data matrix in the `data` slot.
#' @export
setMethod("values", "DataTrack", function(x, all = FALSE) {
    if (sum(dim(x@data)) == 0) {
        x@data
    } else {
        sel <- if (all) rep(TRUE, ncol(x@data)) else seqnames(x) == chromosome(x)
        x@data[, sel, drop = FALSE]
    }
})

#' @describeIn DataTrack-class Replace the data matrix in the `data` slot.
#' @export
setReplaceMethod("values", "DataTrack", function(x, value) {
    if (!is.matrix(value)) {
        if (!is.numeric(value) || length(value) != length(x)) {
            stop("Not numeric or invalid length of replacement vector.")
        } else {
            value <- t(as.matrix(value))
        }
    } else {
        if (!is.numeric(value) || ncol(value) != length(x)) {
            stop("Not numeric or dimensions of replacement value do not match.")
        }
    }
    x@data <- value
    return(x)
})

#' @describeIn DataTrack-class return a vector of strand specifiers for all
#' track items, in the form '+' for the Watson strand, '-' for the Crick
#' strand or '*' for either of the two.
#' @export
setMethod("strand", "DataTrack", function(x) x@strand)

#' @describeIn DataTrack-class replace the strand information for the track items.
#' The replacement value needs to be an appropriate scalar or vector of strand values.
#' @export
setReplaceMethod("strand", "DataTrack", function(x, value) {
    if (!is.character(value) || length(value) != 1 || !value %in% c("+", "-", "*")) {
        stop("Invalid replacement value")
    }
    x@strand <- value
    return(x)
})

#' @describeIn DataTrack-class Split a `DataTrack` object by an appropriate
#' factor vector (or another vector that can be coerced into one).
#' The output of this operation is a list of `DataTrack` objects.
#' @export
setMethod("split", signature("DataTrack"),
    definition = function(x, f, ...) {
        f <- factor(f)
        rs <- as.list(split(ranges(x), f))
        ds <- split(t(values(x)), f)
        nr <- nrow(values(x))
        rnms <- rownames(values(x))
        mapply(function(y, z) {
            x@range <- y
            x@data <- matrix(z, nrow = nr, byrow = TRUE, dimnames = list(rnms, NULL))
            return(x)
        }, rs, ds)
    }
)

## Annotation Accessors ------------------------------------------------------

#' @describeIn DataTrack-class returns NULL since there is no grouping
#' information for the ranges in a `DataTrack`.
#' @export
setMethod("feature", signature(GdObject = "DataTrack"), function(GdObject) NULL)

#' @describeIn DataTrack-class this return the unaltered input object since
#' there is no grouping information for the ranges in a `DataTrack`.
#' @export
setReplaceMethod("feature", signature("DataTrack", "character"), function(GdObject, value) GdObject)

## Consolidate ---------------------------------------------------------------
## Collapse  -----------------------------------------------------------------

#' @importFrom Biobase rowMax rowMin
#' @importFrom matrixStats rowMedians rowMins
.aggregator <- function(GdObject) {
    agFun <- .dpOrDefault(GdObject, "aggregation", "mean")
    if (is.list(agFun)) {
        agFun <- agFun[[1]]
    }
    fun <- if (is.character(agFun)) {
        switch(agFun,
            "mean" = rowMeans,
            "sum" = rowSums,
            "median" = rowMedians,
            "extreme" = function(x) apply(x, 1, .extreme),
            "min" = Biobase::rowMin,
            "max" = Biobase::rowMax,
            rowMeans
        )
    } else {
        if (is.function(agFun)) {
            function(x) apply(x, 1, agFun)
        } else {
            stop(
                "display parameter 'aggregation' has to be a function or a character ",
                "scalar in c('mean', 'median', 'sum', 'min', 'max')"
            )
        }
    }
    return(fun)
}

.runaggregator <- function(GdObject) {
    agFun <- .dpOrDefault(GdObject, "aggregation", "mean")
    if (is.list(agFun)) {
        agFun <- agFun[[1]]
    }
    fun <- if (is.character(agFun)) {
        switch(agFun,
            "mean" = runmean,
            "sum" = runsum,
            "median" = runmed2 <- function(x, k, na.rm = FALSE, ...) {
                na.action <- if (na.rm) {
                    "na.omit"
                } else {
                    "+Big_alternate"
                }
                runmed(x = as.numeric(x), k = k, na.action = na.action)
            },
            "min" = runqmin <- function(x, k, i = 1, ...) {
                runq(x = x, k = k, i = i, ...)
            },
            "max" = runqmax <- function(x, k, i = k, ...) {
                runq(x = x, k = k, i = i, ...)
            },
            runmean
        )
    } else {
        if (is.function(agFun)) {
            # na.rm currently not implemented...
            function(x, k, na.rm = FALSE, endrule = "constant") {
                ans <- vapply(0:(length(x) - k), function(offset) agFun(x[seq_len(k) + offset]), FUN.VALUE = numeric(1))
                ans <- Rle(ans)
                if (endrule == "constant") {
                    j <- (k + 1L) %/% 2L
                    runLength(ans)[1L] <- runLength(ans)[1L] + (j - 1L)
                    runLength(ans)[nrun(ans)] <- runLength(ans)[nrun(ans)] + (j - 1L)
                }
                ans
            }
        } else {
            stop(
                "display parameter 'aggregation' has to be a function or a character ",
                "scalar in c('mean', 'median', 'sum', 'min', 'max')"
            )
        }
    }
    return(fun)
}

## For DataTracks we want to collapse data values using the aggregation function provided by calling .aggregator().
## In addition values can be aggregated over fixed window slices when gpar 'window' is not NULL, and using a sliding
## window approach when 'window' == -1

#' @describeIn DataTrack-class preprocess the track before plotting.
#' This will collapse overlapping track items based on the available resolution
#' and increase the width and height of all track objects to a minimum value
#' to avoid rendering issues. See collapsing for details.
#' @keywords internal
setMethod("collapseTrack", signature(GdObject = "DataTrack"), function(GdObject, diff = .pxResolution(coord = "x"), xrange) {
    if (!length(GdObject)) {
        return(GdObject)
    }
    ## first the data transformation if needed
    values(GdObject) <- score(GdObject)
    collapse <- .dpOrDefault(GdObject, "collapse", FALSE)
    min.width <- .dpOrDefault(GdObject, "min.width", 2)
    min.distance <- max(0, .dpOrDefault(GdObject, "min.distance", 0))
    ## When an averaging window has been set, split the data up into these average chunks
    window <- .dpOrDefault(GdObject, "window")
    windowSize <- .dpOrDefault(GdObject, "windowSize")
    missingAsZero <- .dpOrDefault(GdObject, "missingAsZero", TRUE)
    if (!is.null(window) || collapse) {
        GdObject <- GdObject[, order(range(GdObject))]
    }
    r <- ranges(GdObject)
    drange <- c(floor(xrange[1]), ceiling(xrange[2]))
    if (!is.null(window)) {
        rr <- if (is(r, "GRanges")) ranges(r) else r
        fw <- FALSE
        if (window == "auto") {
            window <- min(
                ncol(values(GdObject)), 1000,
                ceiling(width(range(rr)) / (min.width * diff))
            )
        }
        if (window == "fixed") {
            fw <- TRUE
            window <- 100
        }
        if (!is.numeric(window) || length(window) != 1L) {
            stop("gpar 'window' must be a numeric scalar")
        }
        window <- as.integer(window)
        sc <- values(GdObject)
        agFun <- .aggregator(GdObject)
        if (window == 1) {
            sc <- matrix(agFun(sc), ncol = 1)
            rtmp <- IRanges(start = max(1, drange[1]), end = max(1, drange[2] - 1))
            r <- if (is(r, "GRanges")) GRanges(seqnames = seqnames(r)[1], ranges = rtmp) else rtmp
        } else if (window < 1) {
            if (is.null(windowSize)) {
                windowSize <- (max(GdObject) - min(GdObject)) / 100
            }
            if (windowSize %% 2 != 1) {
                windowSize <- windowSize + 1
            }
            if (missingAsZero) {
                rm <- vector("integer", width(range(range(GdObject))))
            } else {
                rm <- as.integer(rep(NA, width(range(range(GdObject)))))
            }
            ind <- unlist(mapply(function(x, y) x:y, start(GdObject), end(GdObject))) - min(GdObject) + 1
            rm[ind] <- rep(sc[1, ], width(GdObject))
            runAgFun <- .runaggregator(GdObject)
            runwin <- suppressWarnings(runAgFun(Rle(as.numeric(rm)), k = windowSize, endrule = "constant", na.rm = TRUE))
            seqSel <- findRun(as.integer(position(GdObject)) - min(GdObject) + 1, runwin)
            newDat <- matrix(runValue(runwin)[seqSel], nrow = 1)
            if (nrow(sc) > 1) {
                newDat <- rbind(
                    newDat,
                    do.call(rbind, lapply(2:nrow(sc), function(x) {
                        rm[ind] <- rep(sc[x, ], width(GdObject))
                        runwin <- suppressWarnings(runAgFun(Rle(as.numeric(rm)), k = windowSize, endrule = "constant", na.rm = TRUE))
                        seqSel <- findRun(as.integer(position(GdObject)) - min(GdObject) + 1, runwin)
                        runValue(runwin)[seqSel]
                        # suppressWarnings(runValue(runmean(Rle(as.numeric(rm)), k = windowSize, endrule = "constant", na.rm = TRUE)))[seqSel]
                    }))
                )
            }
            sc <- newDat
        } else {
            if (!is.null(window) && window > diff(drange)) {
                window <- diff(drange)
            }
            if (!fw || is.null(windowSize)) {
                windowSize <- diff(drange) %/% window
            } else {
                window <- max(1, diff(drange) %/% windowSize)
            }
            remain <- (diff(drange) - (window * windowSize)) / 2
            ir <- IRanges(
                start = seq(from = drange[1] + remain, to = drange[2] - remain - windowSize, length.out = window),
                width = windowSize
            )
            if (remain > 0) {
                ir <- c(
                    IRanges(start = drange[1], width = ceiling(remain)), ir,
                    IRanges(start = drange[2] - ceiling(remain), width = ceiling(remain))
                )
            }
            ol <- as.matrix(findOverlaps(ir, rr))
            scn <- lapply(split(ol[, 2], ol[, 1]), function(i) agFun(sc[, i, drop = FALSE]))
            scn <- do.call(cbind, scn)
            colnames(scn) <- as.character(unique(ol[, 1]))
            sc <- matrix(NA, ncol = length(ir), nrow = nrow(scn))
            sc[, as.integer(colnames(scn))] <- scn
            r <- if (is(r, "GRanges")) {
                GRanges(
                    seqnames = chromosome(GdObject), ranges = ir,
                    strand = unique(as.character(strand(GdObject)))
                )
            } else {
                ir
            }
        }
        GdObject@range <- r
        GdObject@data <- sc
    }
    ## If groups need to be averaged we have to do it here
    groups <- .dpOrDefault(GdObject, "groups")
    if (!is.null(groups) && .dpOrDefault(GdObject, "aggregateGroups", FALSE)) {
        if (!is.factor(groups)) {
            groups <- factor(groups)
        }
        agFun <- .aggregator(GdObject)
        dat <- values(GdObject)
        rownames(dat) <- groups
        datNew <- do.call(rbind, lapply(levels(groups), function(x) agFun(t(dat[groups == x, , drop = FALSE]))))
        GdObject@data <- datNew
        displayPars(GdObject) <- list(groups = levels(groups))
    }
    ## Compute native coordinate equivalent to 1 pixel and resize
    r <- .resize(r, min.width, diff)
    ## Collapse overlapping ranges (less than minXDist space between them) including the associated attributes using
    ## "|" as separator. For both "strand" and "feature" we take the first available entry, which is not optimal but
    ## seems to be the sanest thing to do here...
    if (collapse) {
        minXDist <- min.distance * diff
        rr <- if (is(r, "GRanges")) ranges(r) else r
        if (minXDist < 1) {
            ## We have to fake smaller ranges because reduce will merge also neighboring ranges
            width(rr) <- width(rr) - 1
            rr <- reduce(rr, min.gapwidth = minXDist)
            width(rr) <- width(rr) + 1
        } else {
            rr <- reduce(r, min.gapwidth = minXDist)
        }
        sc <- values(GdObject)
        if (length(rr) == 1) {
            r <- GRanges(seqnames = 1, strand = strand(GdObject)[1], ranges = rr)
            GdObject@range <- r
            GdObject@data <- matrix(rowMeans(sc, na.rm = TRUE), ncol = 1)
        } else if (length(rr) < length(r)) {
            startInd <- sort(unique(vapply(start(rr), function(x) which(start(r) == x), FUN.VALUE = numeric(1L))))
            st <- strand(GdObject)
            startInd <- if (tail(startInd, 1) == length(r)) c(startInd, length(r) + 1) else c(startInd, length(r))
            vsplit <- split(t(as.data.frame(sc, stringsAsFactors = FALSE)), cut(seq_len(length(r)), startInd, iclude.lowest = TRUE, right = FALSE))
            agFun <- .dpOrDefault(GdObject, "aggregation", "mean")
            if (is.list(agFun)) {
                agFun <- agFun[[1]]
            }
            newScore <- if (is.character(agFun)) {
                switch(agFun,
                    "mean" = lapply(vsplit, function(x) rowMeans(matrix(x, nrow = nrow(sc), byrow = TRUE), na.rm = TRUE)),
                    "sum" = lapply(vsplit, function(x) rowSums(matrix(x, nrow = nrow(sc), byrow = TRUE), na.rm = TRUE)),
                    "median" = lapply(vsplit, function(x) rowMedians(matrix(x, nrow = nrow(sc), byrow = TRUE), na.rm = TRUE)),
                    lapply(vsplit, function(x) rowMeans(matrix(x, nrow = nrow(sc), byrow = TRUE), na.rm = TRUE))
                )
            } else {
                if (is.function(agFun)) {
                    lapply(vsplit, function(x) apply(matrix(x, nrow = nrow(sc), byrow = TRUE), 1, function(y) agFun(y)[1]))
                } else {
                    stop("display parameter 'aggregation' has to be a function or a character ", "scalar in c('mean', 'median', 'sum')")
                }
            }
            newScore <- unlist(newScore)
            r <- GRanges(seqnames = seq_len(length(rr)), strand = st, ranges = rr)
            GdObject@data <- newScore
            GdObject@range <- r
        }
    }
    ## Reconstruct the RangedData object and return
    GdObject@range <- r
    return(GdObject)
})

## Subset --------------------------------------------------------------------

#' @describeIn DataTrack-class subset the items in the `DataTrack` object.
#' This is essentially similar to subsetting of the `GRanges` object in the
#' `range` slot. For most applications, the subset method may be more appropriate.
#' @export
setMethod("[", signature(x = "DataTrack"), function(x, i, j, ..., drop = FALSE) {
    x <- .deepCopyPars(x)
    if (!missing(i)) {
        x@data <- x@data[i, , drop = drop]
        displayPars(x) <- list(groups = as.vector(.dpOrDefault(x, "groups")[i]))
    }
    if (!missing(j)) {
        x@range <- x@range[j, ]
        if (ncol(x@data) > 0) {
            x@data <- x@data[, j, drop = drop]
        }
    }
    return(x)
})

## For DataTracks we cut exactly, and also reduce to the current chromosome unless told explicitly not to

#' @describeIn DataTrack-class Subset a `DataTrack` by coordinates
#' and sort if necessary.
#' @export
setMethod("subset", signature(x = "DataTrack"), function(x, from = NULL, to = NULL, sort = FALSE, drop = TRUE, use.defaults = TRUE, ...) {
    ## Subset to a single chromosome first
    if (drop) {
        csel <- seqnames(x) != chromosome(x)
        if (any(csel)) {
            x <- x[, !csel]
        }
    }
    if (!length(x)) {
        return(x)
    }
    ranges <- if (use.defaults) .defaultRange(x, from = from, to = to) else c(from = ifelse(is.null(from), -Inf, from), to = ifelse(is.null(to), Inf, to))
    x <- x[, start(x) >= ranges["from"] & end(x) <= ranges["to"]]
    if (sort) {
        x <- x[, order(range(x))]
    }
    return(x)
})

## ReferenceDataTracks need to stream the data from file and then pass the results on to the next method

#' @describeIn DataTrack-class Subset a `ReferenceDataTrack` by coordinates
#' and sort if necessary.
#' @export
setMethod("subset", signature(x = "ReferenceDataTrack"), function(x, from, to, chromosome, ...) {
    ## We only need to reach out into the referenced file once if the range is already contained in the object
    if (missing(from) || is.null(from) || missing(to) || is.null(to)) {
        stop("Need both start and end location to subset a ReferenceDataTrack")
    }
    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- Gviz::chromosome(x)
    }
    subRegion <- GRanges(seqnames = chromosome[1], ranges = IRanges(start = from, end = to))
    if (length(ranges(x)) == 0 || !all(overlapsAny(ranges(x), subRegion))) {
        vals <- x@stream(x@reference, subRegion)
        x@range <- vals
        mcols(x@range) <- NULL
        x@data <- .prepareDtData(if (ncol(values(vals))) as.data.frame(values(vals)) else matrix(nrow = 0, ncol = 0), length(vals))
        chromosome(x) <- chromosome[1]
    }
    return(callNextMethod(x = x, from = from, to = to, drop = FALSE, ...))
})

## Position ------------------------------------------------------------------

setMethod("score", signature("DataTrack"), function(x, from = NULL, to = NULL, sort = FALSE, transformation = TRUE, ...) {
    if (!is.null(from) && !is.null(to)) {
        x <- subset(x, from = from, to = to, sort = sort, ...)
    }
    vals <- values(x)
    ## apply data transformation if one is set up
    trans <- .dpOrDefault(x, "transformation")
    if (is.list(trans)) {
        trans <- trans[[1]]
    }
    if (transformation && !is.null(trans)) {
        if (!is.function(trans) || length(formals(trans)) != 1L) {
            stop("gpar 'transformation' must be a function with a single argument")
        }
        test <- trans(vals)
        if (!is.numeric(test) || !is.matrix(test) || !all(dim(test) == dim(vals))) {
            stop(
                "The function in gpar 'transformation' results in invalid output.\n",
                "It has to return a numeric matrix with the same dimensions as the input data."
            )
        }
        vals <- test
    }
    return(vals)
})

## DrawAxis ------------------------------------------------------------------

#' @describeIn DataTrack-class add a y-axis to the title panel of a track.
#' @export
setMethod("drawAxis", signature(GdObject = "DataTrack"), function(GdObject, ...) {
    if (as.logical(.dpOrDefault(GdObject, "legend", FALSE)) && !is.null(.dpOrDefault(GdObject, ".__groupLevels"))) {
        pushViewport(viewport(
            y = 1, height = unit(1, "npc") - unit(.dpOrDefault(GdObject, ".__verticalSpace"), "inches"),
            just = c(0.5, 1)
        ))
        on.exit(popViewport(1))
    }
    type <- match.arg(.dpOrDefault(GdObject, "type", "p"), .PLOT_TYPES, several.ok = TRUE)
    isOnlyHoriz <- length(setdiff(type, "horizon")) == 0
    if (!isOnlyHoriz && .dpOrDefault(GdObject, "showAxis", TRUE)) {
        callNextMethod()
    } else {
        if (.dpOrDefault(GdObject, "showSampleNames", FALSE)) {
            groups <- .dpOrDefault(GdObject, "groups")
            sn <- if (is.null(groups)) rownames(values(GdObject)) else rev(unlist(split(rownames(values(GdObject)), factor(groups))))
            cex.sn <- .dpOrDefault(GdObject, "cex.sampleNames", .dpOrDefault(GdObject, "cex.axis", 1))
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
        }
    }
})

## DrawGD --------------------------------------------------------------------

## Helper function to return the absolute extreme value in a vector
.extreme <- function(x) if (all(is.na(x))) NA else x[which.max(abs(x))]

## Map numeric range to values from 1 to n
.z2icol <- function(z, n, xrange = range(z, na.rm = TRUE)) {
    res <- round((z - xrange[1]) / diff(xrange) * (n - 1)) + 1
    res[res > n] <- n
    res[res < 1] <- 1
    return(res)
}

#' @describeIn DataTrack-class plot the object to a graphics device.
#' The return value of this method is the input object, potentially updated
#' during the plotting operation. Internally, there are two modes in which the
#' method can be called. Either in 'prepare' mode, in which case no plotting is
#' done but the object is preprocessed based on the available space, or in
#' 'plotting' mode, in which case the actual graphical output is created.
#' Since subsetting of the object can be potentially costly, this can be
#' switched off in case subsetting has already been performed before or
#' is not necessary.
#'
#' @importFrom grDevices boxplot.stats
#' @importFrom latticeExtra panel.horizonplot
#' @export
setMethod("drawGD", signature("DataTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, subset = TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if ((is.logical(debug) && debug) || debug == "prepare") {
        browser()
    }
    imageMap(GdObject) <- NULL
    type <- .dpOrDefault(GdObject, "type", "p")
    type <- match.arg(type, .PLOT_TYPES, several.ok = TRUE)
    ## Grouping may be useful for some of the plot types, may be ignored for others
    vals <- values(GdObject)
    groups <- .dpOrDefault(GdObject, "groups")
    if (!is.null(groups) && length(groups) != nrow(vals)) {
        stop("'groups' must be a vector of the same length as the number of rows in the data matrix (", nrow(vals), ")")
    }
    if (!is.null(groups) && !is.factor(groups)) {
        groups <- factor(groups)
    }
    stacked <- .dpOrDefault(GdObject, "stackedBars", FALSE)
    ## The general "col" parameter should be the default for all relevant colors except when there are groups.
    pcols <- .getPlottingFeatures(GdObject)
    ## In prepare mode we collapse the track to allow for aggregation and so on since we need the final data
    ## values to draw the axis.
    if (prepare) {
        if (subset) {
            GdObject <- subset(GdObject, from = minBase, to = maxBase)
        }
        xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
        pushViewport(viewport(xscale = xscale, yscale = c(0, 1), clip = TRUE))
        diff <- .pxResolution(coord = "x")
        GdObject <- collapseTrack(GdObject, diff = diff, xrange = c(minBase, maxBase))
        popViewport(1)
        ## If we have groups and stacked histograms we have to adjust the ylim values, also for regular histograms
        if ("histogram" %in% type) {
            vals <- values(GdObject)
            groups <- rep(groups, ncol(vals))
            ylim <- .dpOrDefault(GdObject, "ylim")
            if (!is.null(groups) && nlevels(groups) > 1) {
                if (ncol(vals)) {
                valsS <- split(vals, groups)
                valsS <- do.call(cbind, lapply(seq_along(valsS), function(i) {
                    tmp <- t(matrix(valsS[[i]], ncol = ncol(vals)))
                    if (ncol(tmp)) {
                        colnames(tmp) <- rep(names(valsS)[i], ncol(tmp))
                    }
                    tmp
                }))
                } else {
                    valsS <- matrix(nrow = nlevels(groups), ncol = 0, dimnames = list(levels(groups)))
                }
                displayPars(GdObject) <- list(".__valsS" = valsS)
                if (stacked == TRUE && is.null(ylim)) {
                    ylim <- suppressWarnings(range(unlist(apply(valsS, 1, function(x) {
                        x <- x[!is.na(x)]
                        sel <- x >= 0
                        tmp <- NULL
                        if (!all(is.na(sel))) {
                            if (any(sel)) {
                                tmp <- c(min(x[sel]), sum(x[sel]))
                            }
                            if (any(!sel)) {
                                tmp <- c(max(x[!sel]), tmp, sum(x[!sel]))
                            }
                        }
                        tmp
                    }))))
                    if (length(type) > 1) {
                        ylim <- range(c(ylim, vals))
                    }
                    displayPars(GdObject) <- list(ylim = ylim)
                }
            } else {
                if (is.null(ylim)) {
                     valsA <- t(vals)
                    ylim <- if (!length(valsA)) c(-1, 1) else c(min(c(0, valsA), na.rm = TRUE), max(valsA, na.rm = TRUE))
                    if (length(type) > 1) {
                        ylim <- range(c(ylim, vals), na.rm = TRUE)
                    }
                    displayPars(GdObject) <- list(ylim = ylim)
                }
            }
        }
        ## If we want a legend we have to figure out how much vertical space is needed
        grps <- .dpOrDefault(GdObject, "groups")
        if (!is.factor(grps)) {
            grps <- factor(grps)
        }
        if (is.null(grps) || nlevels(grps) == 1 || length(setdiff(type, c("gradient", "mountain", "grid", "horizon"))) == 0) {
            displayPars(GdObject) <- list(legend = FALSE)
        }
        if (as.logical(as.logical(.dpOrDefault(GdObject, "legend", FALSE))) && nlevels(grps) > 1) {
            pushViewport(viewport(width = unit(1, "npc") - unit(0.2, "inches"), gp = .fontGp(GdObject, "legend")))
            grps <- levels(grps)
            legInfo <- .legendInfo()[type, , drop = FALSE]
            for (i in colnames(legInfo)) {
                legInfo[, i] <- any(legInfo[, i]) && !any(duplicated(pcols[[i]][seq_along(grps)]))
            }
            legFactors <- sort(names(which(apply(legInfo, 2, any))))
            boxSize <- if (length(setdiff(legFactors, c("col", "cex"))) == 0) 0.1 else 0.3
            spacing <- 0.1
            hspacing <- 0.02
            lengths <- as.numeric(convertUnit(stringWidth(grps), "inches"))
            heights <- as.numeric(convertWidth(stringHeight(grps), "inches"))
            colWidth <- max(lengths + boxSize + spacing * 2)
            availSpace <- vpLocation()$isize
            colNum <- max(1, availSpace["width"] %/% colWidth)
            rowNum <- ceiling(length(grps) / colNum)
            rowHeight <- max(c(heights, 0.1))
            vertSpace <- (rowHeight * rowNum) + (hspacing * (rowNum - 1)) + 0.2
            displayPars(GdObject) <- list(
                ".__verticalSpace" = vertSpace, ".__layoutDims" = c(rowNum, colNum),
                ".__boxSize" = boxSize, ".__spacing" = spacing, ".__groupLevels" = grps,
                ".__legFactors" = legFactors
            )
            popViewport(1)
        }
        return(invisible(GdObject))
    }
    if ((is.logical(debug) && debug) || debug == "draw") {
        browser()
    }
    ## We only proceed if there is something to draw within the ranges, but still may have to add the grid and the legend.
    ## Legend drawing causes another viewport for all the other graphics to be opened and will be called after all other
    ## drawing has finished, hence we call it in on.exit
    if (subset) {
        GdObject <- subset(GdObject, from = minBase, to = maxBase)
    }
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
    ## The optional legend is plotted below the data
    grpLevels <- .dpOrDefault(GdObject, ".__groupLevels")
    if (as.logical(.dpOrDefault(GdObject, "legend", FALSE)) && !is.null(grpLevels)) {
        lSpace <- .dpOrDefault(GdObject, ".__verticalSpace")
        pushViewport(viewport(
            y = 1, height = unit(1, "npc") - unit(lSpace, "inches"),
            just = c(0.5, 1)
        ))
        on.exit({
            popViewport(1)
            cex <- .dpOrDefault(GdObject, "cex.legend", 0.8)
            legFactors <- .dpOrDefault(GdObject, ".__legFactors", character())
            pushViewport(viewport(y = 0, height = unit(lSpace, "inches"), just = c(0.5, 0), gp = .fontGp(GdObject, "legend")))
            pushViewport(viewport(width = unit(1, "npc") - unit(0.1, "inches"), height = unit(1, "npc") - unit(0.1, "inches")))
            boxSize <- .dpOrDefault(GdObject, ".__boxSize")
            spacing <- .dpOrDefault(GdObject, ".__spacing")
            dims <- .dpOrDefault(GdObject, ".__layoutDims")
            for (i in seq_along(grpLevels)) {
                grpLev <- grpLevels[i]
                row <- (((i) - 1) %/% dims[2]) + 1
                col <- (((i) - 1) %% dims[2]) + 1
                pushViewport(viewport(width = 1 / dims[2], height = 1 / dims[1], x = (1 / dims[2]) * (col - 1), y = 1 - ((1 / dims[1]) * (row - 1)), just = c(0, 1)))
                grid.rect(gp = gpar(col = "transparent", fill = .dpOrDefault(GdObject, "background.legend", "transparent")))
                if (length(setdiff(legFactors, c("col"))) == 0) {
                    grid.rect(
                        width = unit(boxSize, "inches"), height = unit(boxSize, "inches"), x = 0, just = c(0, 0.5),
                        gp = gpar(fill = pcols$col[grpLev], col = .DEFAULT_SHADED_COL)
                    )
                } else {
                    if (any(c("pch", "col.symbol") %in% legFactors)) {
                        panel.points(unit(boxSize / 2, "inches"), 0.5, pch = pcols$pch[grpLev], cex = pcols$cex[grpLev], col = pcols$col.symbol[grpLev])
                    }
                    if (any(c("lwd", "lty", "col.lines") %in% legFactors)) {
                        ## panel.lines(unit(c(0,boxSize), "inches"), c(0.5, 0.5), col=pcols$col.line[grpLev], lwd=pcols$lwd[grpLev], lty=pcols$lty[grpLev])
                        grid.lines(unit(c(0, boxSize), "inches"), c(0.5, 0.5), gp = gpar(col = pcols$col.line[grpLev], lwd = pcols$lwd[grpLev], lty = pcols$lty[grpLev]))
                    }
                }
                grid.text(x = unit(boxSize + spacing, "inches"), y = 0.5, just = c(0, 0.5), label = grpLevels[i])
                popViewport(1)
            }
            if (.dpOrDefault(GdObject, "box.legend", FALSE)) {
                grid.rect(width = (1 / dims[2]) * length(grpLevels), x = 0, just = "left", gp = gpar(fill = NA))
            }
            popViewport(2)
        })
    }
    if (!length(GdObject)) {
        if ("g" %in% type) {
            panel.grid(
                h = .dpOrDefault(GdObject, "h", -1), v = .dpOrDefault(GdObject, "v", -1),
                col = .dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty = .dpOrDefault(GdObject, "lty.grid", 1),
                lwd = .dpOrDefault(GdObject, "lwd.grid", 1), alpha = alpha
            )
        }
        return(invisible(GdObject))
    }
    vals <- values(GdObject)
    ylim <- suppressWarnings(.dpOrDefault(GdObject, "ylim", range(vals, na.rm = TRUE, finite = TRUE)))
    if (diff(ylim) == 0) {
        ylim <- ylim + c(-1, 1)
    }
    if (all(is.infinite(ylim))) {
        ylim <- c(0, 1)
    }
    ylimExt <- extendrange(r = ylim, f = 0.05)
    xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
    pushViewport(viewport(xscale = xscale, yscale = ylimExt, clip = TRUE))
    ## The plotting parameters, some defaults from the lattice package first
    plot.symbol <- trellis.par.get("plot.symbol")
    superpose.symbol <- trellis.par.get("superpose.symbol")
    superpose.line <- trellis.par.get("superpose.line")
    groups <- rep(groups, ncol(vals))
    ## For loess calculation we need some settings
    span <- .dpOrDefault(GdObject, "span", 1 / 5)
    degree <- .dpOrDefault(GdObject, "degree", 1)
    family <- .dpOrDefault(GdObject, "family", c("symmetric", "gaussian"))
    evaluation <- .dpOrDefault(GdObject, "evaluation", 50)
    font <- .dpOrDefault(GdObject, "font", if (is.null(groups)) plot.symbol$font else superpose.symbol$font)
    fontface <- .dpOrDefault(GdObject, "fontface", if (is.null(groups)) plot.symbol$fontface else superpose.symbol$fontface)
    fontsize <- .dpOrDefault(GdObject, "fontsize", if (is.null(groups)) plot.symbol$fontsize else superpose.symbol$fontsize)
    ## An optional baseline to be added
    baseline <- .dpOrDefault(GdObject, "baseline")
    lwd.baseline <- .dpOrDefault(GdObject, "lwd.baseline", pcols$lwd[1])
    lty.baseline <- .dpOrDefault(GdObject, "lty.baseline", pcols$lty[1])
    ## The actual plotting values
    pos <- position(GdObject)
    x <- rep(pos + 0.5, each = nrow(vals)) # to align it with ticks position
    y <- as.numeric(vals)
    ## A grid should always be plotted first, so we need to catch this here
    wg <- match("g", type, nomatch = NA_character_)
    if (!is.na(wg)) {
        panel.grid(
            h = .dpOrDefault(GdObject, "h", -1), v = .dpOrDefault(GdObject, "v", -1),
            col = pcols$col.grid, lty = pcols$lty.grid, lwd = pcols$lwd.grid
        )
        type <- type[-wg]
    }
    ## The special type 'mountain' has to be handled separately
    if ("mountain" %in% type) {
        mbaseline <- if (is.null(baseline)) 0 else baseline[1]
        fill.mountain <- .dpOrDefault(GdObject, "fill.mountain", superpose.symbol$fill)[c(1, 2)]
        col.mountain <- .dpOrDefault(GdObject, "col.mountain", pcols$col)[1]
        col.baseline <- .dpOrDefault(GdObject, "col.baseline", col.mountain)[1]
        lwd.mountain <- .dpOrDefault(GdObject, "lwd.mountain", pcols$lwd)[1]
        lty.mountain <- .dpOrDefault(GdObject, "lty.mountain", pcols$lty)[1]
        .panel.mountain(x, y,
            col = col.mountain, fill = fill.mountain, span = span, degree = degree, family = family,
            evaluation = evaluation, lwd = lwd.mountain, lty = lty.mountain, col.line = col.mountain, alpha = alpha,
            baseline = mbaseline
        )
        if (!is.na(mbaseline)) {
            panel.abline(h = mbaseline, col = col.baseline, lwd = lwd.baseline, lty = lty.baseline, alpha = alpha)
        }
    }
    ## The special type 'polygon' has to be handled separately
    if ("polygon" %in% type) {
        mbaseline <- if (is.null(baseline)) 0 else baseline[1]
        fill.mountain <- .dpOrDefault(GdObject, "fill.mountain", superpose.symbol$fill)[c(1, 2)]
        col.mountain <- .dpOrDefault(GdObject, "col.mountain", pcols$col)[1]
        col.baseline <- .dpOrDefault(GdObject, "col.baseline", col.mountain)[1]
        lwd.mountain <- .dpOrDefault(GdObject, "lwd.mountain", pcols$lwd)[1]
        lty.mountain <- .dpOrDefault(GdObject, "lty.mountain", pcols$lty)[1]
        .panel.polygon(x, y,
            col = col.mountain, fill = fill.mountain, lwd = lwd.mountain,
            lty = lty.mountain, col.line = col.mountain, alpha = alpha,
            baseline = mbaseline
        )
        if (!is.na(mbaseline)) {
            panel.abline(h = mbaseline, col = col.baseline, lwd = lwd.baseline, lty = lty.baseline, alpha = alpha)
        }
    }
    ## Also the type 'boxplot' is handled up front
    if ("boxplot" %in% type) {
        diff <- .pxResolution(coord = "x")
        box.ratio <- .dpOrDefault(GdObject, "box.ratio", 1)
        sx <- sort(unique(x))
        sxd <- if (length(sx) == 1) 1 else diff(sx)
        box.width <- .dpOrDefault(GdObject, "box.width", (((min(sxd) * 0.5) / box.ratio) / diff)) * diff
        if (!is.null(groups)) {
            tw <- min(width(GdObject))
            spacer <- diff
            nb <- nlevels(groups)
            bw <- .dpOrDefault(GdObject, "box.width", ((tw - (nb + 2) * spacer) / nb) / diff) * diff
            bcex <- min(pcols$cex[1], (bw / diff) / 20)
            by <- lapply(split(vals, groups), matrix, ncol = ncol(vals))
            for (j in seq_along(by)) {
                nn <- nrow(by[[j]])
                off <- (width(GdObject) - (bw * nb) - ((nb + 2) * spacer)) / 2
                xx <- rep(start(GdObject) + (j * spacer) + (j * bw) + off, each = nn) - (bw / 2)
                .panel.bwplot(xx, as.numeric(by[[j]]),
                    box.ratio = box.ratio, box.width = (bw / 2) / box.ratio, pch = pcols$pch[1],
                    lwd = pcols$lwd[1], lty = pcols$lty[1], fontsize = fontsize,
                    col = pcols$col.histogram, cex = bcex, font = font, fontfamily = font, fontface = fontface,
                    fill = pcols$col[j], varwidth = .dpOrDefault(GdObject, "varwidth", FALSE),
                    notch = .dpOrDefault(GdObject, "notch", FALSE), notch.frac = .dpOrDefault(GdObject, "notch.frac", 0.5),
                    levels.fos = .dpOrDefault(GdObject, "level.fos", sort(unique(xx))),
                    stats = .dpOrDefault(GdObject, "stats", boxplot.stats), coef = .dpOrDefault(GdObject, "coef", 1.5),
                    do.out = .dpOrDefault(GdObject, "do.out", TRUE), alpha = alpha
                )
            }
            diffY <- .pxResolution(coord = "y", 2)
            outline <- apply(vals, 2, range)
            grid.rect(start(GdObject), outline[1, ] - diffY,
                width = width(GdObject), height = abs(outline[2, ] - outline[1, ]) + (2 * diffY),
                gp = gpar(col = .dpOrDefault(GdObject, "col.boxplotFrame", .DEFAULT_SHADED_COL), fill = "transparent", alpha = alpha, lty = "dotted"),
                default.units = "native", just = c("left", "bottom")
            )
        } else {
            bcex <- min(pcols$cex[1], ((box.width * 2) / diff) / 20)
            .panel.bwplot(x, y,
                box.ratio = box.ratio, box.width = box.width, pch = pcols$pch[1],
                lwd = pcols$lwd[1], lty = pcols$lty[1], fontsize = fontsize,
                col = pcols$col.histogram, cex = bcex, font = font, fontfamily = font, fontface = fontface,
                fill = pcols$fill[1], varwidth = .dpOrDefault(GdObject, "varwidth", FALSE),
                notch = .dpOrDefault(GdObject, "notch", FALSE), notch.frac = .dpOrDefault(GdObject, "notch.frac", 0.5),
                levels.fos = .dpOrDefault(GdObject, "level.fos", sort(unique(x))),
                stats = .dpOrDefault(GdObject, "stats", boxplot.stats), coef = .dpOrDefault(GdObject, "coef", 1.5),
                do.out = .dpOrDefault(GdObject, "do.out", TRUE), alpha = alpha
            )
        }
    }
    ## 'histogram' fills up the full range area if its width is > 1
    if ("histogram" %in% type) {
        ylimSort <- sort(ylimExt)
        yy <- if (ylimSort[1] <= 0 && ylimSort[2] >= 0) 0 else ylimSort[1]
        if (!is.null(groups) && nlevels(groups) > 1) {
            valsS <- .dpOrDefault(GdObject, ".__valsS")
            if (stacked) {
                curMinPos <- curMaxPos <- rep(yy, nrow(valsS))
                for (s in seq_len(ncol(valsS))) {
                    if (!all(is.na(valsS[, s]))) {
                        sel <- !is.na(valsS[, s]) & valsS[, s] >= 0
                        yyy <- curMinPos
                        yyy[sel] <- curMaxPos[sel]
                        offset <- yyy
                        offset[offset != yy] <- 0
                        grid.rect(start(GdObject), yyy,
                            width = width(GdObject), height = valsS[, s] - offset,
                            gp = gpar(col = "transparent", fill = pcols$col[colnames(valsS)[s]], lwd = pcols$lwd[1], lty = pcols$lty[1], alpha = alpha), default.units = "native",
                            just = c("left", "bottom")
                        )
                        curMaxPos[sel] <- curMaxPos[sel] + (valsS[sel, s] - offset[sel])
                        curMinPos[!sel] <- curMinPos[!sel] + (valsS[!sel, s] - offset[!sel])
                    }
                }
                diff <- .pxResolution(coord = "x", pcols$lwd[1] + 1)
                tooNarrow <- width(GdObject) < diff
                if (!all(tooNarrow)) {
                    grid.rect(start(GdObject)[!tooNarrow], curMinPos[!tooNarrow],
                        width = width(GdObject)[!tooNarrow],
                        height = (curMaxPos - curMinPos)[!tooNarrow],
                        gp = gpar(fill = "transparent", col = pcols$col.histogram, lwd = pcols$lwd[1], lty = pcols$lty[1], alpha = alpha),
                        default.units = "native", just = c("left", "bottom")
                    )
                }
            } else {
                spacer <- .pxResolution(min.width = 1, coord = "x")
                yOff <- .pxResolution(min.width = 1, coord = "y")
                outline <- apply(valsS, 1, function(x) range(c(yy, x), na.rm = TRUE))
                grid.rect(start(GdObject), outline[1, ] - yOff,
                    width = width(GdObject), height = apply(outline, 2, diff) + (yOff * 2),
                    gp = gpar(col = pcols$col.histogram, fill = pcols$fill.histogram, lwd = pcols$lwd[1], lty = pcols$lty[1], alpha = alpha), default.units = "native",
                    just = c("left", "bottom")
                )
                len <- ncol(valsS)
                subW <- (width(GdObject) - (spacer * (len + 1))) / len
                sel <- subW > spacer
                ## FIXME: how do we treat this if there is not enough space to plot?
                sel <- !logical(length(subW))
                if (any(sel)) {
                    subW <- subW[sel]
                    valsS <- valsS[sel, ]
                    subX <- rep(start(GdObject)[sel], len) + (subW * rep(seq_len(len) - 1, each = sum(sel))) +
                        (spacer * rep(seq_len(len), each = sum(sel)))
                    grid.rect(subX, yy,
                        width = rep(subW, len), height = valsS - yy,
                        gp = gpar(
                            col = "transparent", fill = rep(pcols$col[seq_len(len)], each = sum(sel)),
                            lwd = pcols$lwd[1], lty = pcols$lty[1], alpha = alpha
                        ), default.units = "native",
                        just = c("left", "bottom")
                    )
                }
            }
        } else {
            valsS <- t(vals)
            grid.rect(start(GdObject), yy,
                width = width(GdObject), height = valsS - yy,
                gp = gpar(col = pcols$col.histogram, fill = pcols$fill.histogram, lwd = pcols$lwd[1], lty = pcols$lty[1], alpha = alpha), default.units = "native",
                just = c("left", "bottom")
            )
        }
    }
    ## gradient summarizes the data as a color gradient
    if ("gradient" %in% type) {
        ncolor <- .dpOrDefault(GdObject, "ncolor", 100)
        gradient <- colorRampPalette(.dpOrDefault(GdObject, "gradient", brewer.pal(9, "Blues")))(ncolor)
        valsScaled <- .z2icol(colMeans(vals, na.rm = TRUE), ncolor, sort(ylim))
        grid.rect(start(GdObject), sort(ylim)[1],
            width = width(GdObject), height = abs(diff(ylim)),
            gp = gpar(col = gradient[valsScaled], fill = gradient[valsScaled], alpha = alpha),
            default.units = "native", just = c("left", "bottom")
        )
    }
    ## heatmap does the same, but for each sample individually
    if ("heatmap" %in% type) {
        ncolor <- .dpOrDefault(GdObject, "ncolor", 100)
        valsScaled <- .z2icol(vals, ncolor, sort(ylim))
        nr <- nrow(vals)
        yy <- seq(min(ylim), max(ylim), len = nr + 1)[-1]
        ydiff <- .pxResolution(coord = "y")
        separator <- .dpOrDefault(GdObject, "separator", 0) * ydiff
        if (!is.null(groups)) {
            valsS <- split(vals, groups)
            freq <- table(factor(.dpOrDefault(GdObject, "groups")))
            cmf <- c(0, cumsum(freq))
            for (s in seq_along(valsS))
            {
                gradient <- colorRampPalette(c("white", pcols$col[s]))(ncolor + 5)[-seq_len(5)]
                valsScaled <- .z2icol(valsS[[s]], ncolor, sort(ylim))
                grid.rect(rep(start(GdObject), each = freq[s]), yy[(cmf[s] + 1):cmf[s + 1]],
                    width = rep(width(GdObject), each = freq[s]),
                    height = max(ydiff, abs(diff(ylim)) * (1 / nr) - separator),
                    gp = gpar(col = gradient[valsScaled], fill = gradient[valsScaled], alpha = alpha),
                    default.units = "native", just = c("left", "top")
                )
            }
        } else {
            gradient <- colorRampPalette(.dpOrDefault(GdObject, "gradient", brewer.pal(9, "Blues")))(ncolor)
            grid.rect(rep(start(GdObject), each = nr), rev(yy),
                width = rep(width(GdObject), each = nr),
                height = max(ydiff, abs(diff(ylim)) * (1 / nr) - separator),
                gp = gpar(col = gradient[valsScaled], fill = gradient[valsScaled], alpha = alpha),
                default.units = "native", just = c("left", "top")
            )
        }
    }
    ## For the horizon plot we can use the latticeExtra panel function, but need to reset the y-range
    if ("horizon" %in% type) {
        nband <- 3
        origin <- .dpOrDefault(GdObject, "horizon.origin", 0)
        gr <- if (is.null(groups)) rep(1, nrow(vals)) else factor(.dpOrDefault(GdObject, "groups"))
        yy <- lapply(split(as.data.frame(vals), gr), colMeans, na.rm = TRUE)
        hfill <- .dpOrDefault(GdObject, "fill.horizon", .DEFAULT_HORIZON_COL)
        hcol <- .dpOrDefault(GdObject, "col.horizon", NA)
        separator <- ceiling(.dpOrDefault(GdObject, "separator", 0) / 2)
        pushViewport(viewport(height = 0.95, clip = TRUE))
        for (i in seq_along(yy)) {
            yi <- yy[[i]]
            horizonscale <- .dpOrDefault(GdObject, "horizon.scale", max(abs(yi - origin), na.rm = TRUE) / nband)
            yr <- origin + c(0, horizonscale)
            pushViewport(viewport(y = (i - 1) / length(yy), height = 1 / length(yy), just = c(0.5, 0), clip = TRUE))
            pushViewport(viewport(height = unit(1, "npc") - unit(separator, "points"), clip = TRUE))
            xscale <- if (!.dpOrDefault(GdObject, "reverseStrand", FALSE)) c(minBase, maxBase) else c(maxBase, minBase)
            pushViewport(viewport(xscale = xscale, yscale = yr, clip = TRUE))
            panel.horizonplot(pos, yi, border = hcol, col.regions = hfill)
            popViewport(3)
        }
        popViewport(1)
    }

    ## plot key-value pairs defined here.
    plot_args <- list(
        type = type, groups = groups, pch = pcols$pch,
        col = pcols$col, col.line = pcols$col.line, col.symbol = pcols$col.symbol, fill = pcols$fill,
        font = font, fontfamily = font, fontface = fontface, lty = pcols$lty, cex = pcols$cex, lwd = pcols$lwd, horizontal = FALSE,
        span = span, degree = degree, family = family, evaluation = evaluation,
        jitter.x = .dpOrDefault(GdObject, "jitter.x", FALSE), jitter.y = .dpOrDefault(GdObject, "jitter.y", FALSE),
        factor = .dpOrDefault(GdObject, "factor", 0.5), amount = .dpOrDefault(GdObject, "amount"),
        alpha = alpha
    )

    ## The rest uses the lattice panel function
    na.rm <- .dpOrDefault(GdObject, "na.rm", FALSE)
    sel <- is.na(y)
    if (na.rm && any(sel)) {
        x <- x[!sel]
        y <- y[!sel]
        groups <- groups[!sel]
    }
    plot_args[["x"]] <- x
    plot_args[["y"]] <- y
    plot_args[["groups"]] <- groups
    plot_args[["subscripts"]] <- seq_along(x)

    ## confidence interval bands
    if ("confint" %in% type) {
        ## column-wise SD calculation
        vectorizedSD <- function(mat, na.rm) {
            ssq <- colSums(mat^2, na.rm = na.rm)
            sumel <- colSums(mat, na.rm = na.rm)
            N <- nrow(mat)
            var <- (1 / (N - 1)) * (ssq - (sumel^2) / N)
            return(sqrt(var))
        }
        debugMode <- FALSE

        my.panel.bands <- function(df, col, fill, font, fontface, ...) {
            upper <- df$high
            lower <- df$low
            x <- df$x
            y <- df$y
            na_idx <- which(is.na(upper))
            ## case 1. there are no error bars to plot at all
            if (length(na_idx) == length(upper)) {
                if (debugMode) message("\t Case 1: all empty. returning")
                return(TRUE)
                ## case 2. no missing points
            } else if (length(na_idx) < 1) {
                if (debugMode) message("\t Case 2: one continuous polygon")
                panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                    border = col, col = fill, alpha = alpha, ...
                )
                ## case 3. have complete data with some or no missing points
            } else {
                curr_start <- min(which(!is.na(upper)))
                if (debugMode) message(sprintf("\t Case 3: %i of %i NA", length(na_idx), length(upper)))
                curr_na_pos <- 1
                while (curr_na_pos <= length(na_idx)) {
                    if (debugMode) message(sprintf("\t\tcurr_na_pos = %i, na_idx length= %i", curr_na_pos, length(na_idx)))
                    ## complete the current poly
                    idx <- curr_start:(na_idx[curr_na_pos] - 1)
                    panel.polygon(c(x[idx], rev(x[idx])), c(upper[idx], rev(lower[idx])),
                        col = fill, border = col, alpha = alpha, ...
                    )
                    ## contiguous empty spots - skip
                    while ((na_idx[curr_na_pos + 1] == na_idx[curr_na_pos] + 1) && (curr_na_pos < length(na_idx))) {
                        if (debugMode) message(sprintf("\t\ttight-loop:curr_na_pos = %i", curr_na_pos))
                        curr_na_pos <- curr_na_pos + 1
                    }
                    ## at this point, either we've finished NA spots or the next one is far away.
                    ## In any case start a poly and move to the next NA spot
                    curr_start <- na_idx[curr_na_pos] + 1
                    curr_na_pos <- curr_na_pos + 1
                }
                ## there is one last polygon at the end of the view range
                if (na_idx[length(na_idx)] < length(upper)) {
                    if (debugMode) message("\tWrapping last polygon")
                    idx <- curr_start:length(upper)
                    panel.polygon(c(x[idx], rev(x[idx])), c(upper[idx], rev(lower[idx])),
                        col = fill, border = col, alpha = alpha, ...
                    )
                }
            }
        }

        fill <- .dpOrDefault(GdObject, "fill.confint", pcols$col)
        col <- .dpOrDefault(GdObject, "col.confint", pcols$col)
        alpha <- .dpOrDefault(GdObject, "alpha.confint")
        outg <- NULL

        if (!is.null(groups)) {
            groups <- .dpOrDefault(GdObject, "groups")
            by <- lapply(split(vals, groups), matrix, ncol = ncol(vals))
            mu <- list()
            confint <- list()
            minnie <- Inf
            maxie <- -Inf

            df <- NULL
            outPlot <- NULL
            mu <- list()
            confint <- list()
            xvals <- position(GdObject)

            ## buffer variation to set final limits
            for (j in seq_along(by)) {
                mu[[j]] <- colMeans(by[[j]], na.rm = TRUE)
                locusSD <- vectorizedSD(by[[j]], na.rm)
                confint[[j]] <- 1.96 * (locusSD / sqrt(nrow(by[[j]])))

                curr_low <- mu[[j]] - confint[[j]]
                curr_high <- mu[[j]] + confint[[j]]
                minnie <- min(c(minnie, curr_low))
                maxie <- max(c(maxie, curr_high))
            }

            names(fill) <- NULL
            for (j in seq_along(by)) {
                g <- names(by)[j]
                if (debugMode) message(g)
                df <- data.frame(
                    x = position(GdObject), y = mu[[j]],
                    low = mu[[j]] - confint[[j]], high = mu[[j]] + confint[[j]],
                    groups = factor(g)
                )
                my.panel.bands(df, col[j], fill[j], alpha, ...)
            }
        } else {
            mu <- colMeans(vals, na.rm = TRUE)
            locusSD <- vectorizedSD(vals, na.rm)
            confint <- 1.96 * (locusSD / sqrt(nrow(vals)))

            df <- data.frame(
                x = position(GdObject), y = mu,
                low = mu - confint, high = mu + confint,
                groups = factor(1)
            )

            my.panel.bands(df, col[1], fill[1], alpha, ...)
        }
    }
    do.call("panel.xyplot", plot_args)
    if (!any(c("mountain", "polygon") %in% type) && !is.null(baseline) && !is.na(baseline)) {
        panel.abline(h = baseline, col = pcols$col.baseline, lwd = lwd.baseline, lty = lty.baseline, alpha = alpha)
    }
    popViewport(1)

    return(invisible(GdObject))
})

## SetAs ---------------------------------------------------------------------

setAs("DataTrack", "data.frame", function(from, to) {
    tmp <- as.data.frame(ranges(from))
    colnames(tmp)[1] <- "chromosome"
    tmp <- cbind(tmp, as.data.frame(t(values(from, all = TRUE))))
    return(tmp)
})

setAs("GRanges", "DataTrack", function(from, to) DataTrack(range = from))

## Show ----------------------------------------------------------------------

## A helper function to plot information regarding additional features on other chromosomes
.addFeatInfo <- function(object, addfeat) {
    freqs <- table(seqnames(object))
    freqs <- freqs[setdiff(names(freqs), chromosome(object))]
    nrChr <- length(freqs)
    msg <- sprintf(
        "There %s %s additional annotation feature%s on %s further chromosome%s%s",
        ifelse(addfeat > 1, "are", "is"),
        addfeat,
        ifelse(addfeat > 1, "s", ""),
        nrChr,
        ifelse(nrChr > 1, "s", ""),
        ifelse(nrChr == 1, sprintf(" (%s)", names(freqs)), "")
    )
    if (nrChr > 1) {
        msg <- if (nrChr > 10) {
            c(
                msg, paste("  ", head(names(freqs), 5), ": ", head(freqs, 5), sep = "", collapse = "\n"),
                "  ...", paste("  ", tail(names(freqs), 5), ": ", tail(freqs, 5), sep = "", collapse = "\n")
            )
        } else {
            c(msg, paste("  ", names(freqs), ": ", freqs, " features", sep = "", collapse = "\n"))
        }
        msg <- c(msg, paste(
            "Call seqlevels(obj) to list all available chromosomes",
            "or seqinfo(obj) for more detailed output"
        ))
    }
    return(msg)
}

#' @describeIn DataTrack-class  Show method.
#' @export
setMethod(
    "show", signature(object = "DataTrack"),
    function(object) {
        msg <- sprintf(
            paste("DataTrack '%s'\n| genome: %s\n| active chromosome: %s\n",
                "| positions: %s\n| samples:%s\n| strand: %s",
                sep = ""
            ),
            names(object),
            genome(object),
            chromosome(object),
            length(object),
            nrow(values(object)),
            strand(object)[1]
        )
        addfeat <- ncol(object@data) - length(object)
        if (addfeat > 0) {
            msg <- c(msg, .addFeatInfo(object, addfeat), "Call chromosome(obj) <- 'chrId' to change the active chromosome")
        }
        cat(paste(msg, collapse = "\n"), "\n")
    }
)

#' @describeIn DataTrack-class  Show method.
#' @export
setMethod("show", signature(object = "ReferenceDataTrack"), function(object) {
    .referenceTrackInfo(object, "ReferenceDataTrack")
})
