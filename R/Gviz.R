##----------------------------------------------------------------------------------------------------------------------
## A bunch of package constants
##----------------------------------------------------------------------------------------------------------------------
.DEFAULT_FILL_COL <- "lightgray"
.DEFAULT_OVERPLOT_COL <- "red"
.DEFAULT_LINE_COL <- "black"
.DEFAULT_SHADED_COL <- "#808080"
.DEFAULT_BRIGHT_SHADED_COL <- "#E0E0E0"
.DEFAULT_SYMBOL_COL <- "#0080FF"
.DEFAULT_LINE_COL <- "darkgray"
.PLOT_TYPES <-  c("p", "l", "b", "a", "s", "g", "r", "S",
                  "smooth", "polygon", "horizon", "histogram",
                  "mountain", "h", "boxplot", "gradient", "heatmap", "confint")
.ALIGNMENT_TYPES <- c("coverage", "sashimi", "pileup")
.THIN_BOX_FEATURES <- c("utr", "ncRNA", "utr3", "utr5", "3UTR", "5UTR", "miRNA", "lincRNA",
                        "three_prime_UTR", "five_prime_UTR")
.DEFAULT_HORIZON_COL <- c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3")
##----------------------------------------------------------------------------------------------------------------------



## Check the class and structure of an object
.checkClass <- function (x, class, length = NULL, verbose = FALSE, mandatory = TRUE){
    if (mandatory && missing(x))
        stop("Argument '", substitute(x), "' is missing with no default",
             call. = verbose)
    msg <- paste("'", substitute(x), "' must be an object of class ",
        paste("'", class, "'", sep = "", collapse = " or "),
        sep = "")
    fail <- !any(sapply(class, function(c, y) is(y, c), x))
    if (!is.null(length) && length(x) != length) {
        if (!is.null(x)) {
            fail <- TRUE
            msg <- paste(msg, "of length", length)
        }
    }
    if (fail)
        stop(msg, call. = verbose)
    else invisible(NULL)
}



## We want to deal with chromosomes in a reasonable way. This coerces likely inputs to a unified
## chromosome name as understood by UCSC. Accepted inputs are:
##    - a single integer or a character coercable to one or integer-character combinations
##    - a character, starting with 'chr' (case insensitive)
## Arguments:
##    o x: a character string to be converted to a valid UCSC chromosome name
##    o force: a logical flag, force prepending of 'chr' if missing
## Value: the UCSC character name
.chrName <- function(x, force=FALSE)
{
    if(!getOption("ucscChromosomeNames") || length(x)==0)
        return(as.character(x))
    xu <- unique(x)
    xum <- sapply(xu, function(y){
        xx <- suppressWarnings(as.integer(y))
        if (!is.na(xx))
            y <- xx
        if(is.numeric(y))
            y <- paste("chr", y, sep = "")
        head <- tolower(substring(y, 1, 3)) == "chr"
        if(!head && force){
            y <-  paste("chr", y, sep = "")
            head <- TRUE
        }
        if(!head){
            stop(sprintf(paste("Invalid chromosome identifier '%s'\nPlease consider setting options(ucscChromosomeNames=FALSE)",
                "to allow for arbitrary chromosome identifiers."), y))
        }
        substring(y, 1, 3) <- tolower(substring(y, 1, 3))
        y})
    names(xum) <- xu
    return(as.vector(xum[as.character(x)]))
}


## Make a deep copy of the display parameter environment
.deepCopyPars <- function(GdObject)
{
    oldPars <- displayPars(GdObject, hideInternal=FALSE)
    GdObject@dp <- DisplayPars()
    displayPars(GdObject) <- oldPars
    return(GdObject)
}



## One central place to check which display types result in stacking. This may change at some point for some
## unimplemented types...
## Arguments:
##    o GdObject: an object inheriting from class GdObject
## Value: a logical skalar indicating whether stacking is needed or not
.needsStacking <- function(GdObject) stacking(GdObject) %in% c("squish", "pack", "full")



## Get the coordinates for an HTML image map from the annotationTrack plot.
## Arguments:
##     o coordinates: a numeric matrix of annotation region coordinates (the bounding box if not rectangular)
## Value: valid HTML image map coordinates based on the current device dimensions
.getImageMap <- function(coordinates)
{
    devSize <- devRes()*par("din")
    loc <- vpLocation()
    size <- loc$location[3:4] - loc$location[1:2]
    xscale <- current.viewport()$xscale
    yscale <- current.viewport()$yscale
    fw <- diff(xscale)
    fh <- diff(yscale)
    u2px <- function(x) ((x - xscale[1])/fw *size[1]) + loc$location[1]
    u2py <- function(y) (devSize[2] - loc$location[2]) - ((y - yscale[1])/fh *size[2])
    return(data.frame(x1=u2px(coordinates[,1]), y1=u2py(coordinates[,4]),
                      x2=u2px(coordinates[,3]), y2=u2py(coordinates[,2]),
                      stringsAsFactors=FALSE))
}


.whichStrand <- function(trackList){
    if(!is.list(trackList))
        trackList <- list(trackList)
    str <- unlist(lapply(trackList, function(x)
                         if(is(x, "HighlightTrack") || is(x, "OverlayTrack")) sapply(x@trackList, .dpOrDefault, "reverseStrand") else{
                             .dpOrDefault(x, "reverseStrand")}))
    return(ifelse(str, "reverse", "forward"))
}



## A function returning the amount of vertical space needed for a track
## Arguments:
##    o x: an object inheriting from class GdObject
## Value: the relative vertical space needed for the track
.verticalSpace <- function(x, totalSpace)
{
    if(is(x, "AlignedReadTrack")){
        size <- if(is.null(displayPars(x, "size"))){
            type <- match.arg(.dpOrDefault(x, "detail", "coverage"), c("reads", "coverage"))
            if(type == "read")
                if(stacking(x) %in% c("sqish", "full")) 5 else 1 else 7} else displayPars(x, "size")
        return(size)
    }
    if(is(x, "DataTrack") && is.null(displayPars(x, "size"))){
        type <- match.arg(.dpOrDefault(x, "type", "p"), .PLOT_TYPES, several.ok=TRUE)
        size <- if(length(type)==1L){ if(type=="gradient") 1 else if(type=="heatmap") nrow(values(x)) else 5} else 5
        return(size)
    }
    if(is(x, "GenomeAxisTrack") || is(x, "IdeogramTrack") || is(x, "SequenceTrack"))
    {
        nv <- displayPars(x, "neededVerticalSpace")
        size <- displayPars(x, "size")
        if(is.null(size))
            if(!is.null(nv))
            {
                size <- nv
                attr(size, "absolute") <- TRUE
            } else size <- 1
        return(size)
    }
    size <- .dpOrDefault(x, "size", 1)
    if(is(x, "StackedTrack"))
      size <- max(size, min(floor(vpLocation()$size["height"]/10), size*max(stacks(x))))
    return(size)
}



## Return a particular displayPars value or a default
## Arguments:
##    o GdObject: an object inheriting from class GdObject
##    o par: the name of the displayPar, or a list of alternatives to go though before finally taking
##       the supplied default value
##    o default: a default value for the parameter if it can't be found in GdObject
## Value: the value of the displayPar
.dpOrDefault <- function(GdObject, par, default=NULL, fromPrototype=FALSE)
{
    val <- getPar(x=GdObject, name=par, asIs=TRUE)
    val <- val[!sapply(val, is.null)]
    if(length(val)==0) {
       if (fromPrototype) {
          val <- .parMappings[[GdObject@name]][par]
          val <- val[!sapply(val, is.null)]
          if(length(val)==0)
              val <- NULL
       } else {
          val <- default
       }
    }else{
        val <- val[[1]]
    }
    return(val)
}



## A special version of the above for font settings. This will try to extract the respective parent
## defaults before finally taking the provided default value.
## Arguments:
##    o GdObject: an object inheriting from class GdObject
##    o par: the name of the displayPar. Can also be a vector in which case a number of alternatives
##       is tested, then the parent default for the first element until finally moving to the supplied default
##    o type: the sub-type
##    o default: a default value for the parameter if it can't be found in GdObject
## Value: the value of the displayPar
.dpOrDefaultFont <- function(GdObject, par, type=NULL, default){
    name <- if(is.null(type)) par else sprintf("%s.%s", par, type)
    val <- getPar(GdObject, name[1])
    name <- name[-1]
    while(is.null(val) && length(name)){
        val <- getPar(GdObject, name[1])
        name <- name[-1]
    }
    if(is.null(val))
        val <- .dpOrDefault(GdObject, par[1], default)
    return(val)
}


## Return the font settings for a GdObject
.fontGp <- function(GdObject, subtype=NULL, ...){
    if(is(GdObject, "OverlayTrack"))
        GdObject <- GdObject@trackList[[1]]
    gp <- list(fontsize=as.vector(.dpOrDefaultFont(GdObject, "fontsize", subtype, 12))[1],
               fontface=as.vector(.dpOrDefaultFont(GdObject, "fontface", subtype, 1))[1],
               fontfamily=as.character(as.vector(.dpOrDefaultFont(GdObject, "fontfamily", subtype, 1)))[1],
               col=as.vector(.dpOrDefaultFont(GdObject, "fontcolor", subtype, "black"))[1],
               lineheight=as.vector(.dpOrDefaultFont(GdObject, "lineheight", subtype, 1))[1],
               alpha=as.vector(.dpOrDefaultFont(GdObject, "alpha", subtype, 1))[1],
               cex=as.vector(.dpOrDefaultFont(GdObject, "cex", subtype, 1))[1])
    gp[names(list(...))] <- list(...)
    gp <- gp[!sapply(gp, is.null)]
    return(do.call(gpar, gp))
}


## Check a list of GdObjects whether an axis needs to be drawn for each of them.
## Arguments:
##    o object: a list of GdObjects
## Value: a logical vector of the same length as 'objects'
.needsAxis <- function(objects) {
    if(!is.list(objects))
        objects <- list(objects)
    atrack <- sapply(objects, function(x){
        is(x, "NumericTrack") ||
        (is(x, "AlignmentsTrack") && "coverage" %in% match.arg(.dpOrDefault(x, "type", .ALIGNMENT_TYPES), .ALIGNMENT_TYPES, several.ok=TRUE)) ||
        (is(x, "AlignedReadTrack") && .dpOrDefault(x, "detail", "coverage")=="coverage")
    })
    isOnlyHoriz <- sapply(objects, function(x){
        res <- FALSE
        if(is(x, "DataTrack")){
            type <- match.arg(.dpOrDefault(x, "type", "p"), .PLOT_TYPES, several.ok=TRUE)
            res <- length(setdiff(type, "horizon")) == 0 && !.dpOrDefault(x, "showSampleNames", FALSE)
        }
        res
    })
    return(atrack & sapply(objects, .dpOrDefault, "showAxis", TRUE) & !isOnlyHoriz)
}

## Check a list of GdObjects whether a title needs to be drawn for each of them.
## Arguments:
##    o object: a list of GdObjects
## Value: a logical vector of the same length as 'objects'
.needsTitle <- function(objects) {
    if(!is.list(objects))
        objects <- list(objects)
    sapply(objects, function(x){
        if(is(x, "HighlightTrack") || is(x, "OverlayTrack"))
            any(sapply(x@trackList, .dpOrDefault, "showTitle", TRUE)) else .dpOrDefault(x, "showTitle", TRUE)
    })
}


## Helper function to set up the text size based on the available space
## Arguments:
##    o trackList: a list of GdObjects
##    o sizes: a matching vector of relative vertical sizes
##    o title.width: the available width for the title
## Value: a list with items:
##    o spaceNeeded: the necessary vertical space
##    o cex: the character expansion factor
##    o title.width: the updated available title width
##    o spacing: the amount of spacing between tracks
##    o nwrap: the final (wrapped) title text
.setupTextSize <- function(trackList, sizes, title.width, panelOnly=FALSE, spacing=5)
{
    curVp <- vpLocation()
    trackList <- lapply(trackList, function(x) if(is(x, "OverlayTrack")) x@trackList[[1]] else x)
    spaceNeeded <- if(is.null(sizes)) lapply(trackList, .verticalSpace, curVp$size["height"]) else {
        if(length(sizes) != length(trackList))
            stop("The 'sizes' vector has to match the size of the 'trackList'.")
        rev(sizes)
    }
    whichAbs <- sapply(spaceNeeded, function(x) !is.null(attr(x, "absolute")) && attr(x, "absolute"))
    spaceNeeded <- unlist(spaceNeeded)
    leftVetSpace <- curVp$size["height"]-sum(spaceNeeded[whichAbs])
    spaceNeeded[!whichAbs] <- spaceNeeded[!whichAbs]/sum(spaceNeeded[!whichAbs])*leftVetSpace
    spaceNeeded <- spaceNeeded/sum(spaceNeeded)
    if(!panelOnly)
    {
        ## Figure out the fontsize for the titles based on available space. If the space is too small (<atLeast)
        ## we don't plot any text, and we also limit to 'maximum' to avoid overblown labels. If the displayPars
        ## 'cex.title' or 'cex.axis are not NULL, those override everything else.
        nn <- sapply(trackList, names)
        nwrap <- sapply(nn, function(x) paste(strwrap(x, 10), collapse="\n"))
        needAxis <- .needsAxis(trackList)
        isOnlyHoriz <- sapply(trackList, function(x){
            res <- FALSE
            if(is(x, "DataTrack")){
                type <- match.arg(.dpOrDefault(x, "type", "p"), .PLOT_TYPES, several.ok=TRUE)
                res <- length(setdiff(type, "horizon")) == 0
            }
            res
        })
        nwrap[needAxis] <- nn[needAxis]
        lengths <- as.numeric(convertWidth(stringWidth(nwrap),"inches"))+0.2
        heights <- curVp$isize["height"]*spaceNeeded
        atLeast <- 0.5
        maximum <- 1.2
        allCex <- heights/lengths
        parCex <- sapply(trackList, function(x) if(is.null(displayPars(x, "cex.title"))) NA else displayPars(x, "cex.title"))
        toSmall <- allCex<atLeast
        cex <- rep(max(c(atLeast, min(c(allCex[!toSmall], maximum), na.rm=TRUE))), length(allCex))
        cex[!is.na(parCex)] <- parCex[!is.na(parCex)]
        lengthsNoWrap <- as.numeric(convertWidth(stringWidth(nn),"inches"))*cex+0.4
        needsWrapping <- lengthsNoWrap>heights
        nwrap[allCex<atLeast & is.na(parCex)] <- ""
        nwrap[!needsWrapping] <- nn[!needsWrapping]
        ## Figure out the title width based on the available tracks. If there is an axis for at least one of the tracks,
        ## we add 1.5 time the width of a single line of text to accomodate for it.
        wfac <- curVp$isize["width"]
        leaveSpace <- ifelse(any(needAxis), 0.15, 0.2)
        width <- (as.numeric(convertHeight(stringHeight(paste("g_T", nwrap, "g_T", sep="")), "inches"))*cex+leaveSpace)/wfac
        showtitle <- sapply(trackList, .dpOrDefault, "showTitle", TRUE)
        width[!showtitle & !needAxis] <- 0
        width[!showtitle] <- 0
        twfac <- if(missing(title.width) || is.null(title.width)) 1 else title.width
        title.width <- max(width, na.rm=TRUE)
        if(any(needAxis)){
            cex.axis <- structure(sapply(trackList, .dpOrDefault, "cex.axis", 0.6), names=nn)
            axTicks <- unlist(lapply(trackList, function(GdObject){
                if(!is(GdObject, "NumericTrack") && !is(GdObject, "AlignedReadTrack") && !is(GdObject, "AlignmentsTrack"))
                    return(NULL)
                yvals <- if(is(GdObject, "AlignedReadTrack")) runValue(coverage(GdObject, strand="*")) else values(GdObject)
                ylim <- .dpOrDefault(GdObject, "ylim", if(!is.null(yvals) && length(yvals))
                                     range(yvals, na.rm=TRUE, finite=TRUE) else c(-1,1))
                if(diff(ylim)==0)
                    ylim <- ylim+c(-1,1)
                yscale <- extendrange(r=ylim, f=0.05)
                at <- pretty(yscale)
                at[at>=sort(ylim)[1] & at<=sort(ylim)[2]]
                atSpace <- max(as.numeric(convertWidth(stringWidth(at), "inches"))+0.18)*cex.axis[names(GdObject)]
                if(is(GdObject, "DataTrack")){
                    type <- match.arg(.dpOrDefault(GdObject, "type", "p"), .PLOT_TYPES, several.ok=TRUE)
                    if(any(c("heatmap", "gradient") %in% type)){
                        nlevs <- max(1, nlevels(factor(getPar(GdObject, "groups"))))-1
                        atSpace <- atSpace + 0.3 * atSpace + as.numeric(convertWidth(unit(3, "points"), "inches"))*nlevs
                    }
                    if(type %in% c("heatmap", "horizon") && .dpOrDefault(GdObject, "showSampleNames", FALSE)){
                        sn <- rownames(values(GdObject))
                        axSpace <- ifelse(isOnlyHoriz, 0, 10)
                        wd <- max(as.numeric(convertWidth(stringWidth(sn) + unit(axSpace, "points"), "inches")))
                        atSpace <- atSpace + (wd * .dpOrDefault(GdObject, "cex.sampleNames", 0.5))
                    }
                }
                atSpace
            }))
            hAxSpaceNeeded <- (max(axTicks))/wfac
            title.width <- title.width +  hAxSpaceNeeded
        }
    } else {
        title.width <- nwrap <- cex <- NA
    }
    spacing <- as.numeric(convertWidth(unit(spacing, "points"), "npc"))
    title.width <- title.width * twfac
    return(list(spaceNeeded=spaceNeeded, cex=cex, title.width=title.width, spacing=spacing, nwrap=nwrap))
}



## This coerces likely inputs for the genomic strand  to a unified
## strand name. Accepted inputs are:
##    o a single integer, where values <=0 indicate the plus strand, and values >=1 indicate the minus strand
##    o a character, either "+" or "-"
## If extended=TRUE, the additional values 2, "+-", "-+" and "*" are allowed, indicating to use both strands.
## Value: the validated strand name
.strandName <- function(x, extended=FALSE)
{
  fun <- function(x, extended){
    if(!extended) {
      if(is.numeric(x)) {
        x <- min(c(1,max(c(0, as.integer(x)))))
      } else if(is.character(x)) {
        x <- match(x, c("+", "-"))-1
        if(any(is.na(x)))
          stop("The strand has to be specified either as a character ('+' or '-'), or as an integer value (0 or 1)")
      }
    } else {
      if(is.numeric(x)) {
        x <- min(c(2,max(c(0, as.integer(x)))))
      } else if(is.character(x)) {
        x <- min(c(2, match(x, c("+", "-", "+-", "-+", "*"))-1))
        if(any(is.na(x)))
          stop("The strand has to be specified either as a character ('+' or '-'), or as an integer value (0 or 1)")
      }
    }
    x
  }
  return(sapply(x, fun, extended))
}



## Compute native coordinate equivalent to 'min.width' pixel. This assumes that a graphics device is already
## open, otherwise a new window will pop up, which could be a little annoying.
## Arguments:
##    o min.width: the number of pixels
##    o coord: the axis for which to compute the coordinats, one in c("x","y")
## Value: the equivalent of 'min.width' in native coordinates.
.pxResolution <- function(min.width=1, coord=c("x","y"))
{
    coord <- match.arg(coord, several.ok=TRUE)
    curVp <- vpLocation()
    co <- c(x=as.vector(abs(diff(current.viewport()$xscale))/(curVp$size["width"])*min.width),
            y=as.vector(abs(diff(current.viewport()$yscale))/(curVp$size["height"])*min.width))
    return(co[coord])
}


## Deal with composite exons and provide polygon plotting coordinates for them. For all normal exons simply return
## the original bounding box data
## Arguments:
##    o box: the data frame with the bounding box information
##    o type: the type of coordinates for the composite exons, one in 'box', 'arrow' or 'fixedArrow'
## Value: a list with three elements: box, the non-composite exons, pols, the plotting coordinates for the merged composite exons,
##        and polpars, a data.frame of plotting parameters for the polygons
.handleComposite <- function(box, type="box", W=1/4, H=1/3, min.width=10, max.width=Inf){
    box <- box[order(box$start),]
    boxFinal <- data.frame(stringsAsFactors=FALSE)
    polFinal <- data.frame(stringsAsFactors=FALSE)
    polPars <- data.frame(stringsAsFactors=FALSE)
    bss <- split(box, box$transcript)
    for(b in bss){
        ol <- names(which(table(b$exon)>1))
        boxFinal <- rbind(boxFinal, b[!b$exon %in% ol,])
        if(length(ol)){
            b <- b[b$exon %in% ol,]
            r <- IRanges(start=b$cx1, end=b$cx2)
            rr <- reduce(r)
            brs <- split(b, subjectHits(findOverlaps(r, rr)))
            for(j in seq_along(brs)){
                if(nrow(brs[[j]]) == 1){
                    boxFinal <- rbind(boxFinal, brs[[j]])
                }else{
                    xlocs <- as.vector(t(brs[[j]][, c("cx1", "cx2"), drop=FALSE]))
                    ylocs <- unlist(brs[[j]][, c("cy1", "cy2"), drop=FALSE])
                    fh <- 1:(length(ylocs)/2)
                    sh <- (length(ylocs)/2+1):length(ylocs)
                    ylocs[sh] <- rev(ylocs[sh])
                    ylocs <- rep(ylocs, each=2)
                    if(type == "box"){
                        polFinal <- rbind(polFinal, data.frame(x=c(xlocs, rev(xlocs)), y=ylocs,
                                                               id=paste(brs[[j]][1, "transcript"], brs[[j]][1, "exon"], j),
                                                               stringsAsFactors=FALSE))
                        polPars <- rbind(polPars, data.frame(fill=brs[[j]][1, "fill"], col=brs[[j]][1, "col"], stringsAsFactors=FALSE))
                    }else{
                        polPars <- rbind(polPars, data.frame(fill=brs[[j]][1, "fill"], col=brs[[j]][1, "col"], stringsAsFactors=FALSE))
                        str <- brs[[j]]$strand[1]
                        offset <- (abs(brs[[j]]$cy1 - brs[[j]]$cy2)) * H / 2
                        if(str == "+" && abs(diff(xlocs[(length(xlocs)-1):length(xlocs)])) > min.width){
                            offset <- rep(offset, each=2)
                            asel <- (length(xlocs)-1):length(xlocs)
                            bsel <- 1:(min(asel)-1)
                            d <- abs(diff(xlocs[asel]))
                            afp <- if(type == "arrow") rep(xlocs[asel][1] + max(d*W, d-max.width), 2) else {
                                rep(max(xlocs[asel][1], xlocs[asel][2]-W), 2)}
                            xlocs <- c(xlocs[bsel], xlocs[asel][1], afp, xlocs[asel][2], afp, xlocs[asel][1], rev(xlocs[bsel]))
                            mid <- length(ylocs)/2
                            asel <- (mid-1):(mid+2)
                            bsel <- c(1:(mid-2), (mid+3):length(ylocs))
                            fh <- 1:(length(bsel)/2)
                            sh <- (length(bsel)/2+1):length(bsel)
                            ylocs <- c(ylocs[bsel[fh]] + offset[fh], ylocs[asel][1:2] + tail(offset,1), ylocs[asel][1],
                                       ylocs[asel][1] + abs(diff(ylocs[asel][c(1,3)]))/2, ylocs[asel][4],
                                       ylocs[asel][3:4] - tail(offset,1), ylocs[bsel[sh]] - offset[fh])
                            polFinal <- rbind(polFinal, data.frame(x=c(xlocs, rev(xlocs)), y=ylocs,
                                                                   id=paste(brs[[j]][1, "transcript"], brs[[j]][1, "exon"], j),
                                                                   stringsAsFactors=FALSE))
                        }else if(str == "-" && abs(diff(xlocs[1:2])) > min.width){
                            yoffset <- c(rep(offset[-1], each=2), -rev(rep(offset[-1], each=2)))
                            asel <- 1:2
                            bsel <-  3:length(xlocs)
                            d <- abs(diff(xlocs[asel]))
                            afp <- if(type == "arrow") rep(xlocs[asel][2] - max(d*W, d-max.width), 2) else{
                                rep(min(xlocs[asel][2], xlocs[asel][1] + W), 2)}
                            xlocs <- c(xlocs[asel][1], afp, xlocs[asel][2], xlocs[bsel], rev(xlocs[bsel]), xlocs[asel][2], afp, xlocs[asel][1])
                            asel <- c(1:2, (length(ylocs)-1):length(ylocs))
                            bsel <- 3:(length(ylocs)-2)
                            ylocs <- c(ylocs[asel][1] + abs(diff(ylocs[asel][c(1,3)]))/2, ylocs[asel][1], rep(ylocs[asel][1]+offset[1], 2),
                                       ylocs[bsel]+yoffset, rep(ylocs[asel][3]-offset[1], 2), ylocs[asel][3], ylocs[asel][1] + abs(diff(ylocs[asel][c(1,3)]))/2)
                            polFinal <- rbind(polFinal, data.frame(x=xlocs, y=ylocs,
                                                                   id=paste(brs[[j]][1, "transcript"], brs[[j]][1, "exon"], j),
                                                                   stringsAsFactors=FALSE))
                        }else{
                            offset <- c(rep(offset, each=2), -rev(rep(offset, each=2)))
                            polFinal <- rbind(polFinal, data.frame(x=c(xlocs, rev(xlocs)), y=ylocs+offset,
                                                                   id=paste(brs[[j]][1, "transcript"], brs[[j]][1, "exon"], j),
                                                                   stringsAsFactors=FALSE))
                        }
                    }
                }
            }
        }
    }
    return(list(box=boxFinal, pols=polFinal, polpars=polPars))
}


## Take coordinates for the bounding boxes of annotation regions and plot filled arrows inside.
## Arguments:
##    o the data frame with the bounding box information
##    o W: the proportion of the total box width used for the arrow head
##    o H: the proportion of the total box height used for the arrow head
##    o lwd: the boundary line width
##    o lty: the boundary line type
##    o alpha: the transparency
##    o min.width: the minumum width of the arrow head. Below this size a simple box is drawn
## Note that the last arguments 4-9 all have to be of the same length as number of rows in box.
## Value: the function is called for its side-effects of drawing on the graphics device
.filledArrow <- function(box, W=1/4, H=1/3, lwd, lty, alpha, min.width=10, max.width=Inf, absoluteWidth=FALSE) {
    boxC <-  if("transcript" %in% colnames(box)){
        .handleComposite(box, ifelse(absoluteWidth, "fixedArrow", "arrow"), min.width=min.width, max.width=max.width, W=W, H=H) }else {
        list(box=box, pols=data.frame())}
    xx <- yy <- numeric()
    id <- character()
    pars <- data.frame()
    if(nrow(boxC$box)){
        box <- boxC$box
        A <- box[,1:2,drop=FALSE]
        B <- box[,3:4,drop=FALSE]
        ## First everything that is still a box
        osel <- abs(B[,1]-A[,1]) < min.width | !box$strand %in% c("+", "-")
        xx <-  c(A[osel,1], B[osel,1], B[osel,1], A[osel,1])
        offset <- (abs(B[osel,2]-A[osel,2])*H/2)
        yy <- c(rep(A[osel,2]+offset, 2), rep(B[osel,2]-offset, 2))
        id <- rep(seq_len(sum(osel)), 4)
        pars <- data.frame(fill=box$fill, col=box$col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE)[osel,]
        ## Now the arrows facing right
        sel <- !osel & box$strand=="+"
        id <- c(id, rep(seq(from=if(!length(id)) 1 else max(id)+1, by=1, len=sum(sel)), 7))
        d <- abs(B[sel,1]-A[sel,1])
        alp <- if(!absoluteWidth) rep(A[sel,1]+pmax(d*W, d-rep(max.width, sum(sel))),2) else rep(pmax(B[sel,1]-W, pmin(B[sel,1], A[sel,1])),2)
        xx <- c(xx, A[sel,1], alp, B[sel,1], alp, A[sel,1])
        yy <- c(yy, rep(A[sel,2]+(abs(B[sel,2]-A[sel,2])*H/2),2), A[sel,2], A[sel,2]+(abs(B[sel,2]-A[sel,2])/2), B[sel,2],
                rep(B[sel,2]-(abs(B[sel,2]-A[sel,2])*H/2),2))
        pars <- rbind(pars, data.frame(fill=box$fill, col=box$col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE)[sel,])
        ## And finally those facing left
        sel <- !osel & box$strand=="-"
        id <- c(id, rep(seq(from=if(!length(id)) 1 else max(id)+1, by=1, len=sum(sel)), 7))
        d <- abs(B[sel,1]-A[sel,1])
        alp <- if(!absoluteWidth) rep(B[sel,1]-pmax(d*W, d-rep(max.width, sum(sel))), 2) else rep(pmin(A[sel,1]+W, pmax(B[sel,1], A[sel,1])),2)
        xx <- c(xx, B[sel,1], alp, A[sel,1], alp, B[sel,1])
        yy <- c(yy, rep(A[sel,2]+(abs(B[sel,2]-A[sel,2])*H/2),2), A[sel,2], A[sel,2]+(abs(B[sel,2]-A[sel,2])/2), B[sel,2],
                rep(B[sel,2]-(abs(B[sel,2]-A[sel,2])*H/2),2))
        pars <- rbind(pars, data.frame(fill=box$fill, col=box$col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE)[sel,])
    }
    if(nrow(boxC$pols)){
        xx <- c(xx, boxC$pols$x)
        yy <- c(yy, boxC$pols$y)
        id <- c(id, paste("pols", boxC$pols$id))
        pars <- rbind(pars, data.frame(fill=boxC$polpars$fill, col=boxC$polpars$col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE))
    }
    grid.polygon(x=xx, y=yy, gp=gpar(fill=pars$fill, col=pars$col, alpha=pars$alpha, lwd=pars$lwd, lty=pars$lty), default.units="native",
                 id=factor(id))
}


## Take coordinates for the bounding boxes of annotation regions and plot boxes, also making sure that
## composite exons in GeneRegionTracks (i.e., overlapping coordinates and same exon id) are merged
## appropriately
## Arguments:
##    o box: the data frame with the bounding box information
##    o lwd: the boundary line width
##    o lty: the boundary line type
##    o alpha: the transparency
## Value: the function is called for its side-effects of drawing on the graphics device
.filledBoxes <- function(box, lwd, lty, alpha){
    if("transcript" %in% colnames(box)){
        box <- .handleComposite(box, "box")
        if(nrow(box$box))
            grid.rect(box$box$cx2, box$box$cy1, width=box$box$cx2 - box$box$cx1, height=box$box$cy2 - box$box$cy1,
                      gp=gpar(col=as.character(box$box$col), fill=as.character(box$box$fill), lwd=lwd, lty=lty, alpha=alpha),
                      default.units="native", just=c("right", "bottom"))
        if(nrow(box$pols))
            grid.polygon(x=box$pols$x, y=box$pols$y, id=factor(box$pols$id),
                         gp=gpar(col=as.character(box$polpars$col), fill=as.character(box$polpars$fill), lwd=lwd, lty=lty, alpha=alpha),
                         default.units="native")
    }else{
        grid.rect(box$cx2, box$cy1, width=box$cx2 - box$cx1, height=box$cy2 - box$cy1,
                  gp=gpar(col=as.character(box$col), fill=as.character(box$fill), lwd=lwd, lty=lty, alpha=alpha),
                  default.units="native", just=c("right", "bottom"))
    }
}


## Take start and end coordinates for genemodel-type annotations and draw a featherd line indicating
## the strand direction
## Arguments:
##    o xx1, xx2: integer vectors of equal length indicating the start and end of the gene models.
##    o strand: the strand information for each gene model. Needs to be of the same length as xx1 and xx2
##    o coords: the coordinates of the exon features, needed to avoid overlaps.
##    o y: the y value for the arrow bar, usually not set since it should always be 20
##    o W: the width of the arrow feathers in pixels
##    o D: the distance between arrow feathers in pixels
##    o H: the height of the arrow feathers in native coordinates (the total bounding box is usually 40)
##    o col: the boundary color
##    o lwd: the boundary line width
##    o lty: the boundary line type
##    o alpha: the transparency
##    o barOnly: only plot the bar, not the feathers
##    o diff: the current pixel resolution
##    o min.height: the minimum total height in pixels for the feathers (i.e., min.height/2 in each direction)
## Value: the function is called for its side-effects of drawing on the graphics device
.arrowBar <- function(xx1, xx2, strand, coords, y=20, W=3, D=10, H, col, lwd, lty, alpha, barOnly=FALSE,
                      diff=.pxResolution(coord="y"), min.height=3){
    exons <- IRanges(start=coords[,1], end=coords[,3])
    levels <- split(exons, coords[,2]%/%1)
    if(!barOnly){
        onePx <- diff
        if(missing(H)){
            onePy <- .pxResolution(coord="y")
            H <- onePy*min.height/2
        }
        fx1 <- fx2 <- scol <- fy1 <- fy2 <- NULL
        for(i in seq_along(xx1)){
            x1 <- xx1[i]
            x2 <- xx2[i]
            len <- diff(c(x1,x2))/onePx
            if(len>D+W*2){
                ax1 <- seq(from=x1+(onePx*W), to=x1+(len*onePx)-(onePx*W), by=onePx*D)
                ax2 <- ax1+(onePx*W)
                feathers <- IRanges(start=ax1-onePx, end=ax2+onePx)
                cur.level <- y[i]%/%1
                sel <- queryHits(findOverlaps(feathers, resize(levels[[cur.level]], width=width(levels[[cur.level]])-1)))
                if(length(sel)){
                    ax1 <- ax1[-sel]
                    ax2 <- ax2[-sel]
                }
                fx1 <- c(fx1, rep(if(strand[i]=="-") ax1 else ax2, each=2))
                fx2 <- c(fx2, rep(if(strand[i]=="-") ax2 else ax1, each=2))
                scol <- c(scol, rep(col[i], length(ax1)*2))
                fy1 <- c(fy1, rep(rep(y[i], length(ax1)*2)))
                fy2 <- c(fy2, rep(c(y[i]-H, y[i]+H), length(ax1)))
            }
        }
        if(!is.null(fx1) && length(fx1))
            grid.segments(fx1, fy1, fx2, fy2, default.units="native", gp=gpar(col=scol, lwd=lwd, lty=lty, alpha=alpha))
    }
    bars <- data.frame(x1=xx1, x2=xx2, y=y, col=col, stringsAsFactors=FALSE)
    cutBars <- data.frame()
    for(i in 1:nrow(bars)){
        b <- bars[i,]
        cur.level <- b$y%/%1
        ct <- setdiff(IRanges(start=b$x1, end=b$x2), resize(levels[[cur.level]], width=width(levels[[cur.level]])-1))
        if(length(ct))
            cutBars <- rbind(cutBars, data.frame(x1=start(ct), x2=end(ct), y=b$y, col=b$col, stringsAsFactors=FALSE))
    }
    ## fix bug when no introns are present
    if (nrow(cutBars))
        grid.segments(cutBars$x1, cutBars$y, cutBars$x2, cutBars$y, default.units="native", gp=gpar(col=cutBars$col, lwd=lwd, lty=lty, alpha=alpha, lineend="square"))
    ##grid.segments(xx1, y, xx2, y, default.units="native", gp=gpar(col=col, lwd=lwd, lty=lty, alpha=alpha, lineend="square"))
}



## Extract track color for different subtypes within the track and use the default
## color value if no other is found, lightblue if no colors are set at all
## Arguments:
##    o GdObject: object inheriting from class GdObject
## Value: a color character
.getBiotypeColor <- function(GdObject) {
    defCol <- .dpOrDefault(GdObject, "fill", .DEFAULT_FILL_COL)
    col <- sapply(as.character(values(GdObject)[, "feature"]),
                  function(x) .dpOrDefault(GdObject, x)[1], simplify=FALSE)
    needsDef <- sapply(col, is.null)
    col[needsDef] <- rep(defCol, sum(needsDef))[1:sum(needsDef)]
    return(unlist(col))
}


## Compute pretty tickmark location (code from tilingArray package)
## Arguments:
##    o x: a vector of data values
## Value: the tick mark coordinates
.ticks <- function(x){
  rx <- range(x)
  lz <- log((rx[2]-rx[1])/3, 10)
  fl <- floor(lz)
  if( lz-fl > log(5, 10))
    fl <- fl +  log(5, 10)
  tw <- round(10^fl)
  i0 <- ceiling(rx[1]/tw)
  i1 <- floor(rx[2]/tw)
  seq(i0, i1)*tw
}



## A lattice-style panel function to draw smoothed 'mountain' plots
## Arguments:
##    o x, y: the x and y coordinates form the plot
##    o span, degree, family, evaluation: parameters that are passed on to loess
##    o lwd, lty, col: color, with and type of the plot lines
##    o fill: fill colors for areas above and under the baseline, a vector of length two
##    o col.line: color of the baseline
##    o baseline: the y value of the horizontal baseline
##    o alpha: the transparancy
## Value: the function is called for its side-effect of drawing on the graphics device
.panel.mountain <- function (x, y, span=2/3, degree=1, family=c("symmetric", "gaussian"), evaluation=50,
                             lwd=plot.line$lwd, lty=plot.line$lty, col, col.line=plot.line$col,
                             baseline, fill, alpha=1, ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  fill <- rep(fill,2)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 1)
    return()
  if (!missing(col)) {
    if (missing(col.line))
      col.line <- col
  }
  plot.line <- trellis.par.get("plot.line")
  smooth <- loess.smooth(x[ok], y[ok], span = span, family = family,
                         degree = degree, evaluation = evaluation)
  tmp <- as.integer(smooth$y<baseline)
  changePoint = NULL
  for(i in seq_along(tmp))
    if(i>1 && tmp[i]!= tmp[i-1])
      changePoint = c(changePoint, i)
  m <- (smooth$y[changePoint] - smooth$y[changePoint-1]) / (smooth$x[changePoint] - smooth$x[changePoint-1])
  xCross <- ((baseline-smooth$y[changePoint-1])/m) + smooth$x[changePoint-1]
  newX <- newY <- NULL
  j <- 1
  xx <- smooth$x
  yy <- smooth$y
  smooth$x <- c(smooth$x, tail(smooth$x,1))
  smooth$y <- c(smooth$y, baseline)
  xvals <- smooth$x[1]
  yvals <- baseline
  for(i in seq_along(smooth$x)) {
    if(i==length(smooth$x)) {
      xvals <- c(xvals, smooth$x[i])
      yvals <- c(yvals, baseline)
      fcol <- if(mean(yvals)<baseline) fill[1] else fill[2]
      panel.polygon(xvals, yvals, fill=fcol, col=fcol, border=fcol, alpha=alpha)
    } else if(i %in% changePoint) {
      xvals <- c(xvals, xCross[j])
      yvals <- c(yvals, baseline)
      fcol <- if(mean(yvals)<baseline) fill[1] else fill[2]
      panel.polygon(xvals, yvals, fill=fcol, col=fcol, border=fcol, alpha=alpha)
      xvals <- c(xCross[j], smooth$x[i])
      yvals <- c(baseline, smooth$y[i])
      j <- j+1
    } else {
      xvals <- c(xvals, smooth$x[i])
      yvals <- c(yvals, smooth$y[i])
    }
  }
  grid.lines(x=xx, y=yy, gp=gpar(col=col.line, lty=lty, lwd=lwd, alpha=alpha),
               default.units="native", ...)
}

## A lattice-style panel function to draw polygons (like coverage)
## Arguments:
##    o x, y: the x and y coordinates form the plot
##    o lwd, lty, col: color, with and type of the plot lines
##    o fill: fill colors for areas above and under the baseline, a vector of length two
##    o col.line: color of the baseline
##    o baseline: the y value of the horizontal baseline
##    o alpha: the transparancy
## Value: the function is called for its side-effect of drawing on the graphics device
.panel.polygon <- function (x, y, lwd=plot.line$lwd, lty=plot.line$lty, col,
                            col.line=plot.line$col, baseline, fill, alpha=1, ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  fill <- rep(fill,2)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 1)
    return()
  if (!missing(col)) {
    if (missing(col.line))
      col.line <- col
  }
  x <- x[ok]
  y <- y[ok]
  plot.line <- trellis.par.get("plot.line")
  changePoint = NULL
  tmp <- as.integer(y<baseline)
  for(i in seq_along(tmp))
    if(i>1 && tmp[i]!= tmp[i-1])
      changePoint = c(changePoint, i)
  m <- (y[changePoint] - y[changePoint-1]) / (x[changePoint] - x[changePoint-1])
  xCross <- ((baseline-y[changePoint-1])/m) + x[changePoint-1]
  newX <- newY <- NULL
  j <- 1
  x <- c(x, tail(x,1))
  y <- c(y, baseline)
  xvals <- x[1]
  yvals <- baseline
  for(i in seq_along(x)) {
    if(i==length(x)) {
      xvals <- c(xvals, x[i])
      yvals <- c(yvals, baseline)
      fcol <- if(mean(yvals)<baseline) fill[1] else fill[2]
      panel.polygon(xvals, yvals, fill=fcol, col=fcol, border=fcol, alpha=alpha)
    } else if(i %in% changePoint) {
      xvals <- c(xvals, xCross[j])
      yvals <- c(yvals, baseline)
      fcol <- if(mean(yvals)<baseline) fill[1] else fill[2]
      panel.polygon(xvals, yvals, fill=fcol, col=fcol, border=fcol, alpha=alpha)
      xvals <- c(xCross[j], x[i])
      yvals <- c(baseline, y[i])
      j <- j+1
    } else {
      xvals <- c(xvals, x[i])
      yvals <- c(yvals, y[i])
    }
  }
  grid.lines(x=x, y=y, gp=gpar(col=col.line, lty=lty, lwd=lwd, alpha=alpha),
               default.units="native", ...)
}


## A lattice-style panel function to draw box and whisker plots with groups
## Arguments: see ? panel.bwplot for details
## Value: the function is called for its side-effect of drawing on the graphics device
.panel.bwplot <- function (x, y, box.ratio=1, box.width=box.ratio/(1+box.ratio), lwd, lty, fontsize,
                           pch, col, alpha, cex, font, fontfamily, fontface, fill, varwidth=FALSE, notch=FALSE,
                           notch.frac = 0.5, ..., levels.fos=sort(unique(x)),
                           stats=boxplot.stats, coef=1.5, do.out=TRUE)
{
  if (all(is.na(x) | is.na(y)))
    return()
  x <- as.numeric(x)
  y <- as.numeric(y)
  cur.limits <- current.panel.limits()
  xscale <- cur.limits$xlim
  yscale <- cur.limits$ylim
  if (!notch)
    notch.frac <- 0
  blist <- tapply(y, factor(x, levels = levels.fos), stats,
                  coef = coef, do.out = do.out)
  blist.stats <- t(sapply(blist, "[[", "stats"))
  blist.out <- lapply(blist, "[[", "out")
  blist.height <- box.width
  if (varwidth) {
    maxn <- max(table(x))
    blist.n <- sapply(blist, "[[", "n")
    blist.height <- sqrt(blist.n/maxn) * blist.height
  }
  blist.conf <- if (notch)
    sapply(blist, "[[", "conf")
  else t(blist.stats[, c(2, 4), drop = FALSE])
  ybnd <- cbind(blist.stats[, 3], blist.conf[2, ], blist.stats[, 4], blist.stats[, 4], blist.conf[2, ],
                blist.stats[,3], blist.conf[1, ], blist.stats[, 2], blist.stats[, 2], blist.conf[1, ], blist.stats[, 3])
  xleft <- levels.fos - blist.height/2
  xright <- levels.fos + blist.height/2
  xbnd <- cbind(xleft + notch.frac * blist.height/2, xleft,
                xleft, xright, xright, xright - notch.frac * blist.height/2,
                xright, xright, xleft, xleft, xleft + notch.frac *
                blist.height/2)
  xs <- matrix(NA_real_, nrow = nrow(xbnd) * 2, ncol = ncol(xbnd))
  ys <- matrix(NA_real_, nrow = nrow(xbnd) * 2, ncol = ncol(xbnd))
  xs[seq(along.with = levels.fos, by = 2), ] <- xbnd[seq(along.with = levels.fos),]
  ys[seq(along.with = levels.fos, by = 2), ] <- ybnd[seq(along.with = levels.fos),]
  panel.polygon(t(xs), t(ys), lwd=lwd, lty=lty, col = fill, alpha=alpha, border=col)
  panel.segments(rep(levels.fos, 2), c(blist.stats[, 2], blist.stats[, 4]), rep(levels.fos, 2),
                 c(blist.stats[,1], blist.stats[, 5]), col=col, alpha=alpha,
                 lwd=lwd, lty=lty)
  panel.segments(levels.fos - blist.height/2, c(blist.stats[, 1], blist.stats[, 5]), levels.fos + blist.height/2,
                 c(blist.stats[, 1], blist.stats[, 5]), col=col, alpha=alpha, lwd=lwd, lty=lty)
  if (all(pch == "|")) {
    mult <- if (notch)
      1 - notch.frac
    else 1
    panel.segments(levels.fos - mult * blist.height/2,
                   blist.stats[, 3], levels.fos + mult * blist.height/2,
                   blist.stats[, 3], lwd=lwd, lty=lty,
                   col=col, alpha = alpha)
  }
  else {
    panel.points(x = levels.fos, y = blist.stats[, 3],
                 pch = pch, col = col, alpha = alpha, cex = cex,
                 fontfamily = fontfamily, fontface = .chooseFace(fontface,
                                            font), fontsize = fontsize)
  }
  panel.points(x = rep(levels.fos, sapply(blist.out, length)),
               y = unlist(blist.out), pch=pch, col=col,
               alpha=alpha, cex=cex,
               fontfamily=fontfamily, fontface = .chooseFace(fontface,
                                        font), fontsize = fontsize)
}



## Check which parameters have already been set for a GdObject, and
## update all missing ones from the prototype of the current parent
## class.
## Arguments:
##    o x: an object inheriting from class GdObject
##    o class: the parent class from which to draw the missing parameters
## Value: The updated GdObject
.updatePars <- function(x, class){
    current <- getPar(x, hideInternal=FALSE)
    defaults <- getPar(getClass(class)@prototype@dp, hideInternal=FALSE)
    ## Check whether we need to adjust any of those defaults based on the selected scheme
    sid <- getOption("Gviz.scheme")
    scheme <- if(is.null(sid)) list() else .schemes[[sid]]
    schemePars <- if(is.null(scheme) || is.null(scheme[[class]])) list() else scheme[[class]]
    if(!is.list(schemePars))
        stop(sprintf("Corrupted parameter definition for class '%s' in scheme '%s'", class, sid))
    defaults[names(schemePars)] <- schemePars
    missing <- setdiff(names(defaults), names(current))
    if(is.null(getPar(x, ".__appliedScheme")) && length(schemePars)){
        defaults[[".__appliedScheme"]] <- if(is.null(sid)) NA else sid
        defaults[names(schemePars)] <- schemePars
        missing <- c(union(missing, names(schemePars)), ".__appliedScheme")
    }
    x <- setPar(x, defaults[missing], interactive=FALSE)
    return(x)
}



## The scheme settings registry. This is essentially just a nested list of display parameters, and the values in this list
## will be used to initialize the objects. Please note that a display parameter still needs to be defined in the class
## definition for this to work, and that all parameters that are not set explicitely in the scheme will be taken from
## the defaults in that class definition. Scheme parameters still override everything else, and even parameters that are
## not defined for the class will be added from the scheme settings.
.schemes <- new.env()
.schemes[["default"]] <- list(

                              GdObject=list(
                                            alpha=1,
                                            background.panel="transparent",
                                            background.title="lightgray",
                                            cex.axis=NULL,
                                            cex.title=NULL,
                                            cex=1,
                                            col.axis="white",
                                            col.frame="lightgray",
                                            col.grid=.DEFAULT_SHADED_COL,
                                            col.line=NULL,
                                            col.symbol=NULL,
                                            col.title="white",
                                            col=.DEFAULT_SYMBOL_COL,
                                            collapse=TRUE,
                                            fill=.DEFAULT_FILL_COL,
                                            fontcolor.title="white",
                                            fontcolor="black",
                                            fontface.title=2,
                                            fontface=1,
                                            fontfamily.title="sans",
                                            fontfamily="sans",
                                            fontsize=12,
                                            frame=FALSE,
                                            grid=FALSE,
                                            h=-1,
                                            lineheight=1,
                                            lty.grid="solid",
                                            lty="solid",
                                            lwd.border=1,
                                            lwd.grid=1,
                                            lwd=1,
                                            min.distance=1,
                                            min.height=3,
                                            min.width=1,
                                            rot.title=90,
                                            rotation=0,
                                            showAxis=TRUE,
                                            showTitle=TRUE,
                                            size=1,
                                            v=-1
                                            ),

                              StackedTrack=list(
                                                stackHeight=0.75,
                                                reverseStacking=FALSE
                                                ),


                              AnnotationTrack=list(
                                                   arrowHeadWidth=30,
                                                   arrowHeadMaxWidth=40,
                                                   cex.group=0.6,
                                                   cex=1,
                                                   col.line="darkgray",
                                                   col=.DEFAULT_LINE_COL,
                                                   featureAnnotation=NULL,
                                                   fill="lightblue",
                                                   fontcolor.group=.DEFAULT_SHADED_COL,
                                                   fontcolor.item="white",
                                                   fontface.group=2,
                                                   groupAnnotation=NULL,
                                                   lex=1,
                                                   lineheight=1,
                                                   lty="solid",
                                                   lwd=1,
                                                   mergeGroups=FALSE,
                                                   rotation=0,
                                                   shape="arrow",
                                                   showFeatureId=NULL,
                                                   showId=NULL,
                                                   showOverplotting=FALSE,
                                                   size=1
                                                ),

                              DetailsAnnotationTrack=list(
                                                          details.minWidth=100,
                                                          details.ratio=Inf,
                                                          details.size=0.5,
                                                          detailsBorder.col="darkgray",
                                                          detailsBorder.fill="transparent",
                                                          detailsBorder.lty="solid",
                                                          detailsBorder.lwd=1,
                                                          detailsConnector.cex=1,
                                                          detailsConnector.col="darkgray",
                                                          detailsConnector.lty="dashed",
                                                          detailsConnector.lwd=1,
                                                          detailsConnector.pch=20,
                                                          detailsFunArgs=list(),
                                                          groupDetails=FALSE),

                              GeneRegionTrack=list(
                                                   arrowHeadWidth=10,
                                                   arrowHeadMaxWidth=20,
                                                   col=.DEFAULT_LINE_COL,
                                                   collapseTranscripts=FALSE,
                                                   exonAnnotation=NULL,
                                                   fill="#FFD58A",
                                                   min.distance=0,
                                                   shape=c("smallArrow", "box"),
                                                   showExonId=NULL,
                                                   thinBoxFeature=.THIN_BOX_FEATURES,
                                                   transcriptAnnotation=NULL
                                                   ),

                              BiomartGeneRegionTrack=list(
                                                          C_segment="burlywood4",
                                                          D_segment="lightblue",
                                                          J_segment="dodgerblue2",
                                                          Mt_rRNA="yellow",
                                                          Mt_tRNA="darkgoldenrod",
                                                          Mt_tRNA_pseudogene="darkgoldenrod1",
                                                          V_segment="aquamarine",
                                                          miRNA="cornflowerblue",
                                                          miRNA_pseudogene="cornsilk",
                                                          misc_RNA="cornsilk3",
                                                          misc_RNA_pseudogene="cornsilk4",
                                                          protein_coding="orange",
                                                          pseudogene="brown1",
                                                          rRNA="darkolivegreen1",
                                                          rRNA_pseudogene="darkolivegreen" ,
                                                          retrotransposed="blueviolet",
                                                          scRNA="gold4",
                                                          scRNA_pseudogene="darkorange2",
                                                          snRNA="coral",
                                                          snRNA_pseudogene="coral3",
                                                          snoRNA="cyan",
                                                          snoRNA_pseudogene="cyan2",
                                                          tRNA_pseudogene="antiquewhite3",
                                                          utr3="orange",
                                                          utr5="orange"
                                                          ),

                              GenomeAxisTrack=list(
                                                   add35=FALSE,
                                                   add53=FALSE,
                                                   background.title="transparent",
                                                   cex.id=0.7,
                                                   cex=0.8,
                                                   col.id="white",
                                                   col.range="cornsilk4",
                                                   distFromAxis=1,
                                                   exponent=NULL,
                                                   fill.range="cornsilk3",
                                                   fontcolor="#808080",
                                                   fontsize=10,
                                                   labelPos="alternating",
                                                   littleTicks=FALSE,
                                                   lwd=2,
                                                   scale=NULL,
                                                   showId=FALSE,
                                                   showTitle=FALSE,
                                                   size=NULL,
                                                   col="darkgray"
                                                   ),

                              DataTrack=list(
                                             aggregateGroups=FALSE,
                                             aggregation="mean",
                                             amount=NULL,
                                             baseline=NULL,
                                             box.ratio=1,
                                             box.width=NULL,
                                             cex.legend=0.8,
                                             cex.sampleNames=NULL,
                                             cex=0.7,
                                             coef=1.5,
                                             col.baseline=NULL,
                                             col.histogram=.DEFAULT_SHADED_COL,
                                             col.horizon=NA,
                                             col.mountain=NULL,
                                             col.sampleNames="white",
                                             col=trellis.par.get("superpose.line")[["col"]],
                                             collapse=FALSE,
                                             degree=1,
                                             do.out=TRUE,
                                             evaluation=50,
                                             factor=0.5,
                                             family="symmetric",
                                             fill.histogram=NULL,
                                             fill.horizon=c("#B41414", "#E03231", "#F7A99C", "#9FC8DC", "#468CC8", "#0165B3"),
                                             fill.mountain=c("#CCFFFF", "#FFCCFF"),
                                             fontcolor.legend=.DEFAULT_SHADED_COL,
                                             gradient=brewer.pal(9, "Blues"),
                                             groups=NULL,
                                             horizon.origin=0,
                                             horizon.scale=NULL,
                                             jitter.x=FALSE,
                                             jitter.y=FALSE,
                                             levels.fos=NULL,
                                             lty.baseline=NULL,
                                             lty.mountain=NULL,
                                             lwd.baseline=NULL,
                                             lwd.mountain=NULL,
                                             min.distance=0,
                                             na.rm=FALSE,
                                             ncolor=100,
                                             notch.frac=0.5,
                                             notch=FALSE,
                                             pch=20,
                                             separator=0,
                                             showColorBar=TRUE,
                                             showSampleNames=FALSE,
                                             size=NULL,
                                             span=1/5,
                                             stackedBars=TRUE,
                                             stats=boxplot.stats,
                                             transformation=NULL,
                                             type="p",
                                             varwidth=FALSE,
                                             window=NULL,
                                             windowSize=NULL,
                                             ylim=NULL
                                             ),

                              IdeogramTrack=list(
                                                 background.title="transparent",
                                                 bevel=0.45,
                                                 cex.bands=0.7,
                                                 cex=0.8,
                                                 col="red",
                                                 fill="#FFE3E6",
                                                 fontcolor=.DEFAULT_SHADED_COL,
                                                 fontsize=10,
                                                 showBandId=FALSE,
                                                 showId=TRUE,
                                                 showTitle=FALSE,
                                                 size=NULL
                                                 ),

                              SequenceTrack=list(
                                                 add53=FALSE,
                                                 background.title="transparent",
                                                 cex=1,
                                                 col="darkgray",
                                                 complement=FALSE,
                                                 fontcolor=getBioColor("DNA_BASES_N"),
                                                 fontface=2,
                                                 fontsize=10,
                                                 lwd=2,
                                                 min.width=2,
                                                 noLetters=FALSE,
                                                 rotation=0,
                                                 showTitle=FALSE,
                                                 size=NULL
                                                 ),

                              HighlightTrack=list(
                                                  col="red",
                                                  fill="#FFE3E6"
                                                  )

                              )

.schemes[["default.old"]] <- list(
                                  AnnotationTrack=list(col="transparent"))

## A helper function to be called upon package load that tries to find a stored Gviz scheme in the working directory
.collectSchemes <- function(){
    try({
    	if(".GvizSchemes" %in% base::ls(globalenv(), all.names=TRUE)){
            schemes <- get(".GvizSchemes", globalenv())
            if(is.list(schemes) && !is.null(names(schemes)))
                lapply(names(schemes), function(x) addScheme(schemes[[x]], x))
    	}
    },silent=TRUE)
}


## Helper function to select fontfaces
.chooseFace <- function (fontface = NULL, font = 1)
{
    if (is.null(fontface))
        font
    else fontface
}


## Get a scheme
## Arguments:
##    o name: the name of the scheme to get. Defaults to the current one.
## Value: A list containing the scheme
getScheme <- function(name=getOption("Gviz.scheme")){
    s <- NULL
    if(!is.null(name))
        s <- .schemes[[name]]
    if(is.null(s))
        s <- list()
    return(s)
}

## Add a new scheme
## Arguments:
##    o scheme: the scheme to add
##    o name: the name of the scheme to add
addScheme <- function(scheme, name){
    .schemes[[name]] <- scheme
}




## Helper function to compute native coordinate equivalent to 1 pixel and resize all ranges in a GRanges object accordingly
## Arguments:
##    o r: object inheriting from class GdObject
##    o min.width: the minimal pixel width of a region
##    o diff: the pixel resolution
## Value: the updated GdObject
.resize <- function(r, min.width=2, diff=.pxResolution(coord="x"))
{
    if(min.width>0){
        minXDiff <- ceiling(min.width*diff)
        ## Extend all ranges to at least minXDiff
        xdiff <- width(r)
        xsel <- xdiff < minXDiff
        if(any(xsel))
        {
            rr <- if(is(r, "GRanges")) ranges(r) else r
            start(rr)[xsel] <- pmax(1, start(rr)[xsel]-(minXDiff-xdiff[xsel])/2)
            end(rr)[xsel] <- end(rr)[xsel]+(minXDiff-xdiff[xsel])/2
            if(is(r, "GRanges")) r@ranges <- rr else r <- rr
        }
    }
    return(r)
}

## Tables containing the UCSC to ENSEMBL genome mapping
.biomartCurrentVersionTable <- read.delim(system.file(file.path("extdata", "biomartVersionsNow.txt"), package="Gviz"), as.is=TRUE)
.biomartVersionTable <- read.delim(system.file(file.path("extdata", "biomartVersionsLatest.txt"), package="Gviz"), as.is=TRUE)

## Helper function to map between UCSC and ENSEMBl genome information
## Arguments:
##    o id: character scalar, a UCSC genome identifier
## Value: a list with ENSEMBL the genome information
.ucsc2Ensembl <- function(id){
   mt <- match(tolower(id), tolower(.biomartCurrentVersionTable$ucscId))
   val <- .biomartCurrentVersionTable[mt, ]
   if(is.na(mt)){
       mt <- match(tolower(id), tolower(.biomartVersionTable$ucscId))
       val <- .biomartVersionTable[mt, c("species", "value", "dataset", "ucscId", "speciesShort", "speciesLong", "date", "version")]
   }
   return(as.list(val))
}

## Helper function to get the ENSEMBL biomart given a UCSC identifier
## Arguments:
##    o genome: character scalar, a UCSC genome identifier
## Value: a biomaRt object
.getBiomart <- function(genome){
    map <- .ucsc2Ensembl(genome)
    if(map$date == "head"){
        bm <- useMart("ensembl", dataset=map$dataset)
        ds <- listDatasets(bm)
        mt <- ds[match(map$dataset, ds$dataset), "version"]
        if(is.na(mt)){
            stop(sprintf(paste("Gviz thinks that the UCSC genome identifier '%s' should map to the Biomart data set '%s' which is not correct.",
                               "\nPlease manually provide biomaRt object"), genome, map$dataset))
        }
        if(mt != map$value){
            stop(sprintf(paste("Gviz thinks that the UCSC genome identifier '%s' should map to the current Biomart head as '%s',",
                               "but its current version is '%s'.\nPlease manually provide biomaRt object"),
                         genome, map$value, mt))
        }
    }else{
        bm <- useMart(host=sprintf("%s.archive.ensembl.org", tolower(sub(".", "", map$date, fixed=TRUE))), biomart="ENSEMBL_MART_ENSEMBL", dataset=map$dataset)
        ds <- listDatasets(bm)
        mt <- ds[match(map$dataset, ds$dataset), "version"]
        if(is.na(mt)){
            stop(sprintf(paste("Gviz thinks that the UCSC genome identifier '%s' should map to the Biomart data set '%s' which is not correct.",
                               "\nPlease manually provide biomaRt object"), genome, map$dataset))
        }
        if(mt != map$value){
            stop(sprintf(paste("Gviz thinks that the UCSC genome identifier '%s' should map to Biomart archive %s (version %s) as '%s',",
                               "but its version is '%s'.\nPlease manually provide biomaRt object"),
                         genome, sub(".", " ", map$date, fixed=TRUE), map$version, map$value, mt))
        }
    }
    return(bm)
}


## Helper function to translate from a UCSC genome name to a Biomart data set. This also caches the mart
## object in order to speed up subsequent calls
## Arguments:
##    o genome: character giving the UCSC genome
## Value: A BiomaRt connection object
.genome2Dataset <- function(genome)
{
    map <- .ucsc2Ensembl(genome)
    if(is.na(map$date)){
        stop(sprintf("Unable to automatically determine Biomart data set for UCSC genome identifier '%s'.\nPlease manually provide biomaRt object", genome))
    }
    cenv <- environment()
    bm <- .doCache(paste(map$dataset, genome, sep="_"), expression(.getBiomart(genome)), .ensemblCache, cenv)
    return(bm)
}



## Return the plotting range for a GdObject, either from the contained ranges or from overrides.
## This function is vectorized and should also work for lists of GdObjects.
.defaultRange <- function(GdObject, from=NULL, to=NULL, extend.left=0, extend.right=0, factor=0.01, annotation=FALSE)
{
    if(!is.list(GdObject))
        GdObject <- list(GdObject)
    GdObject <- c(GdObject, unlist(lapply(GdObject, function(x) if(is(x, "HighlightTrack") || is(x, "OverlayTrack")) x@trackList else NULL)))
    if(!length(GdObject) || !all(sapply(GdObject, is, "GdObject")))
        stop("All items in the list must inherit from class 'GdObject'")
    GdObject <- GdObject[!sapply(GdObject, is, "OverlayTrack")]
    tfrom <- lapply(GdObject, function(x){tmp <- start(x); if(is(x, "RangeTrack")) tmp <- tmp[seqnames(x)==chromosome(x)]; tmp})
    tfrom <- if(is.null(unlist(tfrom))) Inf else min(sapply(tfrom[listLen(tfrom)>0], min))
    tto <- lapply(GdObject, function(x){tmp <- end(x); if(is(x, "RangeTrack")) tmp <- tmp[seqnames(x)==chromosome(x)]; tmp})
    tto <- if(is.null(unlist(tto))) Inf else max(sapply(tto[listLen(tto)>0], max))
    if((is.null(from) || is.null(to)) && ((is.infinite(tfrom) || is.infinite(tto)) || is(GdObject, "GenomeAxisTrack")))
        stop("Unable to automatically determine plotting ranges from the supplied track(s).\nPlease provide ",
             "range coordinates through the 'from' and 'to' arguments of the plotTracks function.")
    ## FIX the cases with identical "tfrom" and "tto" (one base-pair plotting) by adding +1 to "tto"
    if (tto == tfrom) {
        tto <- tto + 1
    }
    range <- extendrange(r=c(tfrom, tto), f=factor)
    range[1] <- max(1, range[1])
    wasNull <- rep(FALSE, 2)
    if(is.null(from)){
        wasNull[1] <- TRUE
        from <- range[1]
    }
    if(is.null(to)){
        to <- range[2]
        wasNull[2] <- TRUE
    }
    ## We may need some extra space for annotations
    if(annotation){
        rr <- unlist(lapply(GdObject, function(x){
            gr <- .dpOrDefault(x, ".__groupRanges")
            gw <- .dpOrDefault(x, ".__groupLabelWidths", data.frame(before=0, after=0))
            if(is.null(gr) || length(gr) == 0) NULL else c(min(start(gr) + gw$before), max(end(gr) - gw$after))
        }))
        if(!is.null(rr)){
            rr <- matrix(rr, ncol=2, byrow=TRUE)
            if(wasNull[1])
               from <- min(from, rr[,1])
            if(wasNull[2])
                to <- max(to, rr[,2])
        }
    }
    from <- if(extend.left != 0 && extend.left > -1 && extend.left < 1){
        from - (abs(diff(c(from, to))) * extend.left)
    }else{
        from - extend.left
    }
    to <- if(extend.right != 0 && extend.right > -1 && extend.right < 1){
        to + (abs(diff(c(from, to))) * extend.right)
    }else{
        to + extend.right
    }
    if(from > to){
        warning("'from' range can not be larger than 'to', reversing range coordinates")
        tto <- from
        from <- to
        to <- tto
    }
    return(c(from=as.integer(from), to=as.integer(to)))
}


## Figure out the colors to use for a DataTrack object from the supplied display parameters
.getPlottingFeatures <- function(GdObject)
{
    pch <- .dpOrDefault(GdObject, "pch", 20)
    lty <- .dpOrDefault(GdObject, "lty", 1)
    lwd <- .dpOrDefault(GdObject, "lwd", 1)
    cex <- .dpOrDefault(GdObject, "cex", 0.7)
    groups <- .dpOrDefault(GdObject, "groups")
    col <- .dpOrDefault(GdObject, "col", "#0080ff")
    if(is.null(groups)){ ## When there are no groups we force a single color for all lines and points
        col <- col[1]
        col.line <- .dpOrDefault(GdObject, "col.line", col)[1]
        col.symbol <- .dpOrDefault(GdObject, "col.symbol", col)[1]
        pch <- pch[1]
        lwd <- lwd[1]
        lty <- lty[1]
        cex <- cex[1]
    } else { ## Otherwise colors are being mapped to group factors
        if(!is.factor(groups))
            groups <- factor(groups)
        col <- .dpOrDefault(GdObject, "col", trellis.par.get("superpose.line")$col)
        col <- rep(col, nlevels(groups))
        col.line <- rep(.dpOrDefault(GdObject, "col.line", col), nlevels(groups))
        col.symbol <- rep(.dpOrDefault(GdObject, "col.symbol", col), nlevels(groups))
        lwd <- rep(lwd, nlevels(groups))
        lty <- rep(lty, nlevels(groups))
        pch <- rep(pch, nlevels(groups))
        cex <- rep(cex, nlevels(groups))
    }
    col.baseline <- .dpOrDefault(GdObject, "col.baseline", col)
    col.grid <- .dpOrDefault(GdObject, "col.grid", "#e6e6e6")[1]
    fill <- .dpOrDefault(GdObject, "fill", .DEFAULT_FILL_COL)[1]
    fill.histogram <- .dpOrDefault(GdObject, "fill.histogram", fill)[1]
    col.histogram  <- .dpOrDefault(GdObject, "col.histogram", .dpOrDefault(GdObject, "col", .DEFAULT_SHADED_COL))[1]
    lty.grid <- .dpOrDefault(GdObject, "lty.grid", 1)
    lwd.grid <- .dpOrDefault(GdObject, "lwd.grid", 1)
    return(list(col=col, col.line=col.line, col.symbol=col.symbol, col.baseline=col.baseline,
                col.grid=col.grid, col.histogram=col.histogram, fill=fill, fill.histogram=fill.histogram,
                lwd=lwd, lty=lty, pch=pch, cex=cex, lwd.grid=lwd.grid, lty.grid=lty.grid))
}


.legendInfo <- function()
{
    legInfo <- matrix(FALSE, ncol=7, nrow=17, dimnames=list(c("p", "b", "l", "a", "s", "S", "r", "h", "smooth",
                                                              "histogram", "boxplot", "heatmap", "gradient", "mountain", "g", "horizon","confint"),
                                                            c("lty", "lwd", "pch", "col", "cex", "col.lines", "col.symbol")))
    legInfo[2:9, c("lty", "lwd", "col.lines")] <- TRUE
    legInfo[1:2, c("pch", "cex", "col.symbol")] <- TRUE
    legInfo[1:12, "col"] <- TRUE
    return(legInfo)
}

## A helper function to get the currently active chromosomes, also if the track is one of the collection
## track classes
.recChromosome <- function(GdObject){
    chroms <- if(is(GdObject, "HighlightTrack") || is(GdObject, "OverlayTrack"))
        unlist(lapply(GdObject@trackList, .recChromosome)) else chromosome(GdObject)
    return(unique(chroms))
}


## Plot a list of GdObjects as individual tracks similar to the display on the UCSC genome browser
## Arguments:
##    o trackList: a list of GdObjects
##    o from, to: the plotting range, will be figured out automatically from the tracks if missing
##    o sized: a vector of relative vertical sizes, or NULL to auto-detect
##    o panel.only: don't draw track titles, useful to embed in a lattice-like function, this also implies add=TRUE
##    o extend.right, extend.left: extend the coordinates in 'from' and 'too'
##    o title.width: the expansion factor for the width of the title track
## Value: the function is called for its side-effect of drawing on the graphics device
plotTracks <- function(trackList, from=NULL, to=NULL, ..., sizes=NULL, panel.only=FALSE, extend.right=0,
                       extend.left=0, title.width=NULL, add=FALSE, main, cex.main=2, fontface.main=2,
                       col.main="black", margin=6, chromosome=NULL, innerMargin=3){
    ## If we have to open a new device for this but do not run through the whole function because of errors we want to
    ## clean up in the end
    done <- FALSE
    cdev <- dev.cur()
    on.exit(if(cdev==1 && !done) dev.off())
    ## We only need a new plot for regular calls to the function. Both add==TRUE and panel.only=TRUE will add to an existing grid plot
    if(!panel.only && !add)
        grid.newpage()
    if(!is.list(trackList))
        trackList <- list(trackList)
    ## All arguments in ... are considered to be additional display parameters and need to be attached to each item in the track list
    dps <- list(...)
    trackList <- lapply(trackList, function(x) {
        displayPars(x, recursive=TRUE) <- dps
        return(x)
    })

    ## OverlayTracks and HighlightTracks can be discarded if they are empty
    trackList <- trackList[!sapply(trackList, function(x) (is(x, "HighlightTrack") || is(x, "OverlayTrack")) && length(x) < 1)]
    isHt <- which(sapply(trackList, is, "HighlightTrack"))
    isOt <- which(sapply(trackList, is, "OverlayTrack"))
    ## A mix between forward and reverse strand tracks should trigger an alarm
    strds <- unique(.whichStrand(trackList))
    if(!is.null(strds) && length(strds) > 1)
        warning("Plotting a mixture of forward strand and reverse strand tracks.\n Are you sure this is correct?")
    ## We first run very general housekeeping tasks on the tracks for which we don't really need to know anything about device
    ## size, resolution or plotting ranges.
    ## Chromosomes should all be the same for all tracks, if not we will force them to be set to the first one that can be detected.
    ## If plotting ranges are supplied we can speed up a lot of the downstream operations by subsetting first.
    ## We may want to use alpha blending on those devices that support it, but also fall back to non-transparent colors without causing
    ## warnings.
    hasAlpha <- .supportsAlpha()
    chrms <- unique(unlist(lapply(trackList, .recChromosome)))
    if(is.null(chromosome)){
        chrms <- if(!is.null(chrms)) chrms[gsub("^chr", "", chrms)!="NA"] else chrms
        chromosome <- head(chrms, 1)
        if(length(chromosome)==0)
            chromosome <- "chrNA"
        if(!is.null(chrms) && length(unique(chrms))!=1)
            warning("The track chromosomes in 'trackList' differ. Setting all tracks to chromosome '", chromosome, "'", sep="")
    }
    if(!is.null(from) || !(is.null(to))){
        trackList <- lapply(trackList, function(x){
            chromosome(x) <- chromosome
            subset(x, from=from, to=to, chromosome=chromosome, sort=FALSE, stacks=FALSE, use.defaults=FALSE)
        })
    }
    trackList <- lapply(trackList, consolidateTrack, chromosome=chromosome, any(.needsAxis(trackList)), any(.needsTitle(trackList)),
                                title.width, alpha=hasAlpha, ...)

    ## Now we figure out the plotting ranges. If no ranges are given as function arguments we take the absolute min/max of all tracks.
    ranges <- .defaultRange(trackList, from=from, to=to, extend.left=extend.left, extend.right=extend.right, annotation=TRUE)
    ## Now we can subset all the objects in the list to the current boundaries and compute the initial stacking
    trackList <- lapply(trackList, subset, from=ranges["from"], to=ranges["to"], chromosome=chromosome)
    trackList <- lapply(trackList, setStacks, recomputeRanges=FALSE)
    ## Highlight tracks are just a way to add a common highlighting region to several tracks, but other than that we can treat the containing
    ## tracks a normal track objects, and thus unlist them. We only want to record their indexes in the expanded list for later.
    htList <- list()
    expandedTrackList <- if(length(isHt)){
        j <- 1
        tlTemp <- list()
        for(i in seq_along(trackList)){
            if(! i %in% isHt){
                tlTemp <- c(tlTemp, trackList[[i]])
                j <- j+1
            }else{
                tlTemp <- c(tlTemp, trackList[[i]]@trackList)
                htList[[as.character(i)]] <- list(indexes=j:(j+length(trackList[[i]]@trackList)-1),
                                                  track=trackList[[i]])
                j <- j+length(trackList[[i]]@trackList)
            }
        }
        tlTemp
    }else trackList
    ## If there is a AlignmentsTrack and also a SequenceTrack we can tell the former to use the latter, unless already provided
    isAt <- sapply(expandedTrackList, is, "AlignmentsTrack")
    isSt <- sapply(expandedTrackList, is, "SequenceTrack")
    for(ai in which(isAt)){
        if(is.null(expandedTrackList[[ai]]@referenceSequence) && any(isSt))
            expandedTrackList[[ai]]@referenceSequence <- expandedTrackList[[min(which(isSt))]]
    }
    ## We need to reverse the list to get a top to bottom plotting order
    expandedTrackList <- rev(expandedTrackList)
    map <- vector(mode="list", length=length(expandedTrackList))
    titleCoords <- NULL
    names(map) <- rev(sapply(expandedTrackList, names))
    ## Open a fresh page and set up the bounding box, unless add==TRUE
    if(!panel.only) {
        ## We want a margin pixel border
        borderFacts <- 1-((margin*2)/vpLocation()$size)
        vpBound <- viewport(width=borderFacts[1], height=borderFacts[2])
        pushViewport(vpBound)
        ## If there is a header we have to make some room for it here
        if(!missing(main) && main != "")
        {
            vpHeader <- viewport(width=1, height=0.1, y=1, just=c("center", "top"))
            pushViewport(vpHeader)
            grid.text(main, gp=gpar(col=col.main, cex=cex.main, fontface=fontface.main))
            popViewport(1)
            vpMain <-  viewport(width=1, height=0.9, y=0.9, just=c("center", "top"))
        }else{
            vpMain <-  viewport(width=1, height=1)
        }
        pushViewport(vpMain)
        ## A first guestimate of the vertical space that's needed
        spaceSetup <- .setupTextSize(expandedTrackList, sizes, title.width, spacing=innerMargin)
    } else {
        vpBound <- viewport()
        pushViewport(vpBound)
        spaceSetup <- .setupTextSize(expandedTrackList, sizes, spacing=innerMargin)
    }
    ## First iteration to set up all the dimensions by calling the drawGD methods in prepare mode, i.e.,
    ## argument prepare=TRUE. Nothing is drawn at this point, and this only exists to circumvent the
    ## chicken and egg problem of not knowing how much space we need until we draw, but also not knowing
    ## where to draw until we know the space needed.
    for(i in rev(seq_along(expandedTrackList)))
    {
        fontSettings <- .fontGp(expandedTrackList[[i]], cex=NULL)
        vpTrack <-  viewport(x=0, y=sum(spaceSetup$spaceNeeded[1:i]), just=c(0,1), width=1, height=spaceSetup$spaceNeeded[i],
                             gp=fontSettings)
        pushViewport(vpTrack)
        vpContent <- if(!panel.only) viewport(x=spaceSetup$title.width + spaceSetup$spacing,
                                              width=1 - spaceSetup$title.width - spaceSetup$spacing * 2, just=0) else viewport(width=1)
        pushViewport(vpContent)
        expandedTrackList[[i]] <- drawGD(expandedTrackList[[i]], minBase=ranges["from"], maxBase=ranges["to"], prepare=TRUE, subset=FALSE)
        popViewport(2)
    }
    ## Now lets recalculate the space and draw for real
    spaceSetup <- .setupTextSize(expandedTrackList, sizes, title.width, spacing=innerMargin)
    ## First the highlight box backgrounds
    htBoxes <- data.frame(stringsAsFactors=FALSE)
    for(hlite in htList){
        inds <- setdiff(sort(length(expandedTrackList)-hlite$index+1), which(sapply(expandedTrackList, is, "IdeogramTrack")))
        y <- reduce(IRanges(start=inds, width=1))
        yy <- ifelse(start(y)==1, 0, sum(spaceSetup$spaceNeeded[1:start(y)-1]))
        ht <- sum(spaceSetup$spaceNeeded[start(y):end(y)])
        htBoxes <- rbind(htBoxes, data.frame(y=yy, height=ht, x=start(hlite$track), width=width(hlite$track),
                                             col=.dpOrDefault(hlite$track, "col", "orange"),
                                             fill=.dpOrDefault(hlite$track, "fill", "red"),
                                             lwd=.dpOrDefault(hlite$track, "lwd", 1),
                                             lty=.dpOrDefault(hlite$track, "lty", 1),
                                             alpha=.dpOrDefault(hlite$track, "alpha", 1),
                                             inBackground=.dpOrDefault(hlite$track, "inBackground", TRUE),
                                             stringsAsFactors=FALSE))
    }
    .drawHtBoxes <- function(htBoxes, background=TRUE){
        htBoxes <- htBoxes[htBoxes$inBackground == background, , drop=FALSE]
        if(nrow(htBoxes)){
            vpContent <- if(!panel.only) viewport(x=spaceSetup$title.width + spaceSetup$spacing, xscale=ranges,
                                                  width=1 - spaceSetup$title.width - spaceSetup$spacing*2, just=0) else viewport(width=1, xscale=ranges)
            pushViewport(vpContent)
            grid.rect(x=htBoxes$x, just=c(0,1), width=htBoxes$width, y=htBoxes$y+htBoxes$height, height=htBoxes$height,
                      gp=gpar(col=htBoxes$col, fill=htBoxes$fill, lwd=htBoxes$lwd, lty=htBoxes$lty, alpha=htBoxes$alpha), default.units="native")
            popViewport(1)
        }
    }
    if(nrow(htBoxes))
        .drawHtBoxes(htBoxes)
    ## Now the track content
    for(i in rev(seq_along(expandedTrackList)))
    {
        vpTrack <-  viewport(x=0, y=sum(spaceSetup$spaceNeeded[1:i]), just=c(0,1), width=1, height=spaceSetup$spaceNeeded[i])
        pushViewport(vpTrack)
		fill <- .dpOrDefault(expandedTrackList[[i]], "background.title", .DEFAULT_SHADED_COL)
        thisTrack <- if(is(expandedTrackList[[i]], "OverlayTrack")) expandedTrackList[[i]]@trackList[[1]] else expandedTrackList[[i]]
        if(!panel.only) {
            fontSettings <- .fontGp(expandedTrackList[[i]], subtype="title", cex=NULL)
            vpTitle <- viewport(x=0, width=spaceSetup$title.width, just=0, gp=fontSettings)
            pushViewport(vpTitle)
            lwd.border.title <- .dpOrDefault(thisTrack, "lwd.title", 1)
            col.border.title <- .dpOrDefault(thisTrack, "col.border.title", "transparent")
            grid.rect(gp=gpar(fill=fill, col=col.border.title, lwd=lwd.border.title))
            needAxis <- .needsAxis(thisTrack)
            drawAxis(thisTrack, ranges["from"], ranges["to"], subset=FALSE)
            tit <- spaceSetup$nwrap[i]
            ## FIXME: Do we want something smarted for the image map coordinates?
            titleCoords <- rbind(titleCoords, cbind(.getImageMap(cbind(0,0,1,1)),
                                                    title=names(expandedTrackList[[i]])))
            if(.dpOrDefault(thisTrack, "showTitle", TRUE) && !is.null(tit) && tit!="")
            {
                x <- if(needAxis) 0.075 else 0.4
                just <- if(needAxis) c("center", "top") else "center"
                ## FIXME: We need to deal with this when calculating the space for the title bar
                rot <- .dpOrDefault(thisTrack, "rotation.title", 90)
                gp <- .fontGp(thisTrack, "title", cex=spaceSetup$cex[i])
                suppressWarnings(grid.text(tit, unit(x, "npc"), rot=rot, gp=gp, just=just))
            }
            popViewport(1)
        }
        ## Draw the panel background, grid lines if necessary and the panel content
        vpBackground <- if(!panel.only) viewport(x=spaceSetup$title.width,
                                                 width=1-spaceSetup$title.width, just=0) else viewport(width=1)
        pushViewport(vpBackground)
        grid.rect(gp=gpar(col="transparent", fill=.dpOrDefault(thisTrack, "background.panel", "transparent")))
        drawGrid(thisTrack, ranges["from"], ranges["to"])
        popViewport(1)
        fontSettings <- .fontGp(expandedTrackList[[i]], cex=NULL)
        vpContentOuter <- if(!panel.only)
            viewport(x=spaceSetup$title.width, width=1-spaceSetup$title.width,
                     just=0, gp=fontSettings, clip=TRUE) else viewport(width=1, gp=fontSettings, clip=TRUE)
        pushViewport(vpContentOuter)
        vpContent <- if(!panel.only)
            viewport(x=spaceSetup$spacing, width=1-(spaceSetup$spacing*2), just=0, gp=fontSettings) else viewport(width=1, gp=fontSettings)
        pushViewport(vpContent)
        tmp <- drawGD(expandedTrackList[[i]], minBase=ranges["from"], maxBase=ranges["to"], subset=FALSE)
        if(!is.null(tmp))
            map[[(length(map)+1)-i]] <- tmp
        popViewport(2)
        if(.dpOrDefault(thisTrack, "frame", FALSE))
            grid.rect(gp=gpar(col=.dpOrDefault(thisTrack, "col.frame", .DEFAULT_SHADED_COL), fill="transparent"))
        popViewport(1)
    }
    if(nrow(htBoxes))
        .drawHtBoxes(htBoxes, FALSE)
    popViewport(if(panel.only) 1 else 2)
    tc <- as.character(titleCoords[,5])
    tc[which(tc == "" | is.na(tc) | is.null(tc))] = "NA"
    names(tc) <- tc
    if(!is.null(titleCoords))
    {
        tcoord <- as.matrix(titleCoords[,1:4])
        rownames(tcoord) <- names(tc)
        map$titles <- ImageMap(coords=tcoord, tags=list(title=tc))
    }
    done <- TRUE
    return(invisible(map))
}



## Try to extract the (unique) genome information from a GRanges objects with the possibility to fall back to a default value
.getGenomeFromGRange <- function(range, default=NULL){
    gn <- genome(range)
    if(length(unique(gn))>1)
        warning("Only a single genome is supported for this object. Ignoring additional genome information")
    if(length(gn)==0 || all(is.na(gn))){
        if(is.null(default))
            stop("A genome must be supplied when creating this object.")
        return(default[1])
    }
    return(gn[1])
}










## Write all tracks in a list of tracks into
## a single BED file.
exportTracks <- function(tracks, range, chromosome, file)
{
    if(missing(file))
        file <- "customTracks.bed"
    con <- file(file, open="wt")
    writeLines(sprintf("browser position %s:%i-%i", chromosome, range[1], range[2]),
               con)
    writeLines("browser hide all", con)
    for(t in seq_along(tracks))
    {
        track <- tracks[[t]]
        if(length(track)>0 && (is(track, "AnnotationTrack") || is(track, "GeneRegion")))
        {
            track <- as(track, "UCSCData")
            writeLines(as(track@trackLine, "character"), con)
            ## nextMet <- selectMethod("export.bed", c("RangedData", "characterORconnection"))
            ## nextMet(as(track, "RangedData"), con)
            .expBed(as(track, "RangedData"), con)
        }
    }
    close(con)
}

## This funcion is broken in the rtracklayer package
.expBed <- function (object, con, variant = c("base", "bedGraph", "bed15"), color, append)
{
    variant <- match.arg(variant)
    name <- strand <- thickStart <- thickEnd <- color <- NULL
    blockCount <- blockSizes <- blockStarts <- NULL
    df <- data.frame(chrom(object), start(object) - 1, end(object))
    score <- score(object)
    if (!is.null(score)) {
        if (!is.numeric(score) || any(is.na(score)))
            stop("Scores must be non-NA numeric values")
    }
    if (variant == "bedGraph") {
        if (is.null(score))
            score <- 0
        df$score <- score
    }
    else {
        blockSizes <- object$blockSizes
            blockStarts <- object$blockStarts
        if (variant == "bed15" && is.null(blockSizes))
            blockStarts <- blockSizes <- ""
        if (!is.null(blockSizes) || !is.null(blockStarts)) {
            if (is.null(blockSizes))
                stop("'blockStarts' specified without 'blockSizes'")
            if (is.null(blockStarts))
                stop("'blockSizes' specified without 'blockStarts'")
            lastBlock <- function(x) sub(".*,", "", x)
            lastSize <- lastBlock(blockSizes)
            lastStart <- lastBlock(blockStarts)
            if (any(df[[2]] + as.integer(lastSize) + as.integer(lastStart) != df[[3]]) ||
                any(sub(",.*", "", blockStarts) != 0))
                stop("blocks must span entire feature")
            blockCount <- sapply(strsplit(blockSizes, ","), length)
        }
        if (is.null(color))
            color <- object$itemRgb
        if (is.null(color) && !is.null(blockCount))
            color <- "0"
        else if (!is.null(color)) {
            nacol <- is.na(color)
            colmat <- col2rgb(color)
            color <- paste(colmat[1, ], colmat[2, ], colmat[3,
                                                            ], sep = ",")
            color[nacol] <- "0"
        }
        thickStart <- object$thickStart
        thickEnd <- object$thickEnd
        if (is.null(thickStart) && !is.null(color)) {
            thickStart <- start(object)
            thickEnd <- end(object)
        }
        strand <- object$strand
        if (!is.null(thickStart) && is.null(strand)) {
            strand <- rep(NA, nrow(object))
        }
        if (!is.null(strand) && is.null(score))
            score <- 0
        name <- object$name
        if (is.null(name))
            name <- rownames(object)
        if (!is.null(score) && is.null(name))
            name <- rep(NA, nrow(object))
        df$name <- name
        df$score <- score
        df$strand <- strand
        df$thickStart <- thickStart
        df$thickEnd <- thickEnd
        df$itemRgb <- color
        df$blockCount <- blockCount
        df$blockSizes <- blockSizes
        df$blockStarts <- blockStarts
        if (variant == "bed15") {
            df$expCount <- object$expCount
            df$expIds <- object$expIds
            df$expScores <- object$expScores
        }
    }
    scipen <- getOption("scipen")
    options(scipen = 100)
    on.exit(options(scipen = scipen))
    write.table(df, con, sep = "\t", col.names = FALSE, row.names = FALSE,
                quote = FALSE, na = ".", append = append)
}





## Construct a URL to UCSC showing the custom tracks
ucscUrl <- function(chr, range, spec, gen, open=TRUE)
{
    hgid <- system(sprintf("%s %s %s", system.file("lib/testUCSC.pl", package="Gviz"),
                           "customTracks.bed", spec, gen), intern=TRUE, ignore.stderr=TRUE)

    url <- sprintf(paste("http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=%s&Submit=go+to+genome+browser",
                         "&position=%s%%3A%i-%i", sep=""), hgid, chr, range[1], range[2])
    if(open)
        browseURL(url)
    return(url)
}




.updateObj <- function(object)
{
    availSlots <- getObjectSlots(object)
    availSlotNames <- names(availSlots)
    definedSlotNames <- slotNames(object)
    if(length(availSlotNames)==length(definedSlotNames) && all(sort(availSlotNames) == sort(definedSlotNames)))
        return(object)
    commonSlots <- intersect(definedSlotNames, availSlotNames)
    missingSlots <- setdiff(definedSlotNames, availSlotNames)
    newObject <- new(class(object))
    for (s in commonSlots)
        slot(newObject, s) <- availSlots[[s]]
    return(newObject)
}



vpLocation <- function(){
  xres <- devRes()[1]
  yres <- devRes()[2]
  ## find location and pixel-size of current viewport
  devloc1 <- c(convertX(unit(0, "npc"), "inches"),
              convertY(unit(0, "npc"), "inches"), 1)  %*% current.transform()
  devloc2 <- c(convertX(unit(1, "npc"), "inches"),
              convertY(unit(1, "npc"), "inches"), 1)  %*% current.transform()
  x1 <- (devloc1/devloc1[3])[1]*xres
  y1 <- (devloc1/devloc1[3])[2]*yres
  x2 <- (devloc2/devloc2[3])[1]*xres
  y2 <- (devloc2/devloc2[3])[2]*yres
  loc <- c(x1,y1,x2,y2)
  names(loc) <- c("x1", "y1", "x2", "y2")
  size <- c(x2-x1, y2-y1)
  names(size) <- c("width", "height")
  iloc <- c(x1/xres, y1/yres, x2/yres, y2/yres)
  names(iloc) <- c("x1", "y1", "x2", "y2")
  isize <- size/c(xres,yres)
  names(size) <- c("width", "height")
  return(list(location=loc, size=size, ilocation=iloc,
              isize=isize))
}



devRes <- function(){
  ## find R's resolution for the current device
  if(current.viewport()$name != "ROOT"){
    vpt <- current.vpTree()
    depth <- upViewport(0)
    xres <- abs(as.numeric(convertWidth(unit(1, "inches"), "native")))
    yres <- abs(as.numeric(convertHeight(unit(1, "inches"), "native")))
    downViewport(depth)
  }else{
    xres <- abs(as.numeric(convertWidth(unit(1, "inches"), "native")))
    yres <- abs(as.numeric(convertHeight(unit(1, "inches"), "native")))
  }
  retval <- c(xres, yres)
  names(retval) <- c("xres", "yres")
  return(retval)
}



devDims <- function(width, height, ncol=12, nrow=8, res=72){
 f <- (((ncol+1)*0.1+ncol+1)/((nrow+1)*0.1+nrow+1))
 if((missing(width) & missing(height) || !missing(width) & !missing(height)))
   stop("Need either argument 'width' or argument 'height'")
 if(missing(height))
   return(list(width=width, height=width/f, pwidth=width*res, pheight=width/f*res))
 else
   return(list(width=height*f, height, pwidth=height*f*res, pheight=height*res))
}

## Record the display parameters for each class once
.makeParMapping <- function()
{
    classes <-  c("GdObject", "GenomeAxisTrack", "RangeTrack", "NumericTrack", "DataTrack", "IdeogramTrack", "StackedTrack",
                  "AnnotationTrack", "DetailsAnnotationTrack", "GeneRegionTrack", "BiomartGeneRegionTrack", "AlignedReadTrack")
    defs <- try(sapply(classes, function(x) as(getClassDef(x)@prototype@dp, "list"), simplify=FALSE), silent=TRUE)
    if(!is(defs, "try-error") && is.null(.parMappings))
        assignInNamespace(x=".parMappings", value=defs, ns="Gviz")
}
.parMappings <- NULL

## Show available display parameters for a class and their defaults
availableDisplayPars <- function(class)
{
    if(!is.character(class))
        class <- class(class)
    class <- match.arg(class, c("GdObject", "GenomeAxisTrack", "RangeTrack", "NumericTrack", "DataTrack", "IdeogramTrack", "StackedTrack",
                                "AnnotationTrack", "DetailsAnnotationTrack", "GeneRegionTrack", "BiomartGeneRegionTrack", "AlignedReadTrack",
                                "SequenceTrack", "SequenceBSgenomeTrack", "SequenceDNSStringSetTrack"))
    parents <- names(getClassDef(class)@contains)
    .makeParMapping()
    pars <- .parMappings[c(parents, class)]
    finalPars <- inherited <- list()
    for(p in names(pars))
    {
        finalPars[names(pars[[p]])] <- pars[[p]]
        inherited[names(pars[[p]])] <- p
    }
    finalPars <- finalPars[order(names(finalPars))]
    inherited <- inherited[order(names(inherited))]
    return(new("InferredDisplayPars", name=class, inheritance=unlist(inherited), finalPars))
}

## Compute ellipse outline coordinates for bounding boxes
.box2Ellipse <- function(box, np=50)
{
    t <- seq(0, 2*pi, len=np)
    box$width <- box$cx2-box$cx1
    box$height <- box$cy2-box$cy1
    x <- rep(box$cx2 + -box$width + box$width/2, each=np)
    y <- rep(box$cy1 + box$height/2, each=np)
    a <- rep(box$width/2, each=np)
    b <- rep(box$height/2, each=np)
    tau <- 0
    xt <- x + (a*cos(t)*cos(tau)-b*sin(t)*sin(tau))
    yt <- y + (a*cos(t)*sin(tau)+b*sin(t)*cos(tau))
    return(data.frame(x1=xt, y1=yt, id=rep(seq_len(nrow(box)), each=np)))
}

## We store some preset in the options on package load
.onLoad = function(...){
    .collectSchemes()
    options("ucscChromosomeNames"=TRUE, "Gviz.scheme"="default", "Gviz.ucscUrl"=NULL)
}


## A helper function to replace missing function arguments in a list with NULL values. The function environment needs
## to be passed in as argument 'env' for this to work.
.missingToNull <- function(symbol, env=parent.frame()){
    for(i in symbol){
        mis <- try(do.call(missing, args=list(i), envir=env), silent=TRUE)
        if(!is(mis, "try-error") && mis)
            assign(i, NULL, env)
    }
}


## build a covariates data.frame from a variety of different inputs
.getCovars <- function(x){
    if(is.data.frame(x)){
        x
    }else{
        if(is(x, "GRanges")){
            as.data.frame(mcols(x))
        }else{
            if(is(x, "GRangesList")){
                as.data.frame(mcols(unlist(x)))
            }else{
                data.frame()
                ## stop(sprintf("Don't know how to extract covariates from a %s object", class(x)))
            }
        }
    }
}


## Prepare a data.frame or matrix containing the data for a DataTrack object. This involves trying to coerce
## and dropping non-numeric columns with a warning
.prepareDtData <- function(data, len=0){
    if(ncol(data) && nrow(data)){
        for(i in seq_along(data)){
            if(is.character(data[,i]))
                data[,i] <- type.convert(data[,i], as.is=TRUE)
        }
        isNum <- sapply(data, is.numeric)
        if(any(!isNum))
            warning(sprintf("The following non-numeric data column%s been dropped: %s", ifelse(sum(!isNum)>1, "s have", " has"),
                            paste(colnames(data)[!isNum], collapse=", ")))
        if(sum(dim(data))>0){
            data <- t(data[,isNum, drop=FALSE])
        }
    }else{
        data <- matrix(ncol=len, nrow=0)
    }
    if(all(is.na(data)))
        data <- matrix(ncol=len, nrow=0)
    if(ncol(data) != len)
        stop("The columns in the 'data' matrix must match the genomic regions.")
    return(data)
}



## An import function for gff3 files that tries to resolve the parent-child relationship
## between genes, transcripts and exons
.import.gff3 <- function(file){
    dat <- import.gff3(file)
    res <- try({
        genes <- tolower(dat$type) == "gene"
        ginfo <- mcols(dat[genes, ])
        dat <- dat[!genes]
        transcripts <- tolower(dat$type) == "mrna"
        tinfo <- mcols(dat[transcripts, ])
        dat <- dat[!transcripts]
        mt <- match(as.character(dat$Parent), as.character(tinfo$ID))
        if(!all(is.na(mt))){
            if(!"transcript_id" %in% colnames(mcols(dat))) {
                tid <- rep(NA, length(mt))
                tid[!is.na(mt)] <- tinfo[mt[!is.na(mt)], "ID"]
                mcols(dat)[["transcript_id"]] <- tid
            }
            if(!"transcript_name" %in% colnames(mcols(dat))) {
                tn <- rep(NA, length(mt))
                tn[!is.na(mt)] <- tinfo[mt[!is.na(mt)], "Name"]
                mcols(dat)[["transcript_name"]] <- tn
            }
            mt2 <- rep(NA, dim(mcols(dat))[1])
            mt2[!is.na(mt)] <- match(as.character(tinfo[mt[!is.na(mt)], "Parent"]), as.character(ginfo$ID))

            if(!all(is.na(mt2))){
                if(!"gene_id" %in% colnames(mcols(dat))) {
                    gid <- rep(NA, length(mt2))
                    gid[!is.na(mt2)] <- ginfo[mt[!is.na(mt2)], "ID"]
                    mcols(dat)[["gene_id"]] <- gid
                }
                if(!"gene_name" %in% colnames(mcols(dat))) {
                    gn <- rep(NA, length(mt2))
                    gn[!is.na(mt2)] <- ginfo[mt[!is.na(mt2)], "Name"]
                    mcols(dat)[["gene_name"]] <- gn
                }
            }
        }
        if(all(is.na(mcols(dat)[["ID"]])))
            mcols(dat)[["ID"]] <- paste("item", seq_along(dat), sep="_")
        if(!"exon_id" %in% colnames(mcols(dat)))
            mcols(dat)[["exon_id"]] <- mcols(dat)[["ID"]]
        if(!is.null(mcols(dat)[["gene_name"]]) && all(is.na(mcols(dat)[["gene_name"]])))
            mcols(dat)[["gene_name"]] <- NULL
        if(all(is.na(mcols(dat)[["transcript_name"]])))
            mcols(dat)[["transcript_name"]] <- NULL
        dat
    })
    if(is(res, "try-error")){
        warning(sprintf(paste("File '%s' is not valid according to the GFF3 standard and can not be properly parsed.",
                              "Results may not be what you expected!"), file))
        res <- dat
    }
    return(res)
}

## An import function for bigWig files that knowns how to deal with missing seqnames
.import.bw <- function(file, selection){
    bwf <- BigWigFile(path.expand(file))
    if(missing(selection)){
        rr <- import.bw(con=bwf)
    }else{
        si <- seqinfo(bwf)
        rr <- if(!as.character(seqnames(selection)[1]) %in% seqnames(seqinfo(bwf))){
            GRanges(seqnames(selection)[1], ranges=IRanges(1,2), score=1)[0] }else{
                import.bw(con=bwf, selection=selection)}
    }
    return(rr)
}

## An import function for bam files that distinguishes between DataTracks and AnnotationTracks
## FIXME: We probably want this to be able to deal with Gapped Alignments...
.import.bam <- function(file, selection){
    if(!file.exists(paste(file, "bai", sep=".")) &&
       !file.exists(paste(paste(head(strsplit("xxx.bam", ".", fixed=TRUE)[[1]], -1), collapse="."), "bai", sep=".")))
        stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
             "library(Rsamtools)\n\tindexBam(\"", file, "\")")
    sinfo <- scanBamHeader(file)[[1]]
    if(parent.env(environment())[["._trackType"]] == "DataTrack"){
        res <- if(!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)){
            mcols(selection) <- DataFrame(score=0)
            selection
        }else{
            param <- ScanBamParam(what=c("pos", "qwidth"), which=selection, flag=scanBamFlag(isUnmappedQuery=FALSE))
            x <- scanBam(file, param=param)[[1]]
            cov <- coverage(IRanges(x[["pos"]], width=x[["qwidth"]]))
            if(length(cov)==0){
                mcols(selection) <- DataFrame(score=0)
                selection
            }else{
                GRanges(seqnames=seqnames(selection), ranges=IRanges(start=start(cov), end=end(cov)), strand="*", score=runValue(cov))
            }
        }
    } else {
        res <- if(!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)){
            mcols(selection) <- DataFrame(id="NA", group="NA")
            selection[0]
        }else{
            param <- ScanBamParam(what=c("pos", "qwidth", "strand", "qname"), which=selection, flag=scanBamFlag(isUnmappedQuery=FALSE))
            x <- scanBam(file, param=param)[[1]]
            GRanges(seqnames=seqnames(selection), ranges=IRanges(start=x[["pos"]], width=x[["qwidth"]]), strand=x[["strand"]],
                    id=make.unique(x[["qname"]]), group=x[["qname"]])
        }
    }
    return(res)
}

.import.bam.alignments <- function(file, selection){
    indNames <- c(sub("\\.bam$", ".bai", file), paste(file, "bai", sep="."))
    index <- NULL
    for(i in indNames){
        if(file.exists(i)){
            index <- i
            break
        }
    }
    if(is.null(index))
        stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
             "library(Rsamtools)\n\tindexBam(\"", file, "\")")
    pairedEnd <- parent.env(environment())[["._isPaired"]]
    if(is.null(pairedEnd))
        pairedEnd <- TRUE
    bf <- BamFile(file, index=index, asMates=pairedEnd)
    param <- ScanBamParam(which=selection, what=scanBamWhat(), tag="MD", flag=scanBamFlag(isUnmappedQuery=FALSE))
    reads <- if(as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param=param)[[1]] else list()
    md <- if(is.null(reads$tag$MD)) rep(as.character(NA), length(reads$pos)) else reads$tag$MD
    if(length(reads$pos)){
        layed_seq <- sequenceLayer(reads$seq, reads$cigar)
        region <- unlist(bamWhich(param), use.names=FALSE)
        ans <- stackStrings(layed_seq, start(region), end(region), shift=reads$pos-1L, Lpadding.letter="+", Rpadding.letter="+")
        names(ans) <- seq_along(reads$qname)
    }else{
        ans <- DNAStringSet()
    }
    return(GRanges(seqnames=if(is.null(reads$rname)) character() else reads$rname,
                   strand=if(is.null(reads$strand)) character() else reads$strand,
                   ranges=IRanges(start=reads$pos, width=reads$qwidth),
                   id=reads$qname, cigar=reads$cigar, mapq=reads$mapq, flag=reads$flag, md=md, seq=ans,
                   isize=reads$isize, groupid=if(pairedEnd) reads$groupid else seq_along(reads$pos),
                   status=if(pairedEnd) reads$mate_status else rep(factor("unmated", levels=c("mated", "ambiguous", "unmated")),
                                                                   length(reads$pos))))
}


## An import function for fasta file that supports streaming if an index is present
.import.fasta <- function(file, selection, strict=TRUE){
    ffile <- FastaFile(file)
    if(!file.exists(paste(file, "fai", sep="."))){
        if(strict){
            stop("Unable to find index for fasta file '", file, "'. You can build an index using the following command:\n\t",
                 "library(Rsamtools)\n\tindexFa(\"", file, "\")")
        }else{
            return(readDNAStringSet(file))
        }
    }
    idx <- scanFaIndex(file)
    if(!as.character(seqnames(selection)[1]) %in% as.character(seqnames(idx))){
         return(DNAStringSet())
    }else{
        return(scanFa(file, selection))
    }
}


## An import function for the indexed 2bit format
.import.2bit <- function(file, selection){
    tbf <- TwoBitFile(file)
    if(!as.character(seqnames(selection)[1]) %in% seqnames(seqinfo(tbf))){
        return(DNAStringSet())
    }else{
        tmp <- import(tbf, which=selection)
        names(tmp) <- as.character(seqnames(selection)[1])
        return(tmp)
    }
}



## A mapping of (lower-cased) file extensions to import function calls. Most of those are already implemented in the rtracklayer package.
## If no mapping is found an error will be raised suggesting to provide a user-defined import function.
.registerImportFun <- function(file){
    fileExt <- .fileExtension(file)
    file <- path.expand(file)
    return(switch(fileExt,
                  "gff"=import.gff(file),
                  "gff1"=import.gff1(file),
                  "gff2"=import.gff2(file),
                  "gff3"=.import.gff3(file),
                  "gtf"=import.gff2(file),
                  "bed"=import.bed(file),
                  "bedgraph"=import.bedGraph(file),
                  "wig"=import.wig(file),
                  "bw"=.import.bw,
                  "bigwig"=.import.bw,
                  "bam"=.import.bam,
                  stop(sprintf("No predefined import function exists for files with extension '%s'. Please manually provide an import function.",
                               fileExt))))
}


## Get the file extension for a file, taking into account potential gzipping
.fileExtension <- function(file){
    if(!grepl("\\.", file))
        stop("Unable to identify extension for file '", file, "'")
     ext <- sub(".*\\.", "", sub("\\.gz$|\\.gzip$", "", basename(file)))
    if(ext=="")
        stop("Unable to identify extension for file '", file, "'")
    return(tolower(ext))
}


availableDefaultMapping <- function(file, trackType){
    .checkClass(file, "character", 1)
    .checkClass(trackType, "character", 1)
    ext <- tolower(if(grepl("\\.", file)) .fileExtension(file) else file)
    vm <- .defaultVarMap(ext, trackType, justMap=TRUE)
    vm[[".stream"]] <- NULL
    return(vm)
}


## Helper function to handle defaults function arguments
.covDefault <- function(x, cov, def){
    res <- if(missing(x)){
        if(is.null(cov)){
            def
        }else{
            cov
        }
    }else{x}
    return(res)
}


## A helper function to process alignment information from a GRanges object
.computeAlignments <- function(range){
    res <- list(range=range, stackRanges=GRanges(), stacks=numeric())
    if(length(range)){
        alg <- extractAlignmentRangesOnReference(range$cigar)
        rp <- elementNROWS(alg)
        range <- sort(GRanges(seqnames=rep(seqnames(range), rp), strand=rep("*", sum(rp)), ranges=shift(unlist(alg), rep(start(range), rp)-1),
                              id=rep(range$id, rp), entityId=rep(seq_along(rp), rp), cigar=rep(range$cigar, rp), md=rep(range$md, rp),
                              readStrand=rep(strand(range), rp), mapq=rep(range$mapq, rp), flag=rep(range$flag, rp), isize=rep(range$isize, rp),
                              groupid=rep(range$groupid, rp), status=factor(rep(range$status, rp), levels=c("mated", "ambiguous", "unmated")),
                              uid=seq_len(sum(rp))))
        if(length(range)){
            stTmp <- split(range, range$groupid)
            stackRanges <- unlist(range(stTmp))
            ss <- disjointBins(stackRanges)
            range$stack <- ss[match(range$groupid, names(stackRanges))]
            res <- list(range=range, stackRanges=stackRanges, stacks=range$stack)
        }
    }
    return(res)
}

## Check whether the current device supports alpha channel transparency
.supportsAlpha <- function(){
    d <- dev.cur()
    oldwarn <- getOption("warn")
    on.exit({options(warn=oldwarn)
             if(d==1)
                 dev.off()
         })
    options(warn=2)
    ok <- !is(try({grob <- grid.rect(width=0, height=0, gp=gpar(alpha=0.5))
                   ##grid.remove(grob$name)
               }, silent=TRUE), "try-error")
    return(ok)
}

## Return the alpha display parameter from a GdObject, respecting whether transparency is supported on the device
## or not. This is either drawn from the internal '.__hasAlphaSupport' display parameter, or, if not set, is
## determined dynamically.
.alpha <- function(GdObject, postfix=NULL){
    support <- .dpOrDefault(GdObject, ".__hasAlphaSupport", .supportsAlpha())
    wh <- if(is.null(postfix)) "alpha" else c(paste("alpha", postfix, sep="."), "alpha")
    alpha <- .dpOrDefault(GdObject, wh, 1)
    if(alpha != 1 && !support)
        alpha <- 1
    return(alpha)
}

## Draw horizontal arrows into a viewport indicating cropped read alignments
.moreInd <- function(n=3, direction="up", ...){
    nn <- n * 2 + 1
    x <- rep(seq(1/nn, 1 - (1/nn), len=nn-2)[seq(1, nn-2, by=2)], each=n) + c(-1/nn/2, 0, 1/nn/2)
    y <- rep(if(direction == "up") c(0, 1, 0) else c(1, 0, 1), n)
    grid.polyline(x, y, id=rep(1:n, each=3), gp=gpar(...))
}



## Compute mismatches for AlignmentsTracks based on the read sequences and a reference sequence
.findMismatches <- function(GdObject){
    rgo <- .dpOrDefault(GdObject, ".__plottingRange")
    mmPos <- mmSamp <- mmSeq <- mmStack <- NULL
    if(!is.null(rgo)){
        ref <- as.character(as(subseq(GdObject@referenceSequence, start=rgo["from"], end=rgo["to"]), "Rle"))
        cm <- consensusMatrix(GdObject@sequences, as.prob=FALSE, baseOnly=TRUE)[-5,]
        cmm <- colMaxs(cm)
        css <- colSums(cm)
        cmp <- rbind(t(t(cm)/css), 0)
        rownames(cmp)[5] <- "N"
        sel <- is.na(cmp["A",])
        cmp[,sel] <- 0
        cmp["N", sel] <- 1
        consStr <- strsplit(consensusString(cmp), "")[[1]]
        varRegs <- which(cmm != css | (consStr != "N" & consStr != ref))
        if(length(varRegs)){
            rvg <- ref[varRegs]
            sel <- rvg != "-" & rvg != "N"
            if(any(sel)){
                varRegs <- varRegs[sel]
                rvg <- rvg[sel]
                mmTab <- t(sapply(varRegs, function(x) as.character(subseq(GdObject@sequences, x, width=1))))
                isMm <- t(rvg != "-" & mmTab != "+" & mmTab != "-" & mmTab != rvg)
                mmRelPos <- col(isMm)[which(isMm)]
                mmPos <- varRegs[mmRelPos] + rgo["from"] - 1
                mmSampInd <- row(isMm)[which(isMm)]
                mmSamp <- rownames(isMm)[mmSampInd]
                mmSeq <- mmTab[ncol(isMm) * (mmSampInd - 1) + mmRelPos]
                mmStack <- stacks(GdObject)[match(mmSamp, ranges(GdObject)$entityId)]
            }
        }
    }
    return(data.frame(position=mmPos, stack=mmStack, read=mmSamp, base=as.character(mmSeq), stringsAsFactors=TRUE))
}





## Return the default mappings between the metadata columns of an imported GRanges object and those
## of the track's GRanges object.
.defaultVarMap <- function(inputType, trackType, stream, fromUser=FALSE, justMap=FALSE){
    vm <- list(gtf=list(GeneRegionTrack=list(feature="type",
                                             gene=c("gene_id", "gene_name"),
                                             exon=c("exon_name", "exon_id"),
                                             transcript=c("transcript_name", "transcript_id"),
                                             symbol=c("gene_name", "gene_id"))),
               gff=list(AnnotationTrack=list(feature="type",
                                             group="group"),
                        GeneRegionTrack=list(feature="type",
                                             transcript="group")),
               gff1=list(AnnotationTrack=list(feature="type",
                                              group=group),
                        GeneRegionTrack=list(feature="type",
                                             transcript="group")),
               gff2=list(AnnotationTrack=list(feature="type",
                                              group=c("group", "Parent"),
                                              id=c("ID", "Name", "Alias")),
                        GeneRegionTrack=list(feature="type",
                                             gene=c("gene_id", "gene_name"),
                                             exon=c("exon_name", "exon_id"),
                                             symbol=c("gene_name", "gene_id"))),
               gff3=list(AnnotationTrack=list(feature="type",
                                              id=c("ID", "Name", "Alias"),
                                              group="Parent"),
                         GeneRegionTrack=list(feature="type",
                                              gene=c("gene_id", "gene_name"),
                                              exon=c("exon_name", "exon_id", "ID"),
                                              transcript=c("transcript_name", "transcript_id", "Parent"),
                                              symbol=c("gene_name", "gene_id", "Name", "Alias"))),
               bedgraph=list(DataTrack=list(score="score")),
               wig=list(DataTrack=list(score="score")),
               bed=list(AnnotationTrack=list(feature="itemRgb", id="name")),
               bigwig=list(DataTrack=list(score="score",
                                          .stream=TRUE)),
               bw=list(DataTrack=list(score="score",
                                      .stream=TRUE)),
               bam=list(DataTrack=list(score="score",
                                      .stream=TRUE),
                        AnnotationTrack=list(id="id",
                                             group="group",
                                             .stream=TRUE),
                        AlignmentsTrack=list(id="id",
                                             cigar="cigar",
                                             mapq="mapq",
                                             flag="flag",
                                             isize="isize",
                                             groupid="groupid",
                                             status="status",
                                             md="md",
                                             seq="seq",
                                             .stream=TRUE)))
    if(justMap)
        return(vm[[inputType]][[trackType]])
    if(fromUser){
        vm[[inputType]] <- setNames(list(list(".stream"=stream)), trackType)
    }else{
        if(is.null(vm[[inputType]]) || is.null(vm[[inputType]][[trackType]])){
            warning(sprintf(paste("There are no default mappings from %s files to %s. Please provide a manual mapping",
                                  "in the track constructor if you haven't already done so."),
                            inputType, trackType))
            vm[[inputType]] <- setNames(list(list(".stream"=stream)), trackType)
        }
    }
    return(vm[[inputType]][[trackType]])
}


## Helper function to go through the metadata columns of a DataFrame and match their colnames to a mapping if they
## are available
.resolveColMapping <- function(data, args, defMap){
    colnames(mcols(data)) <- paste(colnames(mcols(data)), "orig", sep="__.__")
    for(i in names(defMap)){
        if(is.character(args[[i]]) && length(args[[i]])==1 && paste(args[[i]], "orig", sep="__.__") %in% colnames(mcols(data))){
            defMap[[i]] <- args[[i]]
            args[[i]] <- NULL
        }
        mt <- match(paste(defMap[[i]], "orig", sep="__.__"), colnames(mcols(data)))
        mt <- mt[!is.na(mt)][1]
        if(!is.na(mt))
            mcols(data)[[i]] <- mcols(data)[,mt]
    }
    mcols(data) <- mcols(data)[, !grepl("__.__", colnames(mcols(data))), drop=FALSE]
    return(list(data=data, args=args, defMap=defMap))
}


## For an AnnotationTrack or a GeneRegionTrack, compute the actual ranges for a complete range element group
## (i.e., a whole transcript or track group) and also add the necessary space for the text label if needed.
## We add this information to the internal display parameters '.__groupRanges', '.__groupLabels' and
## '.__groupLabelWidths'. This function has to be called before stacks are being computed because 'setStacks'
## will use the values in '.__groupRanges'.
## Note that the computed ranges are not quite right because we are only crudely guessing the size of the
## title panels at this stage.
.computeGroupRange <- function(GdObject, hasAxis=FALSE, hasTitle=.dpOrDefault(GdObject, "showTitle", TRUE), title.width=1){
    if(is(GdObject, "AnnotationTrack")){
        finalRanges <- IRanges()
        GdObjectOrig <- GdObject
        GdObject <- GdObject[seqnames(GdObject) == chromosome(GdObject)]
        pr <- .dpOrDefault(GdObject, ".__plottingRange", data.frame(from=min(start(GdObject)), to=max(end(GdObject))))
        if(is.null(title.width))
            title.width <- 1
        if(length(GdObject) > 0){
            gp <- group(GdObject)
            needsGrp <- any(duplicated(gp))
            finalRanges <- if(needsGrp){
                groups <- split(range(GdObject), gp)
                unlist(range(groups))
            }else{
                range(GdObject)
            }
            if(.dpOrDefault(GdObject, ".__hasAnno", FALSE)){
                ## The label justification
                just <- .dpOrDefault(GdObject, "just.group", "left")
                rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
                ## A crude guestimate of the space needed for a title
                twidth <- if(hasTitle){
                    fact <- title.width + (hasAxis * 2)
                    .getStringDims(GdObject, "g_T", unit="npc", subtype="title")$height * fact
                }else 0
                tfact <- ifelse(twidth > 1, 1, 1 / (1 - twidth))
                ## The labels and spacers are plotted in a temporary viewport to figure out their size
                labels <- if(needsGrp)
                    sapply(split(identifier(GdObject), gp), function(x)
                           paste(sort(unique(x)), collapse="/")) else setNames(identifier(GdObject), gp)
                xscale <- c(max(pr["from"], min(start(finalRanges))), min(pr["to"], max(end(finalRanges))))
                if(diff(xscale) == 0)
                    xscale[2] <- xscale[2] + 1
                pushViewport(dataViewport(xscale=xscale, extension=0, yscale=c(0,1), gp=.fontGp(GdObject, "group")))
                labelWidths <- setNames(as.numeric(convertWidth(stringWidth(labels),"native")) * tfact * 1.3, names(labels))
                spaceBefore <- as.numeric(convertWidth(unit(3, "points"),"native")) * tfact
                spaceAfter <- as.numeric(convertWidth(unit(7, "points"),"native")) * tfact
                popViewport(1)
                switch(as.character(just),
                       "left"={
                           if(!rev){
                               start(finalRanges) <- start(finalRanges) - (spaceBefore + labelWidths + spaceAfter)
                           }else{
                               end(finalRanges) <- end(finalRanges) + spaceAfter + labelWidths + spaceBefore
                           }
                           sb <- spaceBefore
                           sa <- spaceAfter
                       },
                       "right"={
                           if(!rev){
                               end(finalRanges) <-  end(finalRanges) + spaceAfter + labelWidths + spaceBefore
                           }else{
                               start(finalRanges) <- start(finalRanges) - (spaceBefore + labelWidths + spaceAfter)
                           }
                           sb <- spaceBefore
                           sa <- spaceAfter
                       },
                       "above"={
                           if(!rev){
                               featureWidths <- end(finalRanges) - start(finalRanges)
                               additionalLabelSpace <- ceiling((labelWidths - featureWidths) / 2)
                               additionalLabelSpace[additionalLabelSpace < 0] <- 0
                               end(finalRanges) <- end(finalRanges) + additionalLabelSpace
                               start(finalRanges) <- start(finalRanges) - additionalLabelSpace
                           }else{
                               featureWidths <- start(finalRanges) - end(finalRanges)
                               additionalLabelSpace <- ceiling((labelWidths - featureWidths) / 2)
                               additionalLabelSpace[additionalLabelSpace < 0] <- 0
                               end(finalRanges) <- end(finalRanges) - additionalLabelSpace
                               start(finalRanges) <- start(finalRanges) + additionalLabelSpace
                           }
                           sa <- sb <- 0
                       },
                       "below"={
                           if(!rev){
                               featureWidths <- end(finalRanges) - start(finalRanges)
                               additionalLabelSpace <- ceiling((labelWidths - featureWidths) / 2)
                               additionalLabelSpace[additionalLabelSpace < 0] <- 0
                               end(finalRanges) <- end(finalRanges) + additionalLabelSpace
                               start(finalRanges) <- start(finalRanges) - additionalLabelSpace
                           }else{
                               featureWidths <- start(finalRanges) - end(finalRanges)
                               additionalLabelSpace <- ceiling((labelWidths - featureWidths) / 2)
                               additionalLabelSpace[additionalLabelSpace < 0] <- 0
                               end(finalRanges) <- end(finalRanges) - additionalLabelSpace
                               start(finalRanges) <- start(finalRanges) + additionalLabelSpace
                           }
                           sa <- sb <- 0
                       },
                       stop(sprintf("Unknown label justification '%s'", just)))
                displayPars(GdObjectOrig) <- list(".__groupLabelWidths"=data.frame(before=sb, label=labelWidths, after=sa),
                                                  ".__groupLabels"=labels)
            }
        }
        displayPars(GdObjectOrig) <- list(".__groupRanges"=finalRanges)
    }
    return(GdObjectOrig)
}

## Calculate the vectorized string dimensions and return the results in a data.frame with columns 'width' and 'height'
## The unit of the return value can be controlled, and additional parameters like font size and expansion factors can
## be passed in as additional arguments (all in ... is passed on to 'gpar'). If needed, the font defaults for a
## subtype can be extracted by providing the subtype argument.
.getStringDims <- function(GdObject, string, unit="native", subtype=NULL, ...){
    gp <- .fontGp(GdObject, subtype, ...)
    pushViewport(viewport(gp=gp, xscale=current.viewport()$xscale, yscale=current.viewport()$yscale))
    res <- data.frame(width=as.numeric(convertWidth(stringWidth(string), unit)),
               height=as.numeric(convertHeight(stringHeight(string), unit)))
    popViewport(1)
    return(res)
}


## Check whether transcripts are to be collapsed for a GeneRegionTrack
.transcriptsAreCollapsed <- function(GdObject){
    res <- FALSE
    if(is(GdObject, "GeneRegionTrack")){
       ctrans <- .dpOrDefault(GdObject, "collapseTranscripts", FALSE)
       res <- (is.logical(ctrans) && ctrans == TRUE) || ctrans %in% c("gene", "shortest", "longest", "meta")
   }
    return(res)
}

## Create list for drawing sashimi-like plots
## using summarizeJunctions on GAlignments
## plotting is done via grid.xspline (requires x, y, id, score)
.sashimi.junctions <- function(range, score=1L, lwd.max=10, strand="*", filter=NULL, filterTolerance=0L) {
    ## summarizeJunctions
    range <- sort(range)
    range <- range[!duplicated(range$entityId)]
    ga <- GAlignments(seqnames=seqnames(range), pos=start(range), cigar=range$cigar,
                      strand=if(is.null(range$readStrand)) strand(range) else range$readStrand,
                      seqlengths=seqlengths(range))
    juns <- summarizeJunctions(ga)
    ## filter junctions
    if (!is.null(filter)) {
        ## if filterTolerance is > 0 than pass it as maxgap parameter in findOverlaps
        ## make sure it is positive value
        if (filterTolerance < 0) {
            filterTolerance <- abs(filterTolerance)
            warning(sprintf("\"sashimiFilterTolerance\" can't be negative, taking absolute value of it: %d",
                    filterTolerance))
        }
        ovs <- findOverlaps(juns, filter, type="start", maxgap=filterTolerance)
        ove <- findOverlaps(juns, filter, type="end", maxgap=filterTolerance)
        ## combine both Hits objects, select junctions present in both
        ovv <- rbind(as.matrix(ovs), as.matrix(ove))
        ovv <- ovv[duplicated(ovv),,drop=FALSE]
        ## create row/col index for selecting the correct strand
        ws <- strand(filter[ovv[,"subjectHits"]])
        levels(ws) <- c("plus_score", "minus_score", "score")
        ws <- cbind(row=ovv[,"queryHits"], col=as.character(ws))
        ## rol/col subseting will only works on matrix
        M <- as.matrix(values(juns))
        rownames(M) <- 1:nrow(M)
        ##
        filter$score <- 0L
        filter$score[sort(unique(ovv[,"subjectHits"]))] <- tapply(M[ws], ovv[,"subjectHits"], sum)
        juns <- filter
    } else {
        ## if no filter ranges were defined
        ## select strand (default both)
        juns$score <- if (strand=="+") juns$plus_score else if (strand=="-") juns$minus_score else juns$score
    }
    ## filter based on evidence (default no filtering, 1 read)
    juns <- juns[juns$score >= score]
    if (length(juns)) {
        ## count how many overlaps to determine the y
        ov <- findOverlaps(juns, reduce(juns, min.gapwidth=0L))
        ov <- split(queryHits(ov), subjectHits(ov))
        juns$y <- unlist(sapply(ov, order))
        ## scale the score to lwd.max
        juns$scaled <- (lwd.max-1)/pmax((max(juns$score)-min(c(1, juns$score))), 1)*(juns$score-max(juns$score))+lwd.max
        ## create list
        juns <- list(x=as.numeric(rbind(start(juns),
                         mid(ranges(juns)), end(juns))),
                     y=as.numeric(rbind(0, juns$y, 0)),
                     id=rep(seq_len(length(juns)), each=3),
                     score=juns$score,
                     scaled=juns$scaled)
    } else {
        juns <- list(x=numeric(),
                     y=numeric(),
                     id=integer(),
                     score=numeric(),
                     scaled=numeric())
    }
    return(juns)
}

