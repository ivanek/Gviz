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
##    - a single integer or a character coercable to one
##    - a character, starting with 'chr' (case insensitive)
## Arguments:
##    o x: a character string to be converted to a valid UCSC chromosome name
## Value: the UCSC character name
.chrName <- function(x)
{
    if(!getOption("ucscChromosomeNames"))
        return(as.character(x))
    xu <- unique(x)
    xum <- sapply(xu, function(y){
        xx <- suppressWarnings(as.integer(y))
        if(!is.na(xx))
            y <- xx
        if(is.numeric(y))
            y <- paste("chr", y, sep="")
        substring(y, 1,3) <- tolower(substring(y, 1,3))
        head <- sapply(y, substring, 1,3) == "chr"
        if(!all(head))
            stop(sprintf(paste("Invalid chromosome identifier%s '%s'\nPlease consider setting options(ucscChromosomeNames=FALSE)",
                               "to allow for arbitrary chromosome identifiers."),
                         ifelse(sum(!head)>1, "s", ""),
                         paste(y[!head], collapse=", ")))
        y})
    names(xum) <- xu
    return(as.vector(xum[as.character(x)]))
}


## Make a deep copy of the display parameter environment
.deepCopyPars <- function(GdObject)
{
    oldPars <- displayPars(GdObject)
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
        type <- match.arg(.dpOrDefault(x, "type", "p"), c("p", "l", "b", "a", "s", "g", "r", "S", "smooth",
                                                          "histogram", "mountain", "h", "boxplot", "gradient", "heatmap", "polygon"),
                          several.ok=TRUE)
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
      size <- max(size, size*max(stacks(x)))
    return(size)
}



## Return a particular displayPars value or a default
## Arguments:
##    o GdObject: an object inheriting from class GdObject
##    o par: the name of the displayPar
##    o default: a default value for the parameter if it can't be found in GdObject
## Value: the value of the displayPar
.dpOrDefault <- function(GdObject, par, default=NULL, fromPrototype=FALSE)
{
    val <- getPar(GdObject, par)
    if(is.null(val)) {
       if (fromPrototype) {
          val <- Gviz:::.parMappings[[GdObject@name]][[par]]
       } else {
          val <- default
       }
    }
    return(val)
}



## Check a list of GdObjects whether an axis needs to be drawn for each of them.
## Arguments:
##    o object: a list of GdObjects
## Value: a logical vector of the same length as 'objects'
.needsAxis <- function(objects) {
    if(!is.list(objects))
        objects <- list(objects)
    atrack <- sapply(objects, function(x){
        type <- match.arg(.dpOrDefault(x, "type", "p"), c("p", "l", "b", "a", "s", "g", "r", "S", "smooth",
                                                         "histogram", "mountain", "h", "boxplot", "gradient", "heatmap", "polygon"),
                          several.ok=TRUE)
        is(x, "NumericTrack") || (is(x, "AlignedReadTrack") && .dpOrDefault(x, "detail", "coverage")=="coverage")})
    return(atrack & sapply(objects, .dpOrDefault, "showAxis", TRUE))
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
.setupTextSize <- function(trackList, sizes, title.width, panelOnly=FALSE)
{
    curVp <- vpLocation()
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
        width <- (as.numeric(convertWidth(stringHeight(paste("g_T", nwrap, "g_T", sep="")), "inches"))*cex+leaveSpace)/wfac
        showtitle <- sapply(trackList, .dpOrDefault, "showTitle", TRUE)
        width[!showtitle & !needAxis] <- 0
        width[!showtitle] <- 0
        twfac <- if(missing(title.width) || is.null(title.width)) 1 else title.width
        title.width <- max(width, na.rm=TRUE)
        if(any(needAxis)){
            cex.axis <- structure(sapply(trackList, .dpOrDefault, "cex.axis", 0.6), names=nn)
            axTicks <- unlist(lapply(trackList, function(GdObject){
                if(!is(GdObject, "NumericTrack") && !is(GdObject, "AlignedReadTrack"))
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
                type <- match.arg(.dpOrDefault(GdObject, "type", "p"), c("p", "l", "b", "a", "s", "g", "r", "S", "smooth", "polygon",
                                                                         "histogram", "mountain", "h", "boxplot", "gradient", "heatmap"),
                                  several.ok=TRUE)
                if(any(c("heatmap", "gradient") %in% type)){
                    nlevs <- max(1, nlevels(factor(getPar(GdObject, "groups"))))-1
                    atSpace <- atSpace + 0.3 * atSpace + as.numeric(convertWidth(unit(3, "points"), "inches"))*nlevs
                }
                if(type=="heatmap" && .dpOrDefault(GdObject, "showSampleNames", FALSE)){
                    sn <- rownames(values(GdObject))
                    wd <- max(as.numeric(convertWidth(stringWidth(sn) + unit(10, "points"), "inches")))
                    atSpace <- atSpace + (wd * .dpOrDefault(GdObject, "cex.sampleNames", 0.5))
                }
                atSpace
            }))
            hAxSpaceNeeded <- (max(axTicks))/wfac
            title.width <- title.width +  hAxSpaceNeeded
        }
    } else {
        title.width <- nwrap <- cex <- NA
    }
    spacing <- 0.02
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



## Take coordinates for the bounding boxes of annotation regions and plot filled arrows inside.
## Arguments:
##    o coords as numeric matrix with 4 columns: x1, y1, x2, y2
##    o W: the proportion of the total box width used for the arrow head
##    o H: the proportion of the total box height used for the arrow head
##    o col: the boundary color
##    o fill: the fill color
##    o lwd: the boundary line width
##    o lty: the boundary line type
##    o alpha: the transparency
##    o strand: the strand information, a character of either "+" or "-"
##    o min.width: the minumum width of the arrow head. Below this size a simple box is drawn
## Note that the last arguments 4-9 all have to be of the same length as number of rows in coords.
## Value: the function is called for its side-effects of drawing on the graphics device
.filledArrow <- function(coords, W=1/4, H=1/3, col, fill, lwd, lty, alpha, strand=0, min.width=10) {
    A <- coords[,1:2,drop=FALSE]
    B <- coords[,3:4,drop=FALSE]
    ## First everything that is still a box
    osel <- abs(B[,1]-A[,1]) < min.width | !strand %in% c("+", "-")
    xx <-  c(A[osel,1], B[osel,1], B[osel,1], A[osel,1])
    ##offset <-  ifelse(strand[osel] %in% c("+", "-"), (abs(B[osel,2]-A[osel,2])*H/2), 0)
    offset <- (abs(B[osel,2]-A[osel,2])*H/2)
    yy <- c(rep(A[osel,2]+offset, 2), rep(B[osel,2]-offset, 2))
    id <- rep(seq_len(sum(osel)), 4)
    pars <- data.frame(fill=fill, col=col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE)[osel,]
    ## Now the arrows facing right
    sel <- !osel & strand=="+"
    id <- c(id, rep(seq(from=if(!length(id)) 1 else max(id)+1, by=1, len=sum(sel)), 7))
    xx <- c(xx, A[sel,1], rep(A[sel,1]+(abs(B[sel,1]-A[sel,1])*W),2), B[sel,1], rep(A[sel,1]+(abs(B[sel,1]-A[sel,1])*W),2), A[sel,1])
    yy <- c(yy, rep(A[sel,2]+(abs(B[sel,2]-A[sel,2])*H/2),2), A[sel,2], A[sel,2]+(abs(B[sel,2]-A[sel,2])/2), B[sel,2],
            rep(B[sel,2]-(abs(B[sel,2]-A[sel,2])*H/2),2))
    pars <- rbind(pars, data.frame(fill=fill, col=col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE)[sel,])
    ## And finally those facing left
    sel <- !osel & strand=="-"
    id <- c(id, rep(seq(from=if(!length(id)) 1 else max(id)+1, by=1, len=sum(sel)), 7))
    xx <- c(xx, B[sel,1], rep(B[sel,1]-((B[sel,1]-A[sel,1])*W),2), A[sel,1], rep(B[sel,1]-((B[sel,1]-A[sel,1])*W),2), B[sel,1])
    yy <- c(yy, rep(A[sel,2]+(abs(B[sel,2]-A[sel,2])*H/2),2), A[sel,2], A[sel,2]+(abs(B[sel,2]-A[sel,2])/2), B[sel,2],
            rep(B[sel,2]-(abs(B[sel,2]-A[sel,2])*H/2),2))
    pars <- rbind(pars, data.frame(fill=fill, col=col, lwd=lwd, lty=lty, alpha=alpha, stringsAsFactors=FALSE)[sel,])
    grid.polygon(x=xx, y=yy, gp=gpar(fill=pars$fill, col=pars$col, alpha=pars$alpha, lwd=pars$lwd, lty=pars$lty),
                 default.units="native", id=id)
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
                      diff=.pxResolution(coord="y"), min.height=3)
  {
      if(!barOnly)
      {
          onePx <- diff
          if(missing(H))
          {
              onePy <- .pxResolution(coord="y")
              H <- onePy*min.height/2
          }
          fx1 <- fx2 <- scol <- fy1 <- fy2 <- NULL
          exons <- IRanges(start=coords[,1], end=coords[,3])
          levels <- split(exons, coords[,2])
          for(i in seq_along(xx1))
          {
              x1 <- xx1[i]
              x2 <- xx2[i]
              len <- diff(c(x1,x2))/onePx
              if(len>D+W*2)
              {
                  ax1 <- seq(from=x1+(onePx*W), to=x1+(len*onePx)-(onePx*W), by=onePx*D)
                  ax2 <- ax1+(onePx*W)
                  feathers <- IRanges(start=ax1-onePx, end=ax2+onePx)
                  cur.level <- which(y[i]==unique(y))
                  sel <- queryHits(findOverlaps(feathers, levels[[cur.level]]))
                  if(length(sel))
                  {
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
      grid.segments(xx1, y, xx2, y, default.units="native", gp=gpar(col=col, lwd=lwd, lty=lty, alpha=alpha, lineend="square"))
  }



## Extract track color for different subtypes within the track and use the default
## color value if no other is found, lightblue if no colors are set at all
## Arguments:
##    o GdObject: object inheriting from class GdObject
## Value: a color character
.getBiotypeColor <- function(GdObject) {
    defCol <- .dpOrDefault(GdObject, "fill", Gviz:::.DEFAULT_FILL_COL)
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
                 fontfamily = fontfamily, fontface = lattice:::chooseFace(fontface,
                                            font), fontsize = fontsize)
  }
  panel.points(x = rep(levels.fos, sapply(blist.out, length)),
               y = unlist(blist.out), pch=pch, col=col,
               alpha=alpha, cex=cex,
               fontfamily=fontfamily, fontface = lattice:::chooseFace(fontface,
                                        font), fontsize = fontsize)
}



## Check which parameters have already been set for a GdObject, and
## update all missing ones from the prototype of the current parent
## class.
## Arguments:
##    o x: an object inheriting from class GdObject
##    o class: the parent class from which to draw the missing parameters
## Value: The updated GdObject
.updatePars <- function(x, class)
  {
    current <- getPar(x)
    defaults <- getPar(getClass(class)@prototype@dp)
    missing <- setdiff(names(defaults), names(current))
    displayPars(x) <- defaults[missing]
    return(x)
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



## Helper function to translate from a UCSC genome name to a Biomart data set. This also caches the mart
## object in order to speed up subsequent calls
## Arguments:
##    o genome: character giving the UCSC genome
## Value: A BiomaRt connection object
.genome2Dataset <- function(genome)
{
    ds <- c("mm9"="mmusculus_gene_ensembl",
            "hg19"="hsapiens_gene_ensembl",
            "felCat4"="fcatus_gene_ensembl",
            "galGal3"="ggallus_gene_ensembl",
            "panTro2"="ptroglodytes_gene_ensembl",
            "bosTau4"="btaurus_gene_ensembl",
            "canFam3"="cfamiliaris_gene_ensembl",
            "loxAfr3"="lafricana_gene_ensembl",
            "fr2"="trubripes_gene_ensembl",
            "cavPor3"="cporcellus_gene_ensembl",
            "equCab2"="ecaballus_gene_ensembl",
            "anoCar1"="acarolinensis_gene_ensembl",
            "calJac3"="cjacchus_gene_ensembl",
            "oryLat2"="olatipes_gene_ensembl",
            "monDom5"="mdomestica_gene_ensembl",
            "susScr2"="sscrofa_gene_ensembl",
            "ornAna1"="oanatinus_gene_ensembl",
            "oryCun2"="ocuniculus_gene_ensembl",
            "rn4"="rnorvegicus_gene_ensembl",
            "gasAcu1"="gaculeatus_gene_ensembl",
            "tetNig2"="tnigroviridis_gene_ensembl",
            "xenTro2"="xtropicalis_gene_ensembl",
            "danRer7"="drerio_gene_ensembl",
            "ci2"="cintestinalis_gene_ensembl",
            "dm3"="dmelanogaster_gene_ensembl",
            "ce6"="celegans_gene_ensembl",
            "sacCer2"="scerevisiae_gene_ensembl")
    if(!tolower(genome) %in% tolower(names(ds)))
        stop("Unable to automatically determine Biomart data set for UCSC genome '", genome, "'")
    thisDs <- ds[match(tolower(genome), tolower(names(ds)))]
    cenv <- environment()
    bm <- .doCache(thisDs, expression(useMart("ensembl", dataset=thisDs)), .ensemblCache, cenv)
    return(bm)
}


.annotationSpace <- function(GdObject, from, to)
{
    if(!length(GdObject))
        return(NULL)
    sel <- chromosome(GdObject) == seqnames(GdObject) & start(GdObject)>=from & start(GdObject)<=to
    if (sum(sel) > 0) {
      ids <- identifier(GdObject, FALSE)[sel]
      hasAnno <- .dpOrDefault(GdObject, "showId", FALSE) & ids!=""
      txt <- paste(ids, " ")
      txt[!hasAnno] <- ""
      space <- (as.numeric(convertWidth(stringWidth(txt),"native"))*1.3)
      newFrom <- min(start(GdObject)[sel]-space*1.6)
    } else {
      txt <- ""
      space <- (as.numeric(convertWidth(stringWidth(txt),"native"))*1.3)
      newFrom <- from
    }
    cex <- .dpOrDefault(GdObject, "cex", 1) * .dpOrDefault(GdObject, "cex.symbol", 0.7)
    fontfamily <- .dpOrDefault(GdObject, "fontfamily", 1)
    fontsize <- .dpOrDefault(GdObject, "fontsize", 12)
    fontface <- .dpOrDefault(GdObject, "fontface.symbol", 2)
    pushViewport(dataViewport(xData=c(from, to), extension=0, yscale=c(0, 1), clip=TRUE,
                              gp=gpar(cex=cex, fontfamily=fontfamily, fonface=fontface, fontsize=fontsize)))
    space <- (as.numeric(convertWidth(stringWidth(txt),"native"))*1.3)
    popViewport(1)
    ## newFrom <- min(start(GdObject)[sel]-space*1.6)
    return(newFrom)
}


## Return the plotting range for a GdObject, either from the contained ranges or from overrides.
## This function is vectorized and should also work for lists of GdObjects.
.defaultRange <- function(GdObject, from=NULL, to=NULL, extend.left=0, extend.right=0, factor=0.01, annotation=FALSE)
{
    if(!is.list(GdObject))
        GdObject <- list(GdObject)
    if(!length(GdObject) || !all(sapply(GdObject, is, "GdObject")))
        stop("All items in the list must inherit from class 'GdObject'")
    tfrom <- lapply(GdObject, function(x){tmp <- start(x); if(is(x, "RangeTrack")) tmp <- tmp[seqnames(x)==chromosome(x)]; tmp})
    tfrom <- if(is.null(unlist(tfrom))) Inf else min(sapply(tfrom[listLen(tfrom)>0], min))
    tto <- lapply(GdObject, function(x){tmp <- end(x); if(is(x, "RangeTrack")) tmp <- tmp[seqnames(x)==chromosome(x)]; tmp})
    tto <- if(is.null(unlist(tto))) Inf else max(sapply(tto[listLen(tto)>0], max))
    if((is.null(from) || is.null(to)) && ((is.infinite(tfrom) || is.infinite(tto)) || is(GdObject, "GenomeAxisTrack")))
        stop("Unable to determine plotting ranges from the supplied track(s)")
    range <- extendrange(r=c(tfrom, tto), f=factor)
    range[1] <- max(1, range[1])
    wasNull <- FALSE
    if(is.null(from))
    {
        wasNull <- TRUE
        from <- range[1]
    }
    if(is.null(to))
        to <- range[2]
    from <- from-extend.left
    to <- to+extend.right
    if(from>to)
        stop("'from' range can not be larger than 'to'")
    ## We may need some extra space for annotations
    if(annotation)
    {
        annStarts <- unlist(lapply(GdObject[sapply(GdObject, is, "AnnotationTrack")], .annotationSpace, from, to))
        ## FIXME: Do we want to add annotation space if from was defined by the user?
        if(!is.null(annStarts) && wasNull)
            from <- min(c(from, annStarts))
    }
    return(c(from=as.vector(from), to=as.vector(to)))
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
        col <- .dpOrDefault(GdObject, "col", trellis.par.get("superpose.line")$col)
        col <- rep(col, length(groups))
        col.line <- rep(.dpOrDefault(GdObject, "col.line", col), length(groups))
        col.symbol <- rep(.dpOrDefault(GdObject, "col.symbol", col), length(groups))
        lwd <- rep(lwd, length(groups))
        lty <- rep(lty, length(groups))
        pch <- rep(pch, length(groups))
        cex <- rep(cex, length(groups))
    }
    col.baseline <- .dpOrDefault(GdObject, "col.baseline", col)
    col.grid <- .dpOrDefault(GdObject, "col.grid", "#e6e6e6")[1]
    fill <- .dpOrDefault(GdObject, "fill", Gviz:::.DEFAULT_FILL_COL)[1]
    fill.histogram <- .dpOrDefault(GdObject, "fill.histogram", fill)[1]
    col.histogram  <- .dpOrDefault(GdObject, "col.histogram", .dpOrDefault(GdObject, "col", Gviz:::.DEFAULT_SHADED_COL))[1]
    lty.grid <- .dpOrDefault(GdObject, "lty.grid", 1)
    lwd.grid <- .dpOrDefault(GdObject, "lwd.grid", 1)
    return(list(col=col, col.line=col.line, col.symbol=col.symbol, col.baseline=col.baseline,
                col.grid=col.grid, col.histogram=col.histogram, fill=fill, fill.histogram=fill.histogram,
                lwd=lwd, lty=lty, pch=pch, cex=cex, lwd.grid=lwd.grid, lty.grid=lty.grid))
}


.legendInfo <- function()
{
    legInfo <- matrix(FALSE, ncol=7, nrow=15, dimnames=list(c("p", "b", "l", "a", "s", "S", "r", "h", "smooth",
                                                              "histogram", "boxplot", "heatmap", "gradient", "mountain", "g"),
                                                            c("lty", "lwd", "pch", "col", "cex", "col.lines", "col.symbol")))
    legInfo[2:9, c("lty", "lwd", "col.lines")] <- TRUE
    legInfo[1:2, c("pch", "cex", "col.symbol")] <- TRUE
    legInfo[1:12, "col"] <- TRUE
    return(legInfo)
}


## Plot a list of GdObjects as individual tracks similar to the display on the UCSC genome browser
## Arguments:
##    o trackList: a list of GdObjects
##    o from, to: the plotting range, will be figured out automatically from the tracks if missing
##    o sized: a vector of relative vertical sizes, or NULL to auto-detect
##    o panel.only: don't draw track titles, useful to embed in a lattice-like function
##    o extend.right, extend.left: extend the coordinates in 'from' and 'too'
##    o title.width: the expansion factor for the width of the title track
## Value: the function is called for its side-effect of drawing on the graphics device
plotTracks <- function(trackList, from=NULL, to=NULL, ..., sizes=NULL, panel.only=FALSE, extend.right=0,
                       extend.left=0, title.width=NULL, add=FALSE, main, cex.main=2, fontface.main=2,
                       col.main="black", margin=6, chromosome=NULL)
{
    if(!is.list(trackList))
        trackList <- list(trackList)
    ## We first run very general housekeeping tasks on the tracks for which we don't really need to know anything about device
    ## size, resolution or plotting ranges. Chromosomes should all be the same for all tracks, if not we will force them to
    ## be set to the first one that can be detected
    chrms <- unlist(lapply(trackList, Gviz::chromosome))
    if(is.null(chromosome)){
        chrms <- if(!is.null(chrms)) chrms[gsub("^chr", "", chrms)!="NA"] else chrms
        chromosome <- head(chrms, 1)
        if(length(chromosome)==0)
            chromosome <- "chrNA"
        if(!is.null(chrms) && length(unique(chrms))!=1)
            warning("The track chromosomes in 'trackList' differ. Setting all tracks to chromosome '", chromosome, "'", sep="")
    }
    ## If plotting ranges are supplied we can speed up a lot of the downstream operations by subsetting first
    if(!is.null(from) || !(is.null(to))){
        trackList <- lapply(trackList, subset, from=from, to=to, chromosome=chromosome, sort=FALSE, stacks=FALSE, use.defaults=FALSE)
    }
    trackList <- lapply(trackList, consolidateTrack, chromosome=chromosome, ...)
    ## Now we figure out the plotting ranges. If no ranges are given as function arguments we take the absolute min/max of all tracks.
    if(!panel.only && !add)
        grid.newpage()
    ranges <- .defaultRange(trackList, from=from, to=to, extend.left=extend.left, extend.right=extend.right, annotation=TRUE)
    ## We need to reverse the list to get a top to bottom plotting order
    trackList <- rev(trackList)
    map <- vector(mode="list", length=length(trackList))
    titleCoords <- NULL
    names(map) <- rev(sapply(trackList, names))
    ## Now we can subset all the objects in the list to the current boundaries and compute the initial stacking
    trackList <- lapply(trackList, subset, from=ranges["from"], to=ranges["to"], chromosome=chromosome)
    trackList <- lapply(trackList, setStacks, from=ranges["from"], to=ranges["to"])
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
        spaceSetup <- .setupTextSize(trackList, sizes, title.width)
    } else {
        vpBound <- viewport()
        pushViewport(vpBound)
        spaceSetup <- .setupTextSize(trackList, sizes)
    }
    ## First iteration to set up all the dimensions by calling the drawGD methods in prepare mode, i.e.,
    ## argument prepare=TRUE. Nothing is drawn at this point, and this only exists to circumvent the
    ## chicken and egg problem of not knowing how much space we need until we draw, but also not knowing
    ## where to draw until we know the space needed.
    for(i in rev(seq_along(trackList)))
    {
        fontSettings <- .fontGp(trackList[[i]], cex=NULL)
        vpTrack <-  viewport(x=0, y=sum(spaceSetup$spaceNeeded[1:i]), just=c(0,1), width=1, height=spaceSetup$spaceNeeded[i],
                             gp=fontSettings)
        pushViewport(vpTrack)
        vpContent <- if(!panel.only) viewport(x=spaceSetup$title.width+spaceSetup$spacing,
                                              width=1-spaceSetup$title.width-spaceSetup$spacing, just=0) else viewport(width=1)
        pushViewport(vpContent)
        trackList[[i]] <- drawGD(trackList[[i]], minBase=ranges["from"], maxBase=ranges["to"], prepare=TRUE, subset=FALSE)
        popViewport(2)
    }
    ## Now lets recalculate the space and draw for real
    spaceSetup <- .setupTextSize(trackList, sizes, title.width)
    for(i in rev(seq_along(trackList)))
    {
        fontSettings <- .fontGp(trackList[[i]], cex=NULL)
        vpTrack <-  viewport(x=0, y=sum(spaceSetup$spaceNeeded[1:i]), just=c(0,1), width=1, height=spaceSetup$spaceNeeded[i],
                             gp=fontSettings)
        pushViewport(vpTrack)
		fill <- .dpOrDefault(trackList[[i]], "background.title", Gviz:::.DEFAULT_SHADED_COL)
        if(!panel.only) {
            vpTitle <- viewport(x=0, width=spaceSetup$title.width, just=0)
            pushViewport(vpTitle)
            lwd.border.title <- .dpOrDefault(trackList[[i]], "lwd.border.title", 1)
            col.border.title <- .dpOrDefault(trackList[[i]], "col.border.title", "transparent")
            grid.rect(gp=gpar(fill=fill, col=col.border.title, lwd=lwd.border.title))
            needAxis <- .needsAxis(trackList[[i]])
            drawAxis(trackList[[i]], ranges["from"], ranges["to"], subset=FALSE)
            tit <- spaceSetup$nwrap[i]
            titleCoords <- rbind(titleCoords, cbind(.getImageMap(cbind(0,0,1,1)),
                                                    title=names(trackList[[i]])))
            if(.dpOrDefault(trackList[[i]], "showTitle", TRUE) && !is.null(tit) && tit!="")
            {
                col <- .dpOrDefault(trackList[[i]], "col.title", "white")
                fontface <- .dpOrDefault(trackList[[i]], "fontface.title", 2)
                fontsize <- .dpOrDefault(trackList[[i]], "fontsize.title", 12)
                lineheight <- .dpOrDefault(trackList[[i]], "lineheight.title", 1)
                fontfamily <- .dpOrDefault(trackList[[i]], "fontfamily.title", "sans")
                lcex <- spaceSetup$cex[i]
                x <- if(needAxis) 0.075 else 0.4
                just <- if(needAxis) c("center", "top") else "center"
                suppressWarnings(grid.text(tit, unit(x, "npc"), rot=90, gp=gpar(col=col, fontface=fontface, cex=lcex,
                                                               fontsize=fontsize, lineheight=lineheight,
                                                               fontfamily=fontfamily), just=just))
            }
            popViewport(1)
        }
        ## Draw the panel background, grid lines if necessary and the panel content
        vpBackground <- if(!panel.only) viewport(x=spaceSetup$title.width,
                                                 width=1-spaceSetup$title.width, just=0) else viewport(width=1)
        pushViewport(vpBackground)
        grid.rect(gp=gpar(col="transparent", fill=.dpOrDefault(trackList[[i]], "background.panel", "transparent")))
        drawGrid(trackList[[i]], ranges["from"], ranges["to"])
        popViewport(1)
        vpContent <- if(!panel.only) viewport(x=spaceSetup$title.width+spaceSetup$spacing,
                                              width=1-spaceSetup$title.width-spaceSetup$spacing, just=0) else viewport(width=1)
        pushViewport(vpContent)
        tmp <- drawGD(trackList[[i]], minBase=ranges["from"], maxBase=ranges["to"], subset=FALSE)
        if(!is.null(tmp))
            map[[(length(map)+1)-i]] <- tmp
        popViewport(1)
        if(.dpOrDefault(trackList[[i]], "frame", FALSE))
            grid.rect(gp=gpar(col=.dpOrDefault(trackList[[i]], "col.frame", Gviz:::.DEFAULT_SHADED_COL)))
        popViewport(1)
    }

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


## Return the font settings for a GdObject
.fontGp <- function(GdObject, ...)
{
    gp <- list(fontsize=as.vector(.dpOrDefault(GdObject, "fontsize", 12))[1],
               fontface=as.vector(.dpOrDefault(GdObject, "fontface", 1))[1],
               lineheight=as.vector(.dpOrDefault(GdObject, "lineheight", 1))[1],
               fontfamily=as.character(as.vector(.dpOrDefault(GdObject, "fontfamily", 1)))[1],
               col=as.vector(.dpOrDefault(GdObject, "fontcolor", "white"))[1],
               alpha=as.vector(.dpOrDefault(GdObject, "alpha", 1))[1],
               cex=as.vector(.dpOrDefault(GdObject, "cex", 1))[1])
    gp[names(list(...))] <- list(...)
    gp <- gp[!sapply(gp, is.null)]
    class(gp) <- "gpar"
    return(gp)
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

.DEFAULT_FILL_COL <- "lightgray"
.DEFAULT_OVERPLOT_COL <- "red"
.DEFAULT_LINE_COL <- "black"
.DEFAULT_SHADED_COL <- "#808080"
.DEFAULT_SYMBOL_COL <- "#0080FF"


## We store some preset in the options on package load
.onLoad = function(...){options("ucscChromosomeNames"=TRUE)}


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
    dat <- import.gff3(file, asRangedData=FALSE)
    res <- try({
        genes <- tolower(dat$type) == "gene"
        ginfo <- mcols(dat[genes, ])
        dat <- dat[!genes]
        transcripts <- tolower(dat$type) == "mrna"
        tinfo <- mcols(dat[transcripts, ])
        dat <- dat[!transcripts]
        mt <- match(as.character(dat$Parent), as.character(tinfo$ID))
        if(!all(is.na(mt))){
            if(!"transcript_id" %in% colnames(mcols(dat)))
                mcols(dat)[["transcript_id"]] <- tinfo[mt, "ID"]
            if(!"transcript_name" %in% colnames(mcols(dat)))
                mcols(dat)[["transcript_name"]] <- tinfo[mt, "Name"]
            mt2 <- match(as.character(tinfo[mt, "Parent"]), as.character(ginfo$ID))
            if(!all(is.na(mt2))){
                if(!"gene_id" %in% colnames(mcols(dat)))
                    mcols(dat)[["gene_id"]] <- ginfo[mt2, "ID"]
                if(!"gene_name" %in% colnames(mcols(dat)))
                    mcols(dat)[["gene_name"]] <- ginfo[mt2, "Name"]
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
        rr <- import.bw(con=bwf, asRangedData=FALSE)
    }else{
        si <- seqinfo(bwf)
        rr <- if(!as.character(seqnames(selection)[1]) %in% seqnames(seqinfo(bwf))){
            GRanges(seqnames(selection)[1], ranges=IRanges(1,2), score=1)[0] }else{
                import.bw(con=bwf, selection=selection, asRangedData=FALSE)}
    }
    return(rr)
}

## An import function for bam files that distinguishes between DataTracks and AnnotationTracks
.import.bam <- function(file, selection){
    if(!file.exists(paste(file, "bai", sep=".")))
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
                  "gff"=import.gff(file, asRangedData=FALSE),
                  "gff1"=import.gff1(file, asRangedData=FALSE),
                  "gff2"=import.gff2(file, asRangedData=FALSE),
                  "gff3"=.import.gff3(file),
                  "gtf"=import.gff2(file, asRangedData=FALSE),
                  "bed"=import.bed(file, asRangedData=FALSE),
                  "bedgraph"=import.bedGraph(file, asRangedData=FALSE),
                  "wig"=import.wig(file, asRangedData=FALSE),
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

## Return the default mappings between elementMetadata slots of an imported GRanges object and the elementMetadata
## slots of the track's GRanges object.
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


## Helper function to go through the elementMetadata columns of a DataFrame and match their colnames to a mapping if they
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
