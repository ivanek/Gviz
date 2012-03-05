##----------------------------------------------------------------------------------------------------------------------------
## Some rather general accessors to extract information from all kinds of GdObjects
##----------------------------------------------------------------------------------------------------------------------------
## Extract the full GRanges object from the range slot of an object inheriting from RangeTrack
setMethod("ranges", "RangeTrack", function(x) x@range)
setMethod("ranges", "GenomeAxisTrack", function(x) x@range)

## Extract the IRanges part of the GRanges object from the range slot of an object inheriting from RangeTrack
setMethod("range", "RangeTrack", function(x) ranges(x@range))
setMethod("range", "GenomeAxisTrack", function(x) ranges(x@range))

## seqnames from the range track
setMethod("seqnames", "RangeTrack", function(x) as.character(seqnames(ranges(x))))

## Min and max ranges
setMethod("min", "RangeTrack", function(x) min(start(x)))
setMethod("max", "RangeTrack", function(x) max(end(x)))

## Extract start and end coordinates
setMethod("start", "RangeTrack", function(x) if(length(x)) as.integer(start(range(x))) else NULL)
setReplaceMethod("start", "RangeTrack", function(x, value) {
    start(x@range) <- value
    return(x)})
setReplaceMethod("start", "GenomeAxisTrack", function(x, value) {
    start(x@range) <- value
    return(x)})
setReplaceMethod("start", "IdeogramTrack", function(x, value) return(x))
setMethod("end", "RangeTrack", function(x) if(length(x)) as.integer(end(range(x))) else NULL)
setReplaceMethod("end", "RangeTrack", function(x, value) {
    end(x@range) <- value
    return(x)})
setReplaceMethod("end", "GenomeAxisTrack", function(x, value) {
    end(x@range) <- value
    return(x)})
setReplaceMethod("end", "IdeogramTrack", function(x, value) return(x))
setMethod("width", "RangeTrack", function(x) if(length(x)) as.integer(width(range(x))) else NULL)
setReplaceMethod("width", "RangeTrack", function(x, value) {
    width(x@range) <- value
    return(x)})
setReplaceMethod("width", "IdeogramTrack", function(x, value) return(x))
setMethod("start", "GenomeAxisTrack", function(x) if(length(x)) start(range(x)) else NULL)
setMethod("end", "GenomeAxisTrack", function(x) if(length(x)) end(range(x)) else NULL)
setMethod("width", "GenomeAxisTrack", function(x) if(length(x)) as.integer(width(range(x))) else NULL)
setMethod("start", "IdeogramTrack", function(x) NULL)
setMethod("end", "IdeogramTrack", function(x) NULL)
setMethod("width", "IdeogramTrack", function(x) NULL)

## Return the number of individual annotation items (independent of any grouping) in a RangeTrack
setMethod("length", "RangeTrack", function(x) sum(seqnames(x) == chromosome(x)))
setMethod("length", "GenomeAxisTrack", function(x) length(ranges(x)))
setMethod("length", "IdeogramTrack", function(x) length(ranges(x)))

## Extract the elementMetadata slot from the GRanges object of an object inheriting from RangeTrack as a data.frame.
## For a DataTrack object these values are stored as a numeric matrix in the data slot, and we return this instead.
setMethod("values", "RangeTrack", function(x) as.data.frame(values(ranges(x)), stringsAsFactors=FALSE))
setMethod("values", "GenomeAxisTrack", function(x) as.data.frame(values(ranges(x)), stringsAsFactors=FALSE))
setMethod("values", "DataTrack", function(x) x@data[,seqnames(x) == chromosome(x), drop=FALSE])
setReplaceMethod("values", "DataTrack", function(x, value){
    if(!is.matrix(value))
    {
        if(!is.numeric(value) || length(value) != length(x))
            stop("Invalid length of replacement vector.")
        if(!is.matrix(value) || !is.numeric(value) || ncol(value)!=length(x))
            stop("Dimensions of replacement value do not match.")
    }
    x@data <- value
    return(x)
})

                 
               

## Set or extract the chromosome from a RangeTrack object
setMethod("chromosome", "RangeTrack", function(GdObject) GdObject@chromosome)
setReplaceMethod("chromosome", "RangeTrack", function(GdObject, value){
    GdObject@chromosome <- .chrName(value[1])
    return(GdObject)
})
setReplaceMethod("chromosome", "IdeogramTrack", function(GdObject, value){
    message("Updating chromosome band information")
    tmp <- IdeogramTrack(genome=genome(GdObject), chromosome=.chrName(value[1]), name=names(GdObject))
    displayPars(tmp) <- displayPars(GdObject)
    return(tmp)
})

## Set or extract the chromosome from a RangeTrack object
setMethod("genome", "RangeTrack", function(x) x@genome)
setReplaceMethod("genome", "RangeTrack", function(x, value){
    x@genome <- value[1]
    return(x)
})
setReplaceMethod("genome", "IdeogramTrack", function(x, value){
    message("Updating chromosome band information")
    tmp <- IdeogramTrack(genome=value[1], chromosome=chromosome(x), name=names(x))
    displayPars(tmp) <- displayPars(x)
    return(tmp)
})

## Set or extract the name slot of a GdObject
setMethod("names", "GdObject", function(x) x@name)
setReplaceMethod("names", signature("GdObject", "character"), function(x, value) {
    x@name <- value[1]
    return(x)
})

## Set or extract the strand information from a RangeTrack object
setMethod("strand", "RangeTrack", function(x) as.character(strand(ranges(x))))
setMethod("strand", "GenomeAxisTrack", function(x) as.character(strand(ranges(x))))
setMethod("strand", "DataTrack", function(x) x@strand)
setReplaceMethod("strand", "RangeTrack", function(x, value){
    if(length(value)!=1 && length(value)!=length(x))
        stop("Length of replacement value for the strand information does not match the ",
             "number of items in the track")
    r <- ranges(x)
    strand(r) <- value
    x@range <- r
    return(x)
})
setReplaceMethod("strand", "DataTrack", function(x, value){
    if(!is.character(value) && length(value)!=1 && !value %in% c("+", "-", "*"))
        stop("Invalid replacement value")
    x@strand <- value
    return(x)
 })

## Allow for subsetting of RangeTrack and DataTrack objects
setMethod("[", signature(x="RangeTrack"), function(x, i) {
    x <- .deepCopyPars(x)
    x@range <- x@range[i,]
    return(x)})
setMethod("[", signature(x="GenomeAxisTrack"), function(x, i) {
    x <- .deepCopyPars(x)
    x@range <- x@range[i,]
    return(x)})
setMethod("[", signature(x="IdeogramTrack"), function(x, i) return(x))
setMethod("[", signature(x="DataTrack"), function(x, i, j) {
    x <- .deepCopyPars(x)
    if(!missing(i))
    {
        x@data <- x@data[i,, drop=FALSE]
        displayPars(x) <- list(groups=as.vector(displayPars(x, "groups")[i]))
    }
    if(!missing(j))
    {
        x@range <- x@range[j,]
        x@data <- x@data[,j, drop=FALSE]
    }
    return(x)})
setMethod("[", signature(x="AlignedReadTrack"), function(x, i) {
    if(x@coverageOnly)
        stop("This AlignedReadTrack object contains coverage information only and can not be subset")
    x@range <- x@range[i,]
    x <- setCoverage(x)
    return(x)})

## Split a RangeTrack or DataTrack by a factor or character
setMethod("split", signature("RangeTrack"),
          definition=function(x, f, ...){
              rs <- split(ranges(x), factor(f))
              lapply(rs, function(y) {x@range <- y; return(x)})
          })
setMethod("split", signature("AlignedReadTrack"),
          definition=function(x, f, ...){
              if(x@coverageOnly)
                  stop("This AlignedReadTrack object contains coverage information only and can not be split")
              rs <- split(ranges(x), factor(f))
              lapply(rs, function(y) {x@range <- y; x <- setCoverage(x); return(x)})
		  })
setMethod("split", signature("DataTrack"),
		  definition=function(x, f, ...){
			  rs <- as.list(split(ranges(x), factor(f)))
			  ds <- split(t(values(x)), f)
			  nr <- nrow(values(x))
			  mapply(function(y, z) {x@range <- y; x@data <- matrix(z, nrow=nr, byrow=TRUE); return(x)}, rs, ds)
		  })

## Extract the coverage information
setMethod("coverage", signature("AlignedReadTrack"),
		definition=function(x, strand="*"){
			str <- c("+", "-", "*")[.strandName(strand, extended=TRUE)+1]
			return(if(!is.null(x@coverage[[str]])) x@coverage[[str]] else Rle())
		})

##----------------------------------------------------------------------------------------------------------------------------
## There are several levels of annotation information for most RangeTrack objects: individual features (e.g. exons, biotype),
## groups (e.g. transcripts) and even groups of groups (e.g. genes). Not all are relevant for all subclasses, however we want
## to have accessors and replacement methods for a clean interface. 
##----------------------------------------------------------------------------------------------------------------------------
## Helper functions to extract or replace the various annotation data of a track.
##   o GdObject: the input GeneRegionTrack track object
##   o type: the annotation type, i.e., a column in the elementMetadata slot of the GRanges object
##   o value: the replacement value, has to be of the same length as length(GdObject)
.getAnn <- function(GdObject, type) return(as.character(values(GdObject)[[type]]))
.setAnn <-  function(GdObject, value, type)
{
    v <- values(GdObject)
    if(length(value)>1 && length(value) != nrow(v))
        stop("The length of the replacement value for the '", type, "' annotation does not match the number ",
             "of features in the track.")
    v[[type]] <- value
    elementMetadata(GdObject@range) <- v
    return(GdObject)
}

## Accessors to the gene-level annotation of a GeneRegionTrack. We actually need two methods here, one for the
## human-readable gene symbols and one for the actual gene id
setMethod("gene", signature(GdObject="GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "gene"))
setMethod("symbol", signature(GdObject="GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "symbol"))
setReplaceMethod("gene", signature("GeneRegionTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "gene"))
setReplaceMethod("symbol", signature("GeneRegionTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "symbol"))

## Accessors to the transcript-level annotation of a GeneRegionTrack. 
setMethod("transcript", signature(GdObject="GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "transcript"))
setReplaceMethod("transcript", signature("GeneRegionTrack", "character"),
                 function(GdObject, value) .setAnn(GdObject, value, "transcript"))

## Accessors to the exon-level annotation of a GeneRegionTrack
setMethod("exon", signature(GdObject="GeneRegionTrack"), function(GdObject) .getAnn(GdObject, "exon"))
setReplaceMethod("exon", signature("GeneRegionTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "exon"))

## Accessors to the biotype annotation of a RangeTrack.
setMethod("feature", signature(GdObject="RangeTrack"), function(GdObject) .getAnn(GdObject, "feature"))
setMethod("feature", signature(GdObject="DataTrack"), function(GdObject) NULL)
setReplaceMethod("feature", signature("RangeTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "feature"))
setReplaceMethod("feature", signature("DataTrack", "character"), function(GdObject, value) GdObject)

## Accessors to the grouping information of a AnnotationTrack and a GeneRegionTrack (for which this is essentially
## an alias to the transcript accessor.
setMethod("group", "AnnotationTrack", function(GdObject) .getAnn(GdObject, "group"))
setReplaceMethod("group", signature("AnnotationTrack", "character"), function(GdObject, value) .setAnn(GdObject, value, "group"))
setMethod("group", "GeneRegionTrack", function(GdObject) transcript(GdObject))
setReplaceMethod("group", signature("GeneRegionTrack", "character"),
                 function(GdObject, value) .setAnn(GdObject, value, "transcript"))
setMethod("group", "GdObject", function(GdObject) NULL)

## extract or replace the content of the imageMap slot
setMethod("imageMap", "GdObject", function(GdObject) GdObject@imageMap)
setReplaceMethod("imageMap", signature("GdObject", "ImageMapOrNULL"), function(GdObject, value){
    GdObject@imageMap <- value
    return(GdObject)})

## Context-dependent meta-accessors to the identifier data of a AnnotationTrack and a GeneRegionTrack. For the former, those will
## be the content of the id column in the elementMetadata slot of the GRanges object, for the latter, either the gene ids
## (if gpar geneSymbols==FALSE) or the human-readable gene symbols (if gpar geneSymbols==TRUE). If lowest==TRUE the lowest-level
## annotation is returned, i.e., exon ids  for GeneRegionTracks and id for AnnotationTracks. 
setMethod("identifier", "AnnotationTrack", function(GdObject, lowest=FALSE)  if(lowest) .getAnn(GdObject, "id") else group(GdObject))
setMethod("identifier", "GeneRegionTrack", function(GdObject, lowest=FALSE)
          if(lowest) exon(GdObject) else if(.dpOrDefault(GdObject, "geneSymbols", TRUE)) symbol(GdObject) else gene(GdObject))
setReplaceMethod("identifier", c("AnnotationTrack", "character"), function(GdObject, value){
    group(GdObject) <- value
    return(GdObject)})
setReplaceMethod("identifier", c("GeneRegionTrack", "character"), function(GdObject, value){
    if(.dpOrDefault(GdObject, "geneSymbols", TRUE)) symbol(GdObject) <- value else gene(GdObject) <- value
    return(GdObject)})
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## Stacking controls what to do with overlapping annotation regions. 
##----------------------------------------------------------------------------------------------------------------------------
setMethod("stacking", "StackedTrack", function(GdObject) GdObject@stacking)
setReplaceMethod("stacking", c("StackedTrack", "character"), 
                 function(GdObject, value) {
                     pt <- getClass("StackedTrack")@prototype
                     if(!all(value %in% pt@stackingValues))
                         stop("Problem initializing StackedTrack,  need the following values for 'stacking':",
                              paste(pt@stackingValues, collapase=", "), "\n")
                     GdObject@stacking <- value
                     displayPars(GdObject) <- list(stacking=value)
                     return(GdObject)
                 })
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## Recompute or return the stacking information for different types of StackedTrack objects. Stacking in needed when
## annotation regions in the objects overlap and when the stacking type is set to squish, full or pack. Since stacking
## can be dependent on the available space (if feature annotation is added) we need to be able to recompute this before
## we start the actual plotting. For the different sub-classes of StackedTracks we need different behaviour of setStacks:
##   o StackedTrack: there are no groups, so each feature can be treated separately
##   o AnnotationTrack: features can be grouped, in which case we have to avoid overlapping of the whole group region,
##      i.e, from the start of the first group item to the end of the last
##   o GeneRegionTrack: in addition to grouping (in this case transcripts) we have to factor in additional space before
##      each group if transcript annotation is enabled (gpar showId==TRUE). To do so we need to figure out the current
##      fontsize on the device, which means that a device already has to be open and the appropriate viewport has been
##      pushed to the stack. Hence we have to call setStacks immediately before the actual plotting.
## stacks should return a vector of stacks, where each factor level of the vector indicates membeship
## to a particular stacking level.
##----------------------------------------------------------------------------------------------------------------------------
setMethod("stacks", "StackedTrack", 
          function(GdObject) if(length(GdObject@stacks)) GdObject@stacks else NULL)

setMethod("setStacks", "GdObject", function(GdObject, ...) GdObject)
setMethod("setStacks", "StackedTrack", function(GdObject, ...) {
    bins <- if(!.needsStacking(GdObject)) rep(1, length(GdObject)) else disjointBins(range(GdObject))
    GdObject@stacks <- bins
    return(GdObject)
})
setMethod("setStacks", "AnnotationTrack", function(GdObject, from, to) {
    if(!.needsStacking(GdObject))
    {
        bins <- rep(1, length(GdObject))
    } else {
        ids <- identifier(GdObject, FALSE)
        hasAnno <- .dpOrDefault(GdObject, "showId", FALSE) & ids!=""
        uid <- make.unique(identifier(GdObject, lowest=TRUE))
        gp <- group(GdObject)
        needsGrp <- any(duplicated(gp))
        if(needsGrp){
            groups <- split(range(GdObject), gp)
            uidSplit <- split(uid, gp)
            ## FIXME: This is the rate-limitting step (about 60% of the time spent for the whole track plotting), would be nice to speed this up somehow...
            groupRanges <- sapply(groups, function(x) c(from=min(start(x)), to=max(end(x)), len=length(x)))
        }
        if(any(hasAnno))
        {
            txt <- paste(ids, " ")
            txt[!hasAnno] <- ""
            identifier(GdObject) <- txt
            cex <- .dpOrDefault(GdObject, "cex", 1) * .dpOrDefault(GdObject, "cex.symbol", 0.7)
            fontfamily <- .dpOrDefault(GdObject, "fontfamily", 1)
            fontsize <- .dpOrDefault(GdObject, "fontsize", 12)
            fontface <- .dpOrDefault(GdObject, "fontface.symbol", 2)
            pushViewport(dataViewport(xData=c(from, to), extension=0, yscale=c(0, 40), clip=TRUE,
                                      gp=gpar(cex=cex, fontfamily=fontfamily, fonface=fontface, fontsize=fontsize)))
            if(needsGrp)
            {
                ids <- sapply(split(identifier(GdObject), gp), function(x) x[1])
                gRanges <- IRanges(start=groupRanges["from",]-(as.numeric(convertWidth(stringWidth(ids),"native"))*1.3),
                                   end=groupRanges["to",], names=names(groups))
            } else {
                gRanges <- range(GdObject)
                start(gRanges) <- start(gRanges)-(as.numeric(convertWidth(stringWidth(identifier(GdObject)),"native"))*1.3)
            }
            popViewport(1)
        } else {
            gRanges <- if(needsGrp){ if(length(groups)) IRanges(start=groupRanges["from",], end=groupRanges["to",],
                                                                names=names(groups)) else IRanges() } else {
                                                                    range(GdObject)}
        }
        bins <- if(needsGrp) rep(disjointBins(gRanges), groupRanges["len",]) else disjointBins(gRanges)
        names(bins) <- if(needsGrp) unlist(uidSplit) else uid
        bins <- bins[uid]
    }
    bins <- if(length(bins)) bins else 0
    GdObject@stacks <- bins
    return(GdObject)
})
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## In some cases we don't need a range but rather a single position. Most of the time this is simply taking the geometric
## mean of the range. If numeric values are associated to positions we have to be able to extract those as well.
##----------------------------------------------------------------------------------------------------------------------------
## The geometric mean of annotation ranges
setMethod("position", signature("RangeTrack"), definition=function(GdObject, from=NULL, to=NULL, sort=FALSE)
      {
          GdObject <- subset(GdObject, from=from, to=to, sort=sort)
          return(apply(cbind(start(GdObject), end(GdObject)), 1, mean))
      })
setMethod("position", signature("IdeogramTrack"), definition=function(GdObject, ...) NULL)
     
## The numeric values of data tracks
setMethod("score", signature("DataTrack"), function(x, from=NULL, to=NULL, sort=FALSE, transformation=TRUE)
      {
          x <- subset(x, from=from, to=to, sort=sort)
          vals <- values(x)
          ## apply data transformation if one is set up
          trans <- .dpOrDefault(x, "transformation")
          if(is.list(trans))
              trans <- trans[[1]]
          if(transformation && !is.null(trans))
          {
              if(!is.function(trans) || length(formals(trans))!=1L)
                  stop("gpar 'transformation' must be a function with a single argument")
              test <- trans(vals)
              if(!is.numeric(test) || !is.matrix(test) || !all(dim(test) == dim(vals)))
                  stop("The function in gpar 'transformation' results in invalid output.\n",
                       "It has to return a numeric matrix with the same dimensions as the input data.")
              vals <- test
          }
          return(vals)
      })



smartCollapse <- function(GdObject)
{
    collapse <- .dpOrDefault(GdObject, "collapse", TRUE)
    min.width <- .dpOrDefault(GdObject, "min.width", 2)
    min.distance <- .dpOrDefault(GdObject, "min.distance", 2)
    diff <- .pxResolution(coord="x")
    r <- ranges(GdObject)
    r <- .resize(r, min.width, diff)
    cv <- coverage(r)[[chromosome(GdObject)]]
    dt <- DataTrack(start=start(cv), end=end(cv), data=runValue(cv), chromosome=chromosome(GdObject),
                    genome=genome(GdObject), name=names(GdObject))
}





## There is a natural limit of what can be plotted as individual features caused by the maximum resolution of the device.
## Essentially no object can be smaller than the equivalent of a single pixel on the screen or whatever a pixel corresponds
## to on other devices.
## Thus we need to adapt the all objects in the track to the current resolution. Things that are closer than a certain limit
## distance can be collapsed (if that is desired), things that are smaller than 1 pixel can be blown up to a minimum size.
## All collapseTrack methods should always return a GdObject instance of similar type as the input to allow us to keep the
## collpasing step optional and for all the downstream plotting operations to still work. The internal (mostly the GRanges)
## objects can be modified to whatever is desired. Please note that this method is called after the drawing viewport has been
## pushed, and hence all the coordinate systems are already in place.
## Available arguments are:
##    o GdObject: the input GeneRegionTrack track object
##    o min.width: the minimum width in pixel, everything that's smaller will be expanded to this size.
##    o min.distance: the minimum distance between two features in pixels below which to start collapsing track items.
##    o collapse: logical, collapse overlapping items into a single meta-item.
##    o diff: the equivalent of 1 pixel in the native coordinate system.
##    o xrange: the data range on the x axis. Can be used for some preliminary subsetting to speed things up
##----------------------------------------------------------------------------------------------------------------------------
## For AnnotationTracks we need to collapse the all regions along with the additional annotation.
## For GeneRegionTracks we essentially need to do the same thing as for AnnotationTracks, however the additional annotation columns
## are quite different. If stacking is turned on
## we don't need to collapse across group boundaries. Those will anyways be separated in different rows. However, since
## collapseTrack is called after the stacks have been split we don't have to deal with it here.
setMethod("collapseTrack", signature(GdObject="AnnotationTrack"),        
          function(GdObject, diff=.pxResolution(coord="x"), xrange) {
              collapse <- .dpOrDefault(GdObject, "collapse", TRUE)
              min.width <- .dpOrDefault(GdObject, "min.width", 2)
              min.distance <- .dpOrDefault(GdObject, "min.distance", 2)
              GdObject <- GdObject[order(range(GdObject))]
              r <- ranges(GdObject)
              ## Collapse overlapping ranges (less than minXDist space between them) including the associated attributes using
              ## "|" as separator. For both "strand" and "feature" we take the first available entry, which is not optimal but
              ## seems to be the sanest thing to do here...
              if(collapse)
              {
                  minXDist <- min.distance*diff
                  rr <- ranges(r)
                  if(minXDist<1)
                  {
                      ## We have to fake smaller ranges because reduce will merge also neigbouring ranges
                      width(rr) <- width(rr)-1
                      rr <- reduce(rr, min.gapwidth=minXDist)
                      width(rr) <- width(rr)+1
                  } else {
                      rr <- reduce(rr, min.gapwidth=minXDist)
                  }
                  if(length(rr) < length(r))
                  {
                      startInd <- sort(unique(sapply(start(rr), function(x) min(which(start(r)==x)))))
                      startInd <-  if(tail(startInd,1) == length(r)) c(startInd, length(r)+1) else c(startInd, length(r))
                      vsplit <- split(cbind(values(GdObject), strand=I(as.character(strand(GdObject)))),
                                      cut(seq_len(length(r)), startInd, iclude.lowest=TRUE, right=FALSE))
                      if(is(GdObject, "GeneRegionTrack"))
                      {
                          newValsId <- as.data.frame(t(sapply(vsplit, function(x)
                                                              c(gene=as.character(paste(paste(gsub("  $", "", unique(x$gene)), collapse="|"),
                                                                                        "  ", sep="")),
                                                                transcript=as.character(paste(unique(x$transcript), collapse="|")),
                                                                symbol=as.character(paste(paste(gsub("  $", "", unique(x$symbol)), collapse="|"),
                                                                                                "  ", sep="")),                     
                                                                exon=as.character(paste(unique(x$exon), collapse="|"))))),
                                                     stringsAsFactors=FALSE)
                          rcnams <- c(setdiff(colnames(values(GdObject)), c("gene", "transcript", "symbol", "exon")), "strand") 
                      } else {
                          newValsId <- as.data.frame(t(sapply(vsplit, function(x)
                                                              c(id=as.character(paste(unique(x$id), collapse="|")),
                                                                group=as.character(paste(paste(gsub("  $", "", unique(x$group)), collapse="|"),
                                                                                         "  ", sep=""))))),
                                                     stringsAsFactors=FALSE)
                          rcnams <- c(setdiff(colnames(values(GdObject)), c("id", "group")), "strand")
                      }
                      newValsRest <- t(sapply(vsplit, function(x) sapply(x[1,rcnams], as.character)))
                      colnames(newValsRest) <- rcnams
                      newVals <- data.frame(newValsId, newValsRest, stringsAsFactors=FALSE)
                      str <- newVals$strand
                      newVals <- newVals[,colnames(values(GdObject))]
                      if(.dpOrDefault(GdObject, "showOverplotting", FALSE))
                      {
                          ## Calculate colors indicating the amount of overplotting
                          dens <- sapply(vsplit, nrow)
                          if(length(dens)>1 && any(dens>1) && length(unique(dens))!=1)
                          {
                              #minSat <- max(0.25, 1/max(dens))
                              #saturation <- minSat+((dens-min(dens))/diff(range(dens))/(1/(1-minSat)))
                              #baseCol <- rgb2hsv(col2rgb(.dpOrDefault(GdObject, "col.overplotting",
                              #                                        Gviz:::.DEFAULT_OVERPLOT_COL)))
                              #desatCols <- sapply(saturation, function(x) hsv(baseCol[1], x, baseCol[3]))
                              #newVals$feature[dens>1] <- paste("__desatCol", dens[dens>1], sep="")
                              #names(desatCols) <- newVals$feature
                              #displayPars(GdObject) <- as.list(desatCols[dens>1])
                              newVals[, "density"] <- dens
                          }
                      }
                      r <- GRanges(seqnames=chromosome(GdObject), strand=str, range=rr)
                      elementMetadata(r)<- newVals
                  }
              }
              ## Compute native coordinate equivalent to 1 pixel and resize
              r <- .resize(r, min.width, diff)
              ## Reconstuct the RangedData object and return
              GdObject@range <- r
              return(GdObject)})


.aggregator <- function(GdObject)
{
    agFun <- .dpOrDefault(GdObject, "aggregation", "mean")
    if(is.list(agFun))
        agFun <- agFun[[1]]
    fun <- if(is.character(agFun))
    {
        switch(agFun,
               "mean"=rowMeans,
               "sum"=rowSums,
               "median"=rowMedians,
               "extreme"=function(x) apply(x, 1, .extreme),
               "min"=rowMin,
               "max"=rowMax,
               rowMeans)
    } else {
        if(is.function(agFun)){
            function(x) apply(x, 1, agFun)
        }else stop("display parameter 'aggregation' has to be a function or a character",
                   "scalar in c('mean', 'median', 'sum', 'extreme')")
    }
    return(fun)
}



## For DataTracks we want to collapse data values using the aggregation function provided by calling .aggregator().
## In addition values can be aggregated over fixed window slices when gpar 'window' is not NULL, and using a sliding
## window approach when 'window' == -1
setMethod("collapseTrack", signature(GdObject="DataTrack"), function(GdObject, diff=.pxResolution(coord="x"), xrange) {
    if(!length(GdObject))
        return(GdObject)
    ## first the data transformation if needed
    values(GdObject) <- score(GdObject)
    collapse <- .dpOrDefault(GdObject, "collapse", FALSE)
    min.width <- .dpOrDefault(GdObject, "min.width", 2)
    min.distance <- max(0, .dpOrDefault(GdObject, "min.distance", 0))
    ## When an averaging window has been set, split the data up into these average chunks
    window <- .dpOrDefault(GdObject, "window", NULL)
    windowSize <- .dpOrDefault(GdObject, "windowSize", NULL)
    if(!is.null(window) || collapse)
        GdObject <- GdObject[,order(range(GdObject))]
    r <- ranges(GdObject)
    drange <- c(floor(xrange[1]), ceiling(xrange[2]))
    if(!is.null(window))
    {
        rr <-  if(is(r, "GRanges")) ranges(r) else r
        fw <- FALSE
        if(window=="auto")
            window <- min(ncol(values(GdObject)), 1000,
                          ceiling(width(range(rr))/(min.width*diff)))
        if(window=="fixed"){
            fw <- TRUE
            window <- 100
        }
        if(!is.numeric(window) || length(window)!=1L)
            stop("gpar 'window' must be a numeric scalar")
        window <- as.integer(window)
        sc <- values(GdObject)
        agFun <- .aggregator(GdObject)
        if(window==1)
        {
            sc <- matrix(agFun(sc), ncol=1)
            rtmp  <- IRanges(start=max(1,drange[1]), end=max(1, drange[2]-1))
            r <- if(is(r, "GRanges"))  GRanges(seqnames=seqnames(r)[1], range=rtmp) else rtmp
        } else if(window<1){
            if(is.null(windowSize))
                windowSize <- (max(GdObject)-min(GdObject))/100
            if(windowSize %% 2 !=1)
                windowSize <- windowSize+1
            rm <- vector("integer", width(range(range(GdObject))))
            ind <- unlist(mapply(function(x, y) x:y, start(GdObject), end(GdObject)))-min(GdObject)+1
            rm[ind] <- rep(sc[1,], width(GdObject))
            runwin <- suppressWarnings(runmean(Rle(as.numeric(rm)), k=windowSize, endrule="constant"))
            seqSel <- findRun(as.integer(position(GdObject))-min(GdObject)+1, runwin)
            newDat <- matrix(runValue(runwin)[seqSel], nrow=1)
            if(nrow(sc)>1)
            {
                newDat <- rbind(newDat, 
                	matrix(sapply(2:nrow(sc), function(x) {
                    	rm[ind] <- rep(sc[x,], width(GdObject))
                    	suppressWarnings(runValue(runmean(Rle(as.numeric(rm)), k=windowSize, endrule="constant")))[seqSel]}), 
                    	nrow=nrow(sc)-1, byrow=TRUE)
            	)
            }
            sc <- newDat
        } else {
            if(!is.null(window) && window > diff(drange))
                window <- diff(drange)
            if(!fw || is.null(windowSize)){
                windowSize <-  diff(drange) %/% window
            }else{
                window <- max(1, diff(drange) %/% windowSize)
            }
            remain <-  (diff(drange) - (window * windowSize))/2
            ir <- IRanges(start=seq(from=drange[1]+remain, to=drange[2]-remain-windowSize, length.out=window),
                          width=windowSize)
            if(remain>0)
                ir <- c(IRanges(start=drange[1], width=ceiling(remain)), ir,
                        IRanges(start=drange[2]-ceiling(remain), width=ceiling(remain)))
            ol <- as.matrix(findOverlaps(ir, rr))
            scn <- sapply(split(ol[,2], ol[,1]), function(i) agFun(sc[,i,drop=FALSE]), USE.NAMES=FALSE)
            if(is.null(dim(scn)))
                scn <- matrix(scn, nrow=nrow(values(GdObject)), dimnames=list(NULL, as.character(unique(ol[,1]))))
            sc <- matrix(NA, ncol=length(ir), nrow=nrow(scn))
            sc[, as.integer(colnames(scn))] <- scn
            r <- if(is(r, "GRanges")) GRanges(seqnames=chromosome(GdObject), range=ir,
                                              strand=unique(as.character(strand(GdObject)))) else ir
        }
        GdObject@range <- r
        GdObject@data <- sc
    } 
    ## Compute native coordinate equivalent to 1 pixel and resize
    r <- .resize(r, min.width, diff)
    ## Collapse overlapping ranges (less than minXDist space between them) including the associated attributes using
    ## "|" as separator. For both "strand" and "feature" we take the first available entry, which is not optimal but
    ## seems to be the sanest thing to do here...
    if(collapse)
    {
        minXDist <- min.distance*diff
        rr <- if(is(r, "GRanges")) ranges(r) else r
        if(minXDist<1)
        {
            ## We have to fake smaller ranges because reduce will merge also neigbouring ranges
            width(rr) <- width(rr)-1
            rr <- reduce(rr, min.gapwidth=minXDist)
            width(rr) <- width(rr)+1
        } else {
            rr <- reduce(r, min.gapwidth=minXDist)
        }
        sc <- values(GdObject)
        if(length(rr)==1){
            r <- GRanges(seqnames=1, strand=strand(GdObject)[1], range=rr)
            GdObject@range <- r
            GdObject@data <- matrix(rowMeans(sc, na.rm=TRUE), ncol=1)
        } else if(length(rr) < length(r)){
            startInd <- sort(unique(sapply(start(rr), function(x) which(start(r)==x))))
            st <- strand(GdObject)
            startInd <- if(tail(startInd,1) == length(r)) c(startInd, length(r)+1) else c(startInd, length(r))
            vsplit <- split(t(as.data.frame(sc, stringsAsFactors=FALSE)), cut(seq_len(length(r)), startInd, iclude.lowest=TRUE, right=FALSE))
            agFun <- .dpOrDefault(GdObject, "aggregation", "mean")
            if(is.list(agFun))
                agFun <- agFun[[1]]
            newScore <- if(is.character(agFun)){
                switch(agFun, "mean"=sapply(vsplit, function(x) rowMeans(matrix(x, nrow=nrow(sc), byrow=TRUE), na.rm=TRUE), USE.NAMES=FALSE),
                       "sum"=sapply(vsplit, function(x) rowSums(matrix(x, nrow=nrow(sc), byrow=TRUE), na.rm=TRUE), USE.NAMES=FALSE),
                       "median"=sapply(vsplit, function(x) rowMedians(matrix(x, nrow=nrow(sc), byrow=TRUE), na.rm=TRUE), USE.NAMES=FALSE),
                       sapply(vsplit, function(x) rowMeans(matrix(x, nrow=nrow(sc), byrow=TRUE), na.rm=TRUE), USE.NAMES=FALSE))
            } else {
                if(is.function(agFun)){
                    sapply(vsplit, function(x) apply(matrix(x, nrow=nrow(sc), byrow=TRUE), 1, function(y) agFun(y)[1]), USE.NAMES=FALSE)
                } else stop("display parameter 'aggregation' has to be a function or a character ", "scalar in c('mean', 'median', 'sum')")
            }
            r <- GRanges(seqnames=seq_len(length(rr)), strand=st, range=rr)
            GdObject@data <- newScore
            GdObject@range <- r
        }
    }
    ## Reconstruct the RangedData object and return
    GdObject@range <- r
    return(GdObject)
})



## For a GenomeAxisTrack all we need to do is collapse the optional ranges
setMethod("collapseTrack", signature(GdObject="GenomeAxisTrack"), 
          function(GdObject, min.width=1, min.distance=0, collapse=TRUE, diff=.pxResolution(coord="x"), xrange) {
             
              ## Collapse overlapping ranges (less than minXDist space between them) including the associated attributes using
              ## "|" as separator. For both "strand" and "feature" we take the first available entry, which is not optimal but
              ## seems to be the sanest thing to do here...
              if(collapse)
              {
                  GdObject <- GdObject[order(range(GdObject))]
                  r <- ranges(GdObject)
                  minXDist <- min.distance*diff
                  r <- reduce(r, min.gapwidth=minXDist)
              }
              r <- .resize(r, min.width, diff)
              GdObject@range <- r
              return(GdObject)})
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## Truncate a GdObject and sort by coordinates if necessary.
##----------------------------------------------------------------------------------------------------------------------------
## The default is not to clip at all
setMethod("subset", signature(x="GdObject"), function(x, ...) x)
setMethod("subset", signature(x="IdeogramTrack"), function(x, ...) x)
## For normal ranges we clip everything outside of the boundaries (keeping one extra item left and right
## in order to assure continuation)
setMethod("subset", signature(x="RangeTrack"), function(x, from=NULL, to=NULL, sort=FALSE, drop=TRUE, ...){
    ## Not needed anymore...
    ## Subset to a single chromosome first
    if(drop){
        csel <- seqnames(x) != chromosome(x)
        if(any(csel))
            x <- x[,!csel]
    }
    if(!length(x))
        return(x)
    ranges <- .defaultRange(x, from=from, to=to)
    lsel <- end(x) < ranges["from"]
    if(any(lsel))
        lsel[max(0, max(which(lsel))-1)] <- FALSE   
    rsel <- start(x) > ranges["to"]
    if(any(rsel))
        rsel[min(length(x), min(which(rsel))+1)] <- FALSE
    if(any(lsel) || any(rsel))
        x <- x[!(lsel | rsel),]
    if(sort)
        x <- x[order(range(x)),]
    return(x)
})
## For data tracks we cut exactly, and also reduce to the current chromosome unless told explicitely not to
setMethod("subset", signature(x="DataTrack"), function(x, from=NULL, to=NULL, sort=FALSE, drop=TRUE, ...){
    ## Subset to a single chromosome first
    if(drop){
        csel <- seqnames(x) != chromosome(x)
        if(any(csel))
            x <- x[,!csel]
    }
    if(!length(x))
        return(x)
    ranges <- .defaultRange(x, from=from, to=to)
    x <- x[,start(x)>=ranges["from"] & end(x)<=ranges["to"]]
    if(sort)
        x <- x[,order(range(x))]
    return(x)
})                    
## Only recompute the stacks here
setMethod("subset", signature(x="StackedTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE){
    x <- callNextMethod(x=x, from=from, to=to, sort=sort)
    if(stacks)
        x <- setStacks(x)
    return(x)
})
## In order to keep the grouping information for track regions in the clipped areas we also have to
## keep, for each group, the min-1 and max+1 items.
setMethod("subset", signature(x="AnnotationTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE){
    ## Subset to a single chromosome first
    csel <- seqnames(x) != chromosome(x)
    if(any(csel))
        x <- x[!csel]
    if(length(x))
    {
        ## Nothing to do if everything is within the range
        ranges <- .defaultRange(x, from=from, to=to)
        if(!any(end(x)<ranges["from"] | start(x)> ranges["to"])){
            if(stacks)
                x <- setStacks(x)
            return(x)
        }
        
        ## Now remove everything except for the overlapping groups by first subselecting all groups in the range...
        gsel <- group(x)[queryHits(findOverlaps(range(x), IRanges(ranges["from"], ranges["to"])))]
        x <- x[group(x) %in% gsel]
        ## check if there is anything that overlaps
        if (length(gsel)) {
          ## ... and then finding everything that overlaps
          coords <- data.frame(cbind(start(x), end(x)))
          grps <- sapply(split(coords, group(x)), range)
          gsel <- group(x) %in% colnames(grps)[subjectHits(findOverlaps(IRanges(ranges["from"], ranges["to"]),
                                                                        IRanges(grps[1,], grps[2,])))]
          x <- x[gsel]
          ## Now only keep the adjancend items if there are any.
          lsel <- end(x) < ranges["from"]
          rsel <- start(x) > ranges["to"]
          ## Nothing to do if everything is within the range
          if(any(lsel | rsel))
            {
              grp <- group(x)
              if(any(table(grp)>1))
                {
                  lsel[lsel][unlist(sapply(split(seq_len(sum(lsel)), grp[lsel]), max))] <- FALSE
                  rsel[rsel][unlist(sapply(split(seq_len(sum(rsel)), grp[rsel]), min))] <- FALSE
                }
              x <- x[!(lsel | rsel),]
            }
          if(sort)
            x <- x[order(range(x)),]
          if(stacks)
            x <- setStacks(x)
        }
      }
    return(x)
})
## For the axis track we may have to clip the highlight ranges on the axis.
setMethod("subset", signature(x="GenomeAxisTrack"), function(x, from=NULL, to=NULL, sort=FALSE, ...){
    if(!length(x))
        return(x)
    ranges <- .defaultRange(x, from=from, to=to)
    lsel <- end(x) < ranges["from"] 
    rsel <- start(x) > ranges["to"]
    x <- x[!(lsel | rsel),]
    if(sort)
        x <- x[order(range(x)),]
    return(x)
})
## If the object only stores coverage we subset that, otherwise we can use the RangeTrack method
setMethod("subset", signature(x="AlignedReadTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE){
    if(x@coverageOnly)
    {
        if(is.null(from))
            from <- min(unlist(lapply(x@coverage, function(y) if(length(y)) min(start(y)))))
        if(is.null(to))
            to <- max(unlist(lapply(x@coverage, function(y) if(length(y)) max(start(y)))))
        x@coverage <- lapply(x@coverage, function(y){runValue(y)[end(y)<from | start(y)>to] <- 0; y})
        from <- min(unlist(lapply(x@coverage, function(y) if(length(y)) head(start(y),2)[2])))
        to <- max(unlist(lapply(x@coverage, function(y) if(length(y)) tail(end(y),2)[1])))
        x@range <- GRanges(range=IRanges(start=from, end=to), strand="*", seqnames=1:2)
    }else{
        x <- callNextMethod(x=x, from=from, to=to, sort=sort, stacks=stacks)
    }
    return(x)
})
##----------------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------------
## The indivdual bits and pieces of a Gviz plot are all drawn by separate renderers. Currently, those are a y-axis,
## a grid, and the actual track panel.
##----------------------------------------------------------------------------------------------------------------------------
## For certain GdObject subclasses we may want to draw a y-axis. For others an axis is meaningless, and the default function
## will return NULL without plotting anything.
setMethod("drawAxis", signature(GdObject="GdObject"), function(GdObject, ...) return(NULL))
setMethod("drawAxis", signature(GdObject="DataTrack"), function(GdObject, ...) {
    type <- match.arg(.dpOrDefault(GdObject, "type", "p"), c("p", "l", "b", "a", "s", "g", "r", "S", "smooth",
                                                            "histogram", "mountain", "h", "boxplot", "gradient", "heatmap"),
                      several.ok=TRUE)
    if(.dpOrDefault(GdObject, "legend", FALSE)){
         pushViewport(viewport(y=1, height=unit(1, "npc") - unit(getPar(GdObject, ".__verticalSpace"), "inches"),
                              just=c(0.5, 1)))
         on.exit(popViewport(1))
    }
    if(length(type)==1L && (type=="gradient" || type=="heatmap")) return(NULL) else callNextMethod()
})
setMethod("drawAxis", signature(GdObject="NumericTrack"), function(GdObject, from, to, ...) {
    yvals <- values(GdObject)
    ylim <- .dpOrDefault(GdObject, "ylim", if(!is.null(yvals) && length(yvals)) 
    	 range(yvals, na.rm=TRUE, finite=TRUE) else c(-1,1))
    if(diff(ylim)==0)
        ylim <- ylim+c(-1,1) 
    yscale <- extendrange(r=ylim, f=0.05)
    vpTitleAxis <- viewport(x=0.95, width=0, yscale=yscale, just=0)
    hSpaceAvail <- vpLocation()$isize["width"]/6
    pushViewport(vpTitleAxis)
    col <- .dpOrDefault(GdObject, "col.axis", "white")
    acex <- .dpOrDefault(GdObject, "cex.axis", NULL)
    acol <- .dpOrDefault(GdObject, "col.axis", "white")
    at <- pretty(yscale)
    at <- at[at>=sort(ylim)[1] & at<=sort(ylim)[2]]
    if(is.null(acex))
    {
        vSpaceNeeded <- max(as.numeric(convertWidth(stringHeight(at), "inches")))*length(at)*1.5
        hSpaceNeeded <- max(as.numeric(convertWidth(stringWidth(at), "inches")))
        vSpaceAvail <- abs(diff(range(at)))/abs(diff(yscale))*vpLocation()$isize["height"]
        acex <- max(0.6, min(vSpaceAvail/vSpaceNeeded, hSpaceAvail/hSpaceNeeded))
    }
    suppressWarnings(grid.yaxis(gp=gpar(col=acol, cex=acex), at=at))
    grid.lines(x=c(0,0), y=ylim, gp=gpar(col=acol), default.units="native")
    popViewport(1)
})
setMethod("drawAxis", signature(GdObject="AlignedReadTrack"), function(GdObject, from, to, subset=TRUE) {
    detail <- match.arg(.dpOrDefault(GdObject, "detail", "coverage"), c("coverage", "reads"))
    if(detail!="coverage") return(NULL) else {
        if(subset)
            GdObject <- subset(GdObject, from=from, to=to)
        cov <- coverage(GdObject, strand="*")
        val <- runValue(coverage(GdObject, strand="*"))
        ## We have to figure out the data range, taking transformation into account
        ylim <- .dpOrDefault(GdObject, "ylim") 
        if(is.null(ylim))
        {
            if(!length(val))
                ylim=c(0,1) else{
                    ylim <- c(0, range(val, finite=TRUE, na.rm=TRUE)[2])
                    trans <- displayPars(GdObject, "transformation")[[1]]
                    if(!is.null(trans))
                        ylim <- c(0, trans(ylim[2]))
                }
        }	
        for(s in c("+", "-"))
        {
            pushViewport(viewport(height=0.5, y=ifelse(s=="-", 0, 0.5), just=c("center", "bottom")))
            dummy <- DataTrack(start=rep(mean(c(from, to)),2), end=rep(mean(c(from, to)),2), data=ylim,
                              genome=genome(GdObject), chromosome=chromosome(GdObject))
            oldDp <- displayPars(GdObject)
            oldDp[["ylim"]] <- if(s=="+") ylim else rev(ylim)
            displayPars(dummy) <- oldDp
            drawAxis(dummy, from=from, to=to)
            popViewport(1)
        }
    }
})

## Draw a grid in the background of a GdObject. For some subclasses this is meaningless, and the default function will
## return NULL without plotting anything.
setMethod("drawGrid", signature(GdObject="GdObject"), function(GdObject, ...) return(NULL))
setMethod("drawGrid", signature(GdObject="NumericTrack"), function(GdObject, from, to){
    if(.dpOrDefault(GdObject, "grid", FALSE)) {
        vals <- score(GdObject)
        ylim <- .dpOrDefault(GdObject, "ylim", range(vals, na.rm=TRUE, finite=TRUE))
        if(diff(ylim))
        {
            pushViewport(dataViewport(xData=c(from, to), yData=ylim, extension=c(0, 0.1), clip=TRUE)) 
            panel.grid(h=.dpOrDefault(GdObject, "h", -1), v=.dpOrDefault(GdObject, "v", -1),
                       col=.dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty=.dpOrDefault(GdObject, "lty.grid", 1),
                       lwd=.dpOrDefault(GdObject, "lwd.grid", 1))
            popViewport(1)
        }
    }})
setMethod("drawGrid", signature(GdObject="AnnotationTrack"), function(GdObject, from, to){
    if(.dpOrDefault(GdObject, "grid", FALSE)) {
        pushViewport(dataViewport(xData=c(from, to), extension=c(0, 0), yData=0:1, clip=TRUE)) 
        panel.grid(h=0, v=.dpOrDefault(GdObject, "v", -1),
                   col=.dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty=.dpOrDefault(GdObject, "lty.grid", 1),
                   lwd=.dpOrDefault(GdObject, "lwd.grid", 1))
        popViewport(1)
    }})
setMethod("drawGrid", signature(GdObject="AlignedReadTrack"), function(GdObject, from, to) {
			detail <- match.arg(.dpOrDefault(GdObject, "detail", "coverage"), c("coverage", "reads"))
			if(detail=="coverage"){
                            GdObject <- subset(GdObject, from=from, to=to)
                            cov <- coverage(GdObject, strand="*")
                            val <- if(length(cov)) runValue(cov) else 1
                            ## We have to figure out the data range, taking transformation into account
                            ylim <- .dpOrDefault(GdObject, "ylim") 
                            if(is.null(ylim))
                            {
                                ylim <- c(0, range(val, finite=TRUE, na.rm=TRUE)[2])
                                trans <- displayPars(GdObject, "transformation")[[1]]
                                if(!is.null(trans))
                                    ylim <- c(0, trans(ylim[2]))
                            }	
                            for(s in c("+", "-"))
                            {
                                pushViewport(viewport(height=0.5, y=ifelse(s=="-", 0, 0.5), just=c("center", "bottom")))
                                dummy <- DataTrack(start=rep(mean(c(from, to)),2), end=rep(mean(c(from, to)),2), data=ylim,
                                                   genome=genome(GdObject), chromosome=chromosome(GdObject))
                                oldDp <- displayPars(GdObject)
                                oldDp[["ylim"]] <- if(s=="+") ylim else rev(ylim)
                                displayPars(dummy) <- oldDp
                                drawGrid(dummy, from=from, to=to)
                                popViewport(1)
                            }
			} 
                        return(NULL)
		})
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## All the drawGD methods should support two modes, triggered by the boolean argument 'prepare':
##    In prepare mode: nothing is plotted but the object is prepared for plotting bases on the available space. The return
##       value of the method in this mode should always be the updated object. If nothing needs to be prepared, i.e., if the
##       plotting is independent from the available space, simply return the original object
##    In plotting mode: the object is plotted. Return value is the nout object with optional HTML image map information
##    added to the imageMap slot
## Since subsetting can be potentially expensive when the data are large we want to minimize this operation. Essentially it
## should be done once right at the start, in wich case the subsetting argument can be set to FALSE
##----------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------
## Although the stacking type is not stored as a displayParameter we still want to check whether it is
## included there and set the actual stacking of the object accordingly
setMethod("drawGD", signature("StackedTrack"), function(GdObject, ...){
    st <- .dpOrDefault(GdObject, "stacking")
    if(!is.null(st))
        stacking(GdObject) <- st
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## Draw gene models as found in AnnotationTracks or GeneRegionTracks
##----------------------------------------------------------------------------------------------------------------------------
## Calculate all coordinates and values for the individual stacks first, and append to the
## input list 'toDraw'. This allows us to plot the whole track at once, making use of grid's
## vectorization. Without this tweak, every stack would have to be drawn individually, which
## can be painfully slow...
## Because all of the y coordinates are calculated in the coordinate system of the current
## viewport for the respective stack we have to add offset values.
.prepareGene <- function(GdObject, toDraw=list(), offset=0, minBase, maxBase, diff=.pxResolution(coord="x"), ...)
{
    ## Usually we only proceed if there is something to draw, otherwise we want to return
    ## the unmodified input list items, hence we rbind NULL values
    td.bar <- td.bartext <- NULL
    if(!length(GdObject) || stacking(GdObject) == "hide")
        return(toDraw)
    if(missing(minBase))
        minBase <- min(start(GdObject))-10
    if(missing(maxBase))
        maxBase <- max(end(GdObject))+10
    ylim <- c(0, 1)
    middle <- mean(ylim)
    space <- diff(ylim)/8
    shape <- .dpOrDefault(GdObject, "shape", "arrow")
    ## The 'exon' structure first
    color <- .getBiotypeColor(GdObject)
    boxes <- data.frame(cx1=start(GdObject), cy1=ylim[1]+space+offset, cx2=end(GdObject), cy2=ylim[2]-space+offset,
                        fill=color, strand=strand(GdObject), text=identifier(GdObject, lowest=TRUE),
                        textX=(start(GdObject)+end(GdObject))/2, textY=middle+offset,
                        ##.getImageMap(cbind(start(GdObject), middle-(space*3)+offset, end(GdObject), middle+(space*3)+offset)),
                        .getImageMap(cbind(start(GdObject), ylim[1]+space+offset, end(GdObject), ylim[2]-space+offset)),
                        start=start(GdObject), end=end(GdObject), values(GdObject), stringsAsFactors=FALSE)
    rownames(boxes) <- make.unique(identifier(GdObject, TRUE))
    ## Combine grouped elements for connecting bars and optional text
    gp <- gsub("  $", "", group(GdObject))
    gTab <- table(gp)
    gSize <- gTab>1
    needsGrp <- any(gSize)
    sel <- gp %in% names(which(gSize))
    ## Prepare a data.frame to hold the group annotation if needed
    if(.dpOrDefault(GdObject, "showId", FALSE)){
        td.bartext <- as.data.frame(matrix(ncol=3, nrow=length(gTab),
                                           dimnames=list(names(gTab), c("txt", "x", "y"))),
                                    stringsAsFactors=FALSE)
        td.bartext$y <- middle+offset
        if(sum(!gSize)>0)
            td.bartext[names(which(!gSize)), c("txt", "x")] <- data.frame(identifier(GdObject)[!sel], start(GdObject)[!sel],
                                                                          stringsAsFactors=FALSE)
    }
    if(needsGrp){
        
        GdRanges <- ranges(GdObject[sel])
        values(GdRanges)$color <- .getBiotypeColor(GdObject)[sel]
        values(GdRanges)$identifier <- identifier(GdObject)[sel]
        groups <- as.list(split(GdRanges, gp[sel]))
        nrs <- ifelse(any(c("arrow", "smallArrow") %in% shape), sum(gSize), sum(pmax(2, gTab[gSize]-1))*2)
        td.bar <- as.data.frame(matrix(ncol=5, nrow=nrs,
                                       dimnames=list(NULL, c("sx1", "sx2", "y", "strand", "col"))),
                                stringsAsFactors=FALSE)
        td.bar$y <- offset+0.5
        ## Some special treatment for collapsed items. Not sure whether we want to keep this here...
        compGrps <- grep("|", names(groups), fixed=TRUE)
        if(stacking(GdObject)=="dense" && length(compGrps))
        {
            cg <- groups[compGrps]
            groups <- groups[-compGrps]
            cgn <- strsplit(names(cg), "|", fixed=TRUE)
        }
        min.swidth <- diff*.dpOrDefault(GdObject, "min.width", 2)*3
        ind <- 1
        for(ng in names(groups)){
            g <- groups[[ng]]
            r <- ranges(g)
            if(stacking(GdObject)=="dense" && length(compGrps))
            {
                intersect <- sapply(cgn, function(x, y) y %in% x, ng)
                r <-  sort(c(r, IRanges(unlist(lapply(cg[intersect], start)), unlist(lapply(cg[intersect], end)))))
            }
            color <- unique(values(g)$color)
            end <- end(r)
            start <- start(r)
            from <- end[-length(r)]
            us <- unique(as.character(strand(g)))
            if(length(us)!=1)
                us <- "*"
            fact <- if(us=="+") 1/4 else 3/4
            yspace <- if(us=="+") space else -space
            to <- from+((start[-1]-end[-length(r)])*fact)
            if(length(from))
            {
                ## Bars connecting the indiviual 'exons'
                if(!"arrow" %in% shape && !"smallArrow" %in% shape) {
                    diffs <- start[-1]-from
                    wsel <- (diffs<min.swidth)+1
                    lf <- length(from)*2
                    td.bar[ind:(ind+lf-1), "sx1"] <- c(from, to)
                    td.bar[ind:(ind+lf-1), "sx2"] <- c(to, start[-1])
                    td.bar[ind:(ind+lf-1), "col"] <- rep(color[1], lf)
                    td.bar[ind:(ind+lf-1), "strand"] <- rep(us, lf)
                    ind <- ind+lf
                } else {
                    td.bar[ind, "sx1"] <- min(start)+(min.swidth/3)
                    td.bar[ind, "sx2"] <- max(end)-(min.swidth/3)
                    td.bar[ind, "col"] <- if(length(grep("__desatCol", values(g)$feature[1])))
                        .dpOrDefault(GdObject, "fill", Gviz:::.DEFAULT_FILL_COL) else color[1]
                    td.bar[ind, "strand"] <- us
                    ind <- ind+1
                }
                
            }
            ## The optional 'exon' annotation
            if(.dpOrDefault(GdObject, "showId", FALSE))
            {
                td.bartext[ng, c("txt", "x")] <- data.frame(paste(values(g)[1, "identifier"], " "), min(start),
                                                            stringsAsFactors=FALSE)
            }
        } ## end for
    } 
    ## Finally we append all new coordinates/values to the input list
    toDraw[["box"]] <- rbind(toDraw[["box"]], boxes)
    toDraw[["bar"]] <- rbind(toDraw[["bar"]], td.bar)
    toDraw[["bartext"]] <- rbind(toDraw[["bartext"]], td.bartext)
    return(toDraw)
}

## The actual drawing method
setMethod("drawGD", signature("AnnotationTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...){
    imageMap(GdObject) <- NULL
    ## In prepare mode all we need to do is set up the stacks, so we can quit after this point.
    if(prepare){
        GdObject <- callNextMethod()
        if(subset)
            GdObject <- subset(GdObject, from=minBase, to=maxBase)
        return(invisible(GdObject))
    }
    ## In plotting mode we have to deal with stacks one by one
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    if(!length(GdObject))
        return(invisible(GdObject))
    bins <- stacks(GdObject)
    stacks <- max(bins)
    ids <- identifier(GdObject)
    hasAnno <- .dpOrDefault(GdObject, "showId", FALSE) & ids!=""
    if(any(hasAnno)) {
        txt <- paste(ids, " ")
        txt[!hasAnno] <- ""
        identifier(GdObject) <- txt
    }
    ## We first calculate all coordinates and values...
    toDraw <- list()
    pushViewport(dataViewport(xData=c(minBase, maxBase), extension=0,
                              yscale=c(1, stacks+1), clip=TRUE))
    res <- .pxResolution(coord="x")
    curVp <- vpLocation()
    if(curVp$size["width"]/stacks < .dpOrDefault(GdObject, "min.height", 3))
        stop("Too many stacks to draw. Either increase the device size or limit the drawing to a smaller region.")
    groups <- rev(split(GdObject, bins))
    ## We need to collapse the track object based on the current screen resolution...
    groups <- lapply(groups, collapseTrack, diff=res, xrange=c(minBase, maxBase))
    ## ... adjust the color saturation to indicate overplotting...
    if(.dpOrDefault(GdObject, "showOverplotting", FALSE))
    {
        dens <- as.numeric(unlist(lapply(groups, function(x) values(x)$density)))
        if(length(unique(dens))!=1)
        {
            minSat <- max(0.25, 1/max(dens))
            minDens <- min(dens)
            rDens <- diff(range(dens))
            groups <- lapply(groups, function(x){
                dens <- as.numeric(values(x)$density)
                saturation <- minSat+((dens-minDens)/rDens/(1/(1-minSat)))
                bc <- unique(.getBiotypeColor(x))
                baseCol <- rgb2hsv(col2rgb(bc))
                desatCols <- unlist(lapply(saturation, function(x) hsv(baseCol[1,], x, baseCol[3,])))
                names(desatCols) <- paste(unique(feature(x)), rep(dens, each=length(bc)), sep="_")
                feature(x) <- paste(feature(x), dens, sep="_")
                desatCols <- desatCols[unique(names(desatCols))]
                displayPars(x) <- as.list(desatCols)
                x
            })
        }
    }
    ## ... precompute all the drawing elements...
    for(i in seq_len(stacks))
        toDraw <- .prepareGene(groups[[i]], toDraw=toDraw, offset=i, minBase=minBase, maxBase=maxBase, diff=res)
    ## ... and then draw whatever is needed
    shape <- .dpOrDefault(GdObject, "shape", "arrow")
    border <- .dpOrDefault(GdObject, "col", "transparent")[1]
    lwd <- .dpOrDefault(GdObject, "lwd", 2)
    lty <- .dpOrDefault(GdObject, "lty", 1)
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
    fontsize <- .dpOrDefault(GdObject, "fontsize", 12)
    fontface <- .dpOrDefault(GdObject, "fontface", 1)
    lineheight <- .dpOrDefault(GdObject, "lineheight", 1)
    fontfamily <- .dpOrDefault(GdObject, "fontfamily", 1)
    rotation <- .dpOrDefault(GdObject, "rotation", 0)
    fontcolor <- .dpOrDefault(GdObject, "fontcolor", "white")[1]
    cex <- .dpOrDefault(GdObject, "cex", 1)
    fontcolor.group <- .dpOrDefault(GdObject, "fontcolor.group", Gviz:::.DEFAULT_SHADED_COL)[1]
    cex.group <- .dpOrDefault(GdObject, "cex", 1) * .dpOrDefault(GdObject, "cex.group", 0.6)
    fontsize.group <- .dpOrDefault(GdObject, "fontsize.group", fontsize)
    fontface.group  <- .dpOrDefault(GdObject, "fontface.group", fontface)
    fontfamily.group <- .dpOrDefault(GdObject, "fontfamily.group", fontfamily)
    box <- toDraw[["box"]]
    if(!is.null(box)){
        bar <- toDraw[["bar"]]
        if(!is.null(bar))
            .arrowBar(bar$sx1, bar$sx2, y=bar$y, bar$strand, box[,1:4, drop=FALSE], col=bar$col, lwd=lwd, lty=lty,
                      alpha=alpha, barOnly=(!"smallArrow" %in% .dpOrDefault(GdObject, "shape", "box") || stacking(GdObject)=="dense"),
                      diff=res, min.height=.dpOrDefault(GdObject, "min.height", 3))
        if("box" %in% shape)
            grid.rect(box$cx2, box$cy1, width=box$cx2-box$cx1, height=box$cy2-box$cy1,
                      gp=gpar(col=border, fill=box$fill, lwd=lwd, lty=lty, alpha=alpha),
                      default.units="native", just=c("right", "bottom"))
        if("ellipse" %in% shape){
            ellCoords <- .box2Ellipse(box)
            grid.polygon(x=ellCoords$x1, y=ellCoords$y1, id=ellCoords$id,
                         gp=gpar(col=border, fill=box$fill, lwd=lwd, lty=lty, alpha=alpha),
                         default.units="native")
        }
        if(any(c("arrow", "smallArrow") %in% shape)){
            .filledArrow(box[,1:4], col=border, fill=box$fill, lwd=lwd, lty=lty, alpha=alpha, strand=box$strand, min.width=6*res)
        }
        if(.dpOrDefault(GdObject, "showFeatureId", FALSE))
            grid.text(box$text, box$textX, box$textY, rot=rotation,
                      gp=gpar(col=fontcolor, cex=cex, fontsize=fontsize, fontface=fontface, lineheight=lineheight,
                              fontfamily=fontfamily), default.units="native", just=c("center", "center"))
        bartext <- toDraw[["bartext"]]
        if(!is.null(bartext) && stacking(GdObject)!="dense")
            grid.text(bartext$txt, bartext$x, bartext$y, gp=gpar(col=fontcolor.group, cex=cex.group, fontsize=fontsize.group,
                                                                 fontface=fontface.group, fontfamily=fontfamily.group),
                      default.units="native", just=c("right", "center"))
    }
    popViewport(1)
    ## Set up the image map
    im <- if(!is.null(toDraw$box)) {
        coords <- as.matrix(toDraw$box[,c("x1", "y1", "x2", "y2"),drop=FALSE])
        restCols <- setdiff(colnames(toDraw$box), c("x1", "x2", "y1", "y2", "cx1", "cx2", "cy1", "cy2", "textX", "textY"))
        tags <- sapply(restCols, function(x){
            tmp <- as.character(toDraw$box[,x])
            names(tmp) <- rownames(coords)
            tmp}, simplify=FALSE)
        tags$title <- gsub(" +", "", identifier(GdObject))
        ImageMap(coords=coords, tags=tags) } else NULL
    imageMap(GdObject) <- im
    return(invisible(GdObject))
})

## For a GeneRegionTrack we just set the showExonId alias
setMethod("drawGD", signature("GeneRegionTrack"), function(GdObject,  ...){
    displayPars(GdObject) <- list(showFeatureId=as.vector(displayPars(GdObject, "showExonId")))
    callNextMethod(GdObject, ...)
})
##----------------------------------------------------------------------------------------------------------------------------






##----------------------------------------------------------------------------------------------------------------------------
## Draw a genome axis
##----------------------------------------------------------------------------------------------------------------------------
setMethod("drawGD", signature("GenomeAxisTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
    ## Nothing to do if the coordinate width is 0, so we can quit right away
    imageMap(GdObject) <- NULL
    if((maxBase-minBase)==0)
        return(invisible(GdObject))
    ## We start by setting up the drawing canvas
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    pushViewport(dataViewport(xData=c(minBase, maxBase), yscale=c(-1, 1), extension=0))
    ## Create a useful data range for the axis
    pres <- .pxResolution()
    curVp <- vpLocation()
    cex <- .dpOrDefault(GdObject, "cex", 0.8)
    lwd <- .dpOrDefault(GdObject, "lwd", 1)
    fontface <- .dpOrDefault(GdObject, "fontface", 1)
    add53 <- .dpOrDefault(GdObject, "add53", FALSE)
    add35 <- .dpOrDefault(GdObject, "add35", FALSE)
    lcex <- cex*0.75
    textYOff <-  pres["y"]*3
    textXOff <-  pres["x"]*2
    endMargin <- if(add53 || add35)
        (as.numeric(convertWidth(stringWidth("5'"),"native"))*lcex)+(textXOff*2) else pres["x"]*5
    axRange <- c(minBase+endMargin, maxBase-endMargin)
    ## We want fixed vertical sizes for axis tracks to avoid akward stretching effects.
    color <- .dpOrDefault(GdObject, "col", "darkgray")[1]
    littleTicks <- .dpOrDefault(GdObject, "littleTicks", FALSE)
    dfact <- max(1, .dpOrDefault(GdObject, "distFromAxis", 1))
    labelPos <- .dpOrDefault(GdObject, "labelPos", "alternating")
    lwdAdd <- (lwd-1)/2
    tickHeight <- (ifelse(littleTicks, 2, 1) * 3 * dfact + lwdAdd) * pres["y"]
    ids <- values(GdObject)$id
    showIds <- .dpOrDefault(GdObject, "showId", FALSE) && !is.null(ids) && !all(ids=="")
    rcex <- .dpOrDefault(GdObject, "cex.id", 0.7)
    rcol <- .dpOrDefault(GdObject, "col.id", "white")
    sep <- (if(length(GdObject)){
        if(showIds)
            max(1.5, ((max(as.numeric(convertHeight(stringHeight(ids),"native"))*rcex)+textYOff)/pres["y"])/2) else
        1.5} else 1)+lwdAdd
    pyOff <- pres["y"]*sep
    ## In prepare mode we just want to figure out the optimal size
    if(prepare)
    {
        nsp <- (sum(tickHeight, pyOff*2, textYOff*2 +
                    (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex)*2*1.3)/pres["y"]
        displayPars(GdObject) <- list("neededVerticalSpace"=nsp)
        popViewport(1)
        return(invisible(GdObject))
    }
    ## Plot range if there is any
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
	
    ## in "scale" mode we just plot a simple scale and return ...
    scaleLen <- .dpOrDefault(GdObject, "scale", NULL)
    if(!is.null(scaleLen))
    {
        ##browser()
        len <- (maxBase-minBase + 1)
        if (scaleLen > len) {
            warning(paste("scale (", scaleLen,
                          ") cannot be larger than plotted region",
                          len, " - setting to ~5%\n", sep=""))
            scaleLen = 0.05
        }
        xoff <- len * 0.03 + minBase
        labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
        v <- scaleLen
        ## TODO: for consistency acknowledge the 'exponent' argument if not NULL
        ## currently the exponent for the scale is calculated automatically
        ex <- floor(log10(scaleLen))
        u <- paste(scaleLen, "b")
        if(scaleLen <= 1 && scaleLen > 0) {
            scaleLen <- len * scaleLen 
            ex <- floor(log10(scaleLen))
            u <- paste(round(scaleLen, -ex), "b")
            ##v <- round(scaleLen, -(ex+round(scaleLen/(10^(ex+1))))) #Êrounds to nearest exponent
            v <- round(scaleLen, -ex)
        }
        if(ex >=9) {
            u <- paste(v/1e9, "gb")
        } else if(ex >= 6) {
            u <- paste(v/1e6, "mb")
        } else if(ex >= 3) {
            u <- paste(v/1e3, "kb")
        }
        grid.lines(x=c(xoff, v+xoff), y=c(0,0), default.units="native", gp=gpar(col=color, lwd=lwd, alpha=alpha))
        grid.segments(x0=c(xoff, v+xoff), y0=c(0-tickHeight, 0-tickHeight),
                      x1=c(xoff, v+xoff), y1=c(tickHeight,tickHeight, tickHeight),
                      default.units="native", gp=gpar(col=color, lwd=lwd, alpha=alpha))
        z <- len * 0.01
        if(labelPos=="below"){
            grid.text(label=u, x=xoff+v/2, y=0-(tickHeight/1.5*dfact), just=c("center", "top"),
                      gp=gpar(alpha=alpha, col=color, cex=cex, fontface=fontface), default.units="native")
        } else if(labelPos=="above"){
            grid.text(label=u, x=xoff+v/2, y=tickHeight/1.5*dfact, just=c("center", "bottom"),
                      gp=gpar(alpha=alpha, col=color, cex=cex, fontface=fontface), default.units="native")
        } else {
            grid.text(label=u, x=v+xoff+z, y=0, just=c("left", "center"),
                      gp=gpar(alpha=alpha, col=color, cex=cex, fontface=fontface), default.units="native")
        }
        popViewport(1)
        return(invisible(GdObject))
    }
    
    if(length(GdObject))
    {
        rfill <- .dpOrDefault(GdObject, "fill.range", "cornsilk3")
        rcolor <- .dpOrDefault(GdObject, "col.range", "cornsilk4")[1]
        diff <- .pxResolution(coord="x")
        GdObject <- collapseTrack(GdObject, diff=diff, xrange=c(minBase, maxBase))
        start(GdObject) <- pmax(axRange[1], start(GdObject))
        end(GdObject) <- pmin(axRange[2], end(GdObject))
        coords <- cbind(start(GdObject), -0.1, end(GdObject), 0.1)
        grid.rect(x=start(GdObject), y=-pyOff, width=width(GdObject), height=pyOff*2,
                  default.units="native", just=c("left", "bottom"), gp=gpar(col=rcolor, fill=rfill, alpha=alpha))
        vals <- values(GdObject)
        if(showIds)
            grid.text(ids, x=start(GdObject) + width(GdObject)/2, y=0,
                      gp=gpar(col=rcol, cex=rcex, fontface=fontface),
                      default.units="native", just=c("center", "center"))
        ## Calculate the coordinates for the image map
        map <- as.matrix(.getImageMap(coords))
        rownames(map) <- paste("region", seq_len(nrow(map)), sep="_")
        tags <- lapply(list(title=rownames(map), start=as.character(start(GdObject)), end=as.character(end(GdObject))),
                       function(x){ names(x) <- rownames(map); x})
        imageMap(GdObject) <- ImageMap(coords=map, tags=tags) 
    }
    ## width<1, we can return here, no need for tick marks
    if(abs(diff(axRange))<1){
        popViewport()
        return(invisible(GdObject))
    }
    ## We want two parallel lines with little hooks on the ends
    pyHook <- pres["y"]*(sep+2+lwdAdd)
    pxOff <- pres["x"]*5
    grid.segments(x0=rep(axRange[1], 2), y0=c(-1,1)*pyOff, x1=rep(axRange[2], 2),  y1=c(-1,1)*pyOff,
                  default.units="native", gp=gpar(col=color, alpha=alpha, lwd=lwd))
    grid.segments(x0=c(axRange[2]-pxOff, axRange[1]),  y0=c(pyHook, -pyOff),
                  x1=c(axRange[2], axRange[1]+pxOff), y1=c(pyOff, -pyHook),
                  default.units="native", gp=gpar(col=color, alpha=alpha, lwd=lwd))
    ## Here we plot the top level ticks
    tck <- .ticks(axRange)
    tck <- tck[tck<axRange[2]-pxOff*2 & tck>axRange[1]+pxOff*2]
    y0t <- rep(c(1,-1)*pyOff, length(tck))[1:length(tck)]
    y1t <- y0t + rep(c(tickHeight, -tickHeight), length(tck))[1:length(tck)]
    labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
    y0t <- switch(labelPos, "alternating"=y0t, "revAlternating"=-y0t, "above"=abs(y0t), "below"=-abs(y0t), "beside"=y0t)
    y1t <- switch(labelPos, "alternating"=y1t, "revAlternating"=-y1t, "above"=abs(y1t), "below"=-abs(y1t), "beside"=y1t)
    grid.segments(x0=tck, x1=tck, y0=y0t, y1=y1t,  default.units="native", gp=gpar(col=color, alpha=alpha, lwd=lwd, lineend="square"))
    ## The top level tick labels
    tckText <- tck
    formatValue <- "d"
    exponent <- max(0, .dpOrDefault(GdObject, "exponent", NULL))
    exponent <- if(is.null(.dpOrDefault(GdObject, "exponent", NULL)))
    {
        exp <- 0
        while(all(tck[tck>0]/10^exp > 1))
            exp <- exp+3
        exp-3
    } else  max(0, getPar(GdObject, "exponent"))
    if(exponent > 0){
        tckText <- tckText/(10^exponent)
        if(!exponent %in% c(3, 6, 9))
            formatValue <- "g"
    }
    label <- switch(as.character(exponent),
                    "0"=sprintf("%i", as.integer(tckText)),
                    "3"=sprintf("%s kb", tckText),
                    "6"=sprintf("%s mb", tckText),
                    "9"=sprintf("%s gb", tckText),
                    sapply(tckText, function(x) bquote(paste(.(x), " ",10^.(exponent)))))
    ylabs <- y1t + (ifelse(y1t>0, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex))
    if(is.character(label))
        grid.text(label=label, x=tck, y=ylabs, just=c("centre", "centre"),
                  gp=gpar(cex=cex, fontface=fontface), default.units="native")
    else
        for(i in seq_along(label))
             grid.text(label=label[[i]], x=tck[i], y=ylabs[i], just=c("centre", "centre"),
                       gp=gpar(cex=cex, fontface=fontface), default.units="native")
    ## The scecond level ticks and labels if necessary
    if (.dpOrDefault(GdObject, "littleTicks", FALSE) && length(tck)>1)
    {
        avSpace <- min(diff(tck))
        spaceFac <- 1.8
        spaceNeeded <- min(as.numeric(convertWidth(stringWidth(if(is.character(label)) label else "000000000"),"native"))/2)*lcex*spaceFac
        nTcks <- (avSpace %/% spaceNeeded)
        if(nTcks%%2 == 0)
            nTcks <- nTcks-1
        btck <- tck
        if (!(minBase %in% btck))
            btck <- c(minBase, btck)
        if (!(maxBase %in% btck))
            btck <- c(btck, maxBase)
        y0lt <- y1lt <- ltck <- NULL
        for(i in seq_len(length(btck)-1))
        {
            toFill <- btck[i:(i+1)]
            ttck <- if(i==1) rev(toFill[2]-(avSpace/nTcks)*seq_len(nTcks-1)) else toFill[1]+(avSpace/nTcks)*seq_len(nTcks-1)
            ltck <- c(ltck, ttck)
            ord <- if(i==1){ if(y0t[1]>0) c(1,-1) else c(-1,1) } else if(y0t[i-1]<0) c(1,-1) else c(-1,1)
            y0 <- rep(ord*pyOff, length(ttck))[1:length(ttck)]
            y1 <- y0 + rep(ord*tickHeight/2, length(ttck))[1:length(ttck)]
            y0lt <- c(y0lt, switch(labelPos, "alternating"=y0, "revAlternating"=y0, "above"=abs(y0), "below"=-abs(y0)))
            y1lt <- c(y1lt, switch(labelPos, "alternating"=y1, "revAlternating"=y1, "above"=abs(y1), "below"=-abs(y1)))
        }
        endPadding <- pres["x"]*15
        sel <- ltck > min(tck, axRange+endPadding) & ltck < max(tck, axRange-endPadding)
        if(length(ltck[sel]) && min(diff(tck))>nTcks)
        {
            grid.segments(x0=ltck[sel], x1=ltck[sel], y0=y0lt[sel], y1=y1lt[sel],  default.units="native",
                          gp=gpar(col=color, alpha=alpha, lwd=lwd, lineend="square"))
            ltckText <- ltck[sel]
            if(exponent > 0)
                ltckText <- ltckText/(10^exponent)
            llabel <- switch(as.character(exponent),
                             "0"=sprintf("%i", as.integer(ltckText)),
                             sprintf("%g", ltckText))
            ytlabs <- y1lt + (ifelse(y1lt>0, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*lcex))
            if(is.character(label))
                grid.text(label=llabel, x=ltck[sel], y=ytlabs[sel], just=c("centre", "centre"),
                          gp=gpar(cex=lcex, fontface=fontface), default.units="native")
            else
                for(i in seq_along(llabel))
                    grid.text(label=llabel[[i]], x=ltck[sel][i], y=ytlabs[sel][i], just=c("centre", "centre"),
                              gp=gpar(cex=lcex, fontface=fontface), default.units="native")  
        }
    }
    ## The direction indicators
    if(add53)
    {
        grid.text(label=expression("5'"), x=axRange[1]-textXOff, y=pyOff,
                  just=c("right", "bottom"), gp=gpar(cex=cex*.75, fontface=fontface),
                  default.units="native")
        grid.text(label=expression("3'"), x=axRange[2]+textXOff, y=pyOff,
                  just=c("left", "bottom"), gp=gpar(cex=cex*.75, fontface=fontface),
                  default.units="native")
    }
    if(add35)
    {
        grid.text(label=expression("3'"), x=axRange[1]-textXOff, y=-pyOff,
                  just=c("right", "top"), gp=gpar(cex=cex*.75, fontface=fontface),
                  default.units="native")
        grid.text(label=expression("5'"), x=axRange[2]+textXOff, y=-pyOff,
                  just=c("left", "top"), gp=gpar(cex=lcex, fontface=fontface),
                  default.units="native")
    }
    popViewport()
    return(invisible(GdObject))})
##----------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------
## Draw DetailsAnnotationTrack
##----------------------------------------------------------------------------------------------------------------------------
setMethod("drawGD", signature("DetailsAnnotationTrack"),
		function(GdObject,  minBase, maxBase, prepare=FALSE, ...){	
			
			#browser()   
			if ( prepare ) {
				GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
				return(invisible(GdObject))
			}
			
			col <- .dpOrDefault(GdObject, "detailsConnector.col", fromPrototype=T)
			lty <- .dpOrDefault(GdObject, "detailsConnector.lty", fromPrototype=T)
			lwd <- .dpOrDefault(GdObject, "detailsConnector.lwd", fromPrototype=T)
			pch <- .dpOrDefault(GdObject, "detailsConnector.pch", fromPrototype=T)
			cex <- .dpOrDefault(GdObject, "detailsConnector.cex", fromPrototype=T)
			border.lty <- .dpOrDefault(GdObject, "detailsBorder.lty", fromPrototype=T)
			border.lwd <- .dpOrDefault(GdObject, "detailsBorder.lwd", fromPrototype=T)
			border.col <- .dpOrDefault(GdObject, "detailsBorder.col", fromPrototype=T)
			border.fill <- .dpOrDefault(GdObject, "detailsBorder.fill", fromPrototype=T)
			minwidth <- .dpOrDefault(GdObject, "details.minWidth", fromPrototype=T)
			size <- .dpOrDefault(GdObject, "details.size", fromPrototype=T)
			xyratio <- .dpOrDefault(GdObject, "details.ratio", fromPrototype=T)
			if ( 0 >= size || size > 1 ) {
				warning("details.size must be >0 and <1 - reset to 0.5")
				size = 0.5
			}
			
			len = length(GdObject)
			if( ((maxBase-minBase)/len)/.pxResolution(coord="x") < minwidth ) {
				warning("too much detail for availalbe space (plot fewer annotation or increase details.minWidth)!")
				popViewport(1)
				GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
				return(GdObject)
			}
			#GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
			GdObject <- GdObject[order(start(GdObject))]
			bins <- stacks(GdObject)
			stacks <- max(bins)
			xloc1 <- (end(GdObject) - start(GdObject))/2+start(GdObject)
			#ord <- order(xloc1, decreasing=FALSE)
			yloc1 <- (stacks - (bins - 0.5)+1)
			xloc2 <- ((1/len*seq_len(len))-1/len + (1/len*0.5))
			yloc2 <- rep(1, len) 
			
			## draw details plots (via user supplied function 'fun')
			pushViewport(viewport(height=size, y=1-size, just=c(0.5, 0)))
			w = 1
			v = 0
			vpl = vpLocation()
			r =  vpl$size["width"]/len/vpl$size["height"]
			if ( r > xyratio ) {
				w = xyratio/r
				v = ((1/len) - (1/len*w))/2
			}
			args <- .dpOrDefault(GdObject, "detailsFunArgs", fromPrototype=T)
			for ( i in 1:length(GdObject) ) {
				#pushViewport(viewport(width=1/len, x=(1/len*i)-1/len, just=c(0, 0.5)))
				pushViewport(viewport(width=1/len*w, x=((1/len*i)-1/len)+(v), just=c(0, 0.5)))
				grid.rect(gp=gpar(col=border.col, lwd=border.lwd, lty=border.lty, fill=border.fill))
				args$start = start(GdObject[i])
				args$end = end(GdObject[i])
				args$strand = strand(GdObject[i])
				args$chromosome = chromosome(GdObject[i])
				args$identifier = identifier(GdObject[i], lowest=TRUE)
				args$index = i
				do.call(GdObject@fun, args)
				popViewport(1)
			}
			popViewport(1)
			
			## plot AnnotationTrack and connectors to details
			pushViewport(viewport(xscale=c(minBase, maxBase),
							yscale=c(1, stacks+1), clip=FALSE,
							height=1-size, y=0, just=c(.5, 0)))	
			GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
			grid.segments(x0=unit(xloc1, "native"), x1=xloc2, y0=unit(yloc1, "native"),
					y1=yloc2, gp=gpar(col=col, lwd=lwd, lty=lty, cex=cex))
			grid.points(x=unit(xloc2, "npc"), y=unit(yloc2, "npc"), gp=gpar(col=col, cex=cex), pch=pch)
			grid.points(x=unit(xloc1, "native"), y=unit(yloc1, "native"), gp=gpar(col=col, cex=cex), pch=pch)
			popViewport(1)
			
			return(invisible(GdObject))
		})

##----------------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------------
## Draw a data track
##----------------------------------------------------------------------------------------------------------------------------
## Helper function to return the absolute extreme value in a vector
.extreme <- function(x) if(all(is.na(x))) NA else x[which.max(abs(x))]

## Map numeric range to values from 1 to n
.z2icol <- function(z, n, xrange=range(z, na.rm=TRUE))
{
    res <- round((z - xrange[1])/diff(xrange) * (n - 1)) + 1
    res[res > n] = n
    res[res < 1] = 1
    return(res)
}

setMethod("drawGD", signature("DataTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
    imageMap(GdObject) <- NULL
    type <- .dpOrDefault(GdObject, "type", "p")
    type <- match.arg(type, c("p", "l", "b", "a", "s", "g", "r", "S", "smooth",
                              "histogram", "mountain", "h", "boxplot", "gradient", "heatmap"),
                      several.ok=TRUE)
    ## Grouping may be useful for some of the plot types, may be ignored for others
    vals <- values(GdObject)
    groups <- .dpOrDefault(GdObject, "groups")
    if(!is.null(groups) && length(groups) != nrow(vals))
        stop("'groups' must be a vector of similar length as the number of rows in the data matrix (", nrow(vals), ")")
    if(!is.null(groups) && !is.factor(groups))
        groups <- factor(groups)
    stacked  <- .dpOrDefault(GdObject, "stackedBars", FALSE)
    ## The general "col" parameter should be the default for all relevant colors except when there are groups.
    pcols <- .getPlottingFeatures(GdObject)
    ## In prepare mode we collapse the track to allow for aggregation and so on since we need the final data
    ## values to draw the axis. 
    if(prepare)
    {
        if(subset)
            GdObject <- subset(GdObject, from=minBase, to=maxBase)
        pushViewport(viewport(xscale=c(minBase, maxBase), yscale=c(0,1), clip=TRUE))
        diff <- .pxResolution(coord="x")
        GdObject <- collapseTrack(GdObject, diff=diff, xrange=c(minBase, maxBase))
        popViewport(1)
        ## If we have groups and stacked histograms we have to adjust the ylim values, also for regular histograms
        if("histogram" %in% type)
        {
            vals <- values(GdObject)
            groups <- rep(groups, ncol(vals))
            ylim <- displayPars(GdObject, "ylim")
            agFun <- .aggregator(GdObject)
            if(!is.null(groups) && nlevels(groups)>1)
            {
                valsS <- matrix(sapply(split(vals, groups), function(x) agFun(t(matrix(x, ncol=ncol(vals))))), nrow=ncol(vals))
                displayPars(GdObject) <- list(".__valsS"=valsS)
                if(stacked==TRUE && is.null(ylim))
                {
                    ylim <- range(unlist(apply(valsS, 1, function(x){
                        x <- x[!is.na(x)]
                        sel <- x>=0
                        tmp <- NULL
                        if(!all(is.na(sel)))
                        {
                            if(any(sel))
                                tmp <- c(min(x[sel]), sum(x[sel]))
                            if(any(!sel))
                                tmp <- c(max(x[!sel]), tmp, sum(x[!sel]))
                        }
                        tmp})))
                    if(length(type)>1)
                        ylim <- range(c(ylim, vals))
                    displayPars(GdObject) <- list(ylim=ylim)
                }
            } else {
                if(is.null(ylim))
                {
                    valsA <- agFun(t(vals))
                    ylim <- if(!length(valsA)) c(-1,1) else c(min(c(0,valsA), na.rm=TRUE), max(valsA, na.rm=TRUE))
                    if(length(type)>1)
                        ylim <- range(c(ylim, vals), na.rm=TRUE)
                    displayPars(GdObject) <- list(ylim=ylim)
                }
            }
        }
        ## If we want a legend we have to figure out how much vertical space is needed
        grps <- .dpOrDefault(GdObject, "groups")
        if(is.null(grps) || length(setdiff(type, c("gradient", "mountain", "grid"))) == 0)
            displayPars(GdObject) <- list(legend=FALSE)
        if(.dpOrDefault(GdObject, "legend", FALSE)){
            cex <- .dpOrDefault(GdObject, "cex.legend", 0.8)
            fontsize <- .dpOrDefault(GdObject, "fontsize.legend", 12)
            fontface <- .dpOrDefault(GdObject, "fontface.legend", 1)
            lineheight <- .dpOrDefault(GdObject, "lineheight.legend", 1)
            fontfamily <- .dpOrDefault(GdObject, "fontfamily.legend", 1)
            pushViewport(viewport(width=unit(1, "npc")-unit(0.2,"inches"),
                                  gp=gpar(cex=cex, fontsize=fontsize, fontface=fontface,
                                          lineheight=lineheight)))
            grps <- levels(factor(grps))
            legInfo <- .legendInfo()[type,, drop=FALSE]
            for(i in colnames(legInfo))
                legInfo[,i] <- any(legInfo[,i]) && !any(duplicated(pcols[[i]][1:length(grps)]))
            legFactors <- sort(names(which(apply(legInfo, 2, any))))
            boxSize <-  if(length(setdiff(legFactors, c("col", "cex")))==0) 0.1 else 0.3
            spacing <- 0.1
            hspacing <- 0.02
            lengths <- as.numeric(convertWidth(stringWidth(grps),"inches"))
            heights <- as.numeric(convertWidth(stringHeight(grps),"inches"))
            colWidth <- max(lengths + boxSize + spacing*2)
            availSpace <- vpLocation()$isize
            colNum <- max(1, availSpace["width"] %/% colWidth)
            rowNum <- ceiling(length(grps)/colNum)
            rowHeight <- max(c(heights, 0.1))
            vertSpace <- (rowHeight * rowNum) + (hspacing * (rowNum-1)) + 0.2
            displayPars(GdObject) <- list(".__verticalSpace"=vertSpace, ".__layoutDims"=c(rowNum, colNum),
                                          ".__boxSize"=boxSize, ".__spacing"=spacing, ".__groupLevels"=grps,
                                          ".__legFactors"=legFactors)
            popViewport(1)
        }
        #popViewport(1)
        return(invisible(GdObject))
    }
    ## We only proceed if there is something to draw within the ranges, but still may have to add the grid
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
    if(!length(GdObject))
    {
        if ("g" %in% type)
            panel.grid(h=.dpOrDefault(GdObject, "h", -1), v=.dpOrDefault(GdObject, "v", -1),
                       col=.dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty=.dpOrDefault(GdObject, "lty.grid", 1),
                       lwd=.dpOrDefault(GdObject, "lwd.grid", 1), alpha=alpha)
        return(invisible(GdObject))
    }
    vals <- values(GdObject)
    ylim <- .dpOrDefault(GdObject, "ylim", range(vals, na.rm=TRUE, finite=TRUE))
    if(diff(ylim)==0)
        ylim <- ylim+c(-1,1)
    ylimExt <- extendrange(r=ylim, f=0.05)
    ## The optional legend is plotted underneath the data
    grpLevels <- .dpOrDefault(GdObject, ".__groupLevels")
    if(.dpOrDefault(GdObject, "legend", FALSE)){
        lSpace <- getPar(GdObject, ".__verticalSpace")
        pushViewport(viewport(y=1, height=unit(1, "npc") - unit(lSpace, "inches"),
                              just=c(0.5, 1)))
        on.exit({popViewport(1)
                 cex <- .dpOrDefault(GdObject, "cex.legend", 0.8)
                 legFactors <- .dpOrDefault(GdObject, ".__legFactors", character())
                 fontsize <- .dpOrDefault(GdObject, "fontsize.legend", 12)
                 fontface <- .dpOrDefault(GdObject, "fontface.legend", 1)
                 lineheight <- .dpOrDefault(GdObject, "lineheight.legend", 1)
                 fontfamily <- .dpOrDefault(GdObject, "fontfamily.legend", 1)
                 fontcolor <- .dpOrDefault(GdObject, "fontcolor.legend", Gviz:::.DEFAULT_SHADED_COL)
                 pushViewport(viewport(y=0, height=unit(lSpace, "inches"), just=c(0.5, 0),
                                       gp=gpar(cex=cex, fontsize=fontsize, fontface=fontface, fontcolor=fontcolor,
                                               lineheight=lineheight)))
                 pushViewport(viewport(width=unit(1, "npc") - unit(0.1, "inches"), height=unit(1, "npc") - unit(0.1, "inches")))
                 boxSize <- getPar(GdObject, ".__boxSize")
                 spacing <- getPar(GdObject, ".__spacing")
                 dims <- getPar(GdObject, ".__layoutDims")
                 for(i in seq_along(grpLevels)){
                     row <- (((i)-1) %/% dims[2])+1
                     col <- (((i)-1) %% dims[2])+1
                     pushViewport(viewport(width=1/dims[2], height=1/dims[1], x=(1/dims[2])*(col-1), y=1-((1/dims[1])*(row-1)), just=c(0,1)))
                     if(length(setdiff(legFactors, c("col")))==0){
                         grid.rect(width=unit(boxSize, "inches"), height=unit(boxSize, "inches"), x=0, just=c(0, 0.5),
                                   gp=gpar(fill=pcols$col[i], col=Gviz:::.DEFAULT_SHADED_COL))
                     } else {
                         if(any(c("pch", "col.symbol") %in% legFactors))
                             panel.points(unit(boxSize/2, "inches"), 0.5, pch=pcols$pch[i], cex=pcols$cex[i], col=pcols$col.symbol[i])
                         if(any(c("lwd", "lty", "col.lines") %in% legFactors))
                             panel.lines(unit(c(0,boxSize), "inches"), c(0.5, 0.5), col=pcols$col.line[i], lwd=pcols$lwd[i], lty=pcols$lty[i])
                     }
                     grid.text(x=unit(boxSize+spacing, "inches"), just=c(0, 0.5), label=grpLevels[i], gp=gpar(col=fontcolor))
                     popViewport(1)
                 }
                 popViewport(2)
             })
    } else {
        #on.exit(popViewport(1))
    }
    pushViewport(viewport(xscale=c(minBase, maxBase), yscale=ylimExt, clip=TRUE))
    ## The plotting parameters, some defaults from the lattice package first
    plot.symbol <- trellis.par.get("plot.symbol")
    superpose.symbol <- trellis.par.get("superpose.symbol")
    superpose.line <- trellis.par.get("superpose.line")
    groups <- rep(groups, ncol(vals))
    ## For loess calculation we need some settings
    span <- .dpOrDefault(GdObject, "span", 1/5)
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
    x <- rep(pos, each=nrow(vals))
    y <- as.numeric(vals)
    ## A grid should always be plotted first, so we need to catch this here
    wg <- match("g", type, nomatch = NA_character_)
    if (!is.na(wg)) {
        panel.grid(h=.dpOrDefault(GdObject, "h", -1), v=.dpOrDefault(GdObject, "v", -1),
                   col=pcols$col.grid, lty=pcols$lty.grid, lwd=pcols$lwd.grid)
        type <- type[-wg]
    }
    ## The special type 'mountain' has to be handled separately
    if("mountain" %in% type) {
        mbaseline <- if(is.null(baseline)) 0 else baseline[1]
        fill.mountain <- .dpOrDefault(GdObject, "fill.mountain", superpose.symbol$fill)[1:2]
        col.mountain <- .dpOrDefault(GdObject, "col.mountain", pcols$col)[1]
        col.baseline <- .dpOrDefault(GdObject, "col.baseline", col.mountain)[1]
        lwd.mountain <- .dpOrDefault(GdObject, "lwd.mountain", pcols$lwd)[1]
        lty.mountain <- .dpOrDefault(GdObject, "lty.mountain", pcols$lty)[1]
        .panel.mountain(x, y, col=col.mountain, fill=fill.mountain, span=span, degree=degree, family=family,
                        evaluation=evaluation, lwd=lwd.mountain, lty=lty.mountain, col.line=col.mountain, alpha=alpha,
                        baseline=mbaseline)
        if(!is.na(mbaseline))
            panel.abline(h=mbaseline, col=col.baseline, lwd=lwd.baseline, lty=lty.baseline, alpha=alpha)
    }
    ## Also the type 'boxplot' is handled up front
    if("boxplot" %in% type)
    {
        box.ratio <- .dpOrDefault(GdObject, "box.ratio", 1)
        box.width <- .dpOrDefault(GdObject, "box.width", (min(diff(unique(sort(x))))*0.5)/box.ratio)
        diff <- .pxResolution(coord="x")
        if(!is.null(groups))
        {
            tw <- min(width(GdObject))
            spacer <- diff
            nb <- nlevels(groups)
            bw <- .dpOrDefault(GdObject, "box.width", (tw-(nb+2)*spacer)/nb)
            bcex <- min(pcols$cex[1], (bw/diff)/20)
            by <- lapply(split(vals, groups), matrix, ncol=ncol(vals))
            for(j in seq_along(by))
            {
                xx <- rep(start(GdObject)+(j*spacer)+(j*bw), each=nrow(by[[j]]))-(bw/2)
                .panel.bwplot(xx, as.numeric(by[[j]]), box.ratio=box.ratio, box.width=(bw/2)/box.ratio, pch=pcols$pch[1],
                              lwd=pcols$lwd[1], lty=pcols$lty[1], fontsize=fontsize,
                              col=pcols$col.histogram, cex=bcex, font=font, fontfamily=font, fontface=fontface,
                              fill=pcols$col[j], varwidth=.dpOrDefault(GdObject, "varwidth", FALSE), 
                              notch=.dpOrDefault(GdObject, "notch", FALSE), notch.frac=.dpOrDefault(GdObject, "notch.frac", 0.5),
                              levels.fos=.dpOrDefault(GdObject, "level.fos", sort(unique(xx))), 
                              stats=.dpOrDefault(GdObject, "stats", boxplot.stats), coef=.dpOrDefault(GdObject, "coef", 1.5), 
                              do.out=.dpOrDefault(GdObject, "do.out", TRUE), alpha=alpha)
            }
            diffY <- .pxResolution(coord="y", 2)
            outline <- apply(vals, 2, range)
            grid.rect(start(GdObject), outline[1,]-diffY, width=width(GdObject), height=abs(outline[2,]-outline[1,])+(2*diffY),
                  gp=gpar(col=pcols$col.histogram, fill="transparent", alpha=alpha, lty="dotted"),
                  default.units="native", just=c("left", "bottom"))
        } else {
            bcex <- min(pcols$cex[1], ((box.width*2)/diff)/20)
            .panel.bwplot(x, y, box.ratio=box.ratio, box.width=box.width, pch=pcols$pch[1],
                          lwd=pcols$lwd[1], lty=pcols$lty[1], fontsize=fontsize,
                          col=pcols$col.histogram, cex=bcex, font=font, fontfamily=font, fontface=fontface,
                          fill=pcols$fill[1], varwidth=.dpOrDefault(GdObject, "varwidth", FALSE), 
                          notch=.dpOrDefault(GdObject, "notch", FALSE), notch.frac=.dpOrDefault(GdObject, "notch.frac", 0.5),
                          levels.fos=.dpOrDefault(GdObject, "level.fos", sort(unique(x))), 
                          stats=.dpOrDefault(GdObject, "stats", boxplot.stats), coef=.dpOrDefault(GdObject, "coef", 1.5), 
                          do.out=.dpOrDefault(GdObject, "do.out", TRUE), alpha=alpha)
        }
    }
    ## 'histogram' fills up the full range area if its width is > 1
    if("histogram" %in% type)
    {
        ylimSort <- sort(ylimExt)
        yy <- if(ylimSort[1]<=0 && ylimSort[2]>=0) 0 else ylimSort[1]
        if(!is.null(groups) && nlevels(groups)>1)
        {
            valsS <- displayPars(GdObject, ".__valsS")
            if(stacked)
            {
                curMinPos <- curMaxPos <- rep(yy, nrow(valsS))
                for(s in seq_len(ncol(valsS)))
                {
                    if(!all(is.na(valsS[,s])))
                    {
                        sel <- !is.na(valsS[,s]) & valsS[,s]>=0
                        yyy <- curMinPos
                        yyy[sel] <- curMaxPos[sel]
                        offset <- yyy
                        offset[offset!=yy] <- 0
                        grid.rect(start(GdObject), yyy, width=width(GdObject), height=valsS[,s]-offset,
                                  gp=gpar(col="transparent", fill=pcols$col[s], lwd=pcols$lwd[1], lty=pcols$lty[1], alpha=alpha), default.units="native",
                                  just=c("left", "bottom"))
                        curMaxPos[sel] <- curMaxPos[sel]+(valsS[sel,s]-offset[sel])
                        curMinPos[!sel] <- curMinPos[!sel]+(valsS[!sel,s]-offset[!sel])
                    }
                }
                diff <- .pxResolution(coord="x", pcols$lwd[1]+1)
                tooNarrow <- width(GdObject)<diff
                if(!all(tooNarrow))
                    grid.rect(start(GdObject)[!tooNarrow], curMinPos[!tooNarrow], width=width(GdObject)[!tooNarrow],
                              height=(curMaxPos-curMinPos)[!tooNarrow],
                              gp=gpar(fill="transparent", col=pcols$col.histogram, lwd=pcols$lwd[1], lty=pcols$lty[1], alpha=alpha),
                              default.units="native", just=c("left", "bottom"))
            } else {
                spacer <- .pxResolution(min.width=1, coord="x")
                yOff <- .pxResolution(min.width=1, coord="y")
                outline <- apply(valsS, 1, function(x) range(c(yy, x), na.rm=TRUE))
                grid.rect(start(GdObject), outline[1,]-yOff, width=width(GdObject), height=apply(outline, 2, diff)+(yOff*2),
                          gp=gpar(col=pcols$col.histogram, fill=pcols$fill.histogram, lwd=pcols$lwd[1], lty=pcols$lty[1], alpha=alpha), default.units="native",
                          just=c("left", "bottom"))
                len <- ncol(valsS)
                subW <- (width(GdObject)-(spacer*(len+1)))/len
                sel <- subW > spacer
                ## FIXME: how do we treat this if there is not enough space to plot?
                sel <- !logical(length(subW))
                if(any(sel))
                {
                    subW <- subW[sel]
                    valsS <- valsS[sel,]
                    subX <- rep(start(GdObject)[sel], len) + (subW * rep(seq_len(len)-1, each=sum(sel))) +
                        (spacer * rep(seq_len(len), each=sum(sel)))
                    grid.rect(subX, yy, width=rep(subW, len), height=valsS-yy,
                              gp=gpar(col="transparent", fill=rep(pcols$col[1:len], each=sum(sel)),
                                      lwd=pcols$lwd[1], lty=pcols$lty[1], alpha=alpha), default.units="native",
                              just=c("left", "bottom"))
                }
            }
        } else {
            agFun <- .aggregator(GdObject)
            valsS <- agFun(t(vals))
            grid.rect(start(GdObject), yy, width=width(GdObject), height=valsS-yy,
                      gp=gpar(col=pcols$col.histogram, fill=pcols$fill.histogram, lwd=pcols$lwd[1], lty=pcols$lty[1], alpha=alpha), default.units="native",
                      just=c("left", "bottom"))
        }
    }
    ## gradient summarizes the data as a color gradient
    if("gradient" %in% type)
    {
        ncolor <- .dpOrDefault(GdObject, "ncolor", 100)
        gradient <- colorRampPalette(.dpOrDefault(GdObject, "gradient", brewer.pal(9, "Blues")))(ncolor)
        valsScaled <- .z2icol(colMeans(vals, na.rm=TRUE), ncolor, sort(ylim))
        grid.rect(start(GdObject), sort(ylim)[1], width=width(GdObject), height=abs(diff(ylim)),
                  gp=gpar(col=gradient[valsScaled], fill=gradient[valsScaled], alpha=alpha),
                  default.units="native", just=c("left", "bottom"))
    }
    ## heatmap does the same, but for each sample individually
    if("heatmap" %in% type)
    {
        ncolor <- .dpOrDefault(GdObject, "ncolor", 100)
        valsScaled <- .z2icol(vals, ncolor, sort(ylim))
        nr <- nrow(vals)
        yy <- seq(min(ylim), max(ylim), len=nr+1)[-1]
        ydiff <- .pxResolution(coord="y")
        separator <- .dpOrDefault(GdObject, "separator", 0)*ydiff
        if(!is.null(groups))
        {
            valsS <- split(vals, groups)
            freq <- table(factor(displayPars(GdObject, "groups")))
            cmf <- c(0, cumsum(freq))
            for(s in seq_along(valsS))
            {
                gradient <- colorRampPalette(c("white", pcols$col[s]))(ncolor+5)[-(1:5)]
                valsScaled <- .z2icol(valsS[[s]], ncolor, sort(ylim))
                grid.rect(rep(start(GdObject), each=freq[s]), yy[(cmf[s]+1):cmf[s+1]], width=rep(width(GdObject), each=freq[s]),
                  height=max(ydiff, abs(diff(ylim))*(1/nr)-separator),
                  gp=gpar(col=gradient[valsScaled], fill=gradient[valsScaled], alpha=alpha),
                  default.units="native", just=c("left", "top"))
            }
        } else {
            gradient <- colorRampPalette(.dpOrDefault(GdObject, "gradient", brewer.pal(9, "Blues")))(ncolor)
            grid.rect(rep(start(GdObject), each=nr), yy, width=rep(width(GdObject), each=nr),
                      height=max(ydiff, abs(diff(ylim))*(1/nr)-separator),
                      gp=gpar(col=gradient[valsScaled], fill=gradient[valsScaled], alpha=alpha),
                      default.units="native", just=c("left", "top"))
        }
    }
    ## The rest uses the lattice panel function
    na.rm <- .dpOrDefault(GdObject, "na.rm", FALSE)
    sel <- is.na(y)
    if(na.rm && any(sel))
    {
        x <- x[!sel]
        y <- y[!sel]
        groups <- groups[!sel]
    }
    panel.xyplot(x, y, type=type, groups=groups, pch=pcols$pch, col=pcols$col, col.line=pcols$col.line, col.symbol=pcols$col.symbol,
                 font=font, fontfamily=font, fontface=fontface, lty=pcols$lty, cex=pcols$cex, fill=pcols$fill, lwd=pcols$lwd, horizontal=FALSE,
                 span=span, degree=degree, family=family, evaluation=evaluation,
                 jitter.x=.dpOrDefault(GdObject, "jitter.x", FALSE), jitter.y=.dpOrDefault(GdObject, "jitter.y", FALSE),
                 factor=.dpOrDefault(GdObject, "factor", 0.5), amount=.dpOrDefault(GdObject, "amount"),
                 subscripts=seq_along(x), alpha=alpha)
    if(!"mountain" %in% type && !is.null(baseline) && !is.na(baseline))
        panel.abline(h=baseline, col=pcols$col.baseline, lwd=lwd.baseline, lty=lty.baseline, alpha=alpha)
    popViewport(1)
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------------
## Draw a AlignedRead track
##----------------------------------------------------------------------------------------------------------------------------
setMethod("drawGD", signature("AlignedReadTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
    imageMap(GdObject) <- NULL
    detail <- match.arg(.dpOrDefault(GdObject, "detail", "coverage"), c("reads", "coverage"))
    ## Nothing to do in prepare mode if detail is not 'reads', so we can quit right away, else we need to set the stacking info
    if(prepare){
        if(detail=="read"){
            if(subset)
                GdObject <- subset(GdObject, from=minBase, to=maxBase)
            ##GdObject <- setStacks(GdObject)
        }
        return(invisible(GdObject))
    }
    ## In plotting mode we either show all the reads (time-consuming), or the coverage only
    rad <- 0.015
    xx <- -0.01
    loc <- vpLocation()$size
    diff <- .pxResolution(coord="x")
    radv <- rad / if(loc["width"] < loc["height"]) c(1,loc[2]/loc[1]) else c(loc[1]/loc[2],1)
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    if(!length(GdObject))
    {
        ## No reads, but we still need the strand separator
        panel.abline(h=0.5, col="lightgray", lwd=2)
        grid.circle(xx, c(0.25, 0.75), rad, gp=gpar(fill="lightgray", col="lightgray"))
        grid.segments(c(rep(xx-radv[1]+(radv[1]/2),2), xx), c(0.25, 0.75, 0.75-(radv[2]/2)),
                      c(rep(xx+radv[1]-(radv[1]/2),2), xx), c(0.25, 0.75, 0.75+radv[2]/2),
                      gp=gpar(col="white", lwd=2, lineend="square"), default.units="native")
        return(invisible(GdObject))
    }
    ## If type is 'coverage' all we need to do is compute a coverage vector, create dummy DataTracks and pass everything on
    if(detail=="coverage")
    {
        ## We want to distinguish between strands, so an extra spitting step is needed for this to work
        val <- c(0, max(unlist(sapply(c("+", "-"), function(x) if(length(coverage(GdObject, strand=x)))
                               max(coverage(GdObject, strand=x)) else NULL))))
        trans <- displayPars(GdObject, "transformation")[[1]]
        if(!is.null(trans))
            val[2] <- trans(val[2])
        ylim <- .dpOrDefault(GdObject, "ylim", val)
        for(s in c("+", "-"))
        {
            cov <- coverage(GdObject, strand=s)
            pushViewport(viewport(height=0.5, y=ifelse(s=="-", 0, 0.5), just=c("center", "bottom")))
            sel <- suppressWarnings(runValue(cov)>0)
            dtr <- if(any(sel)) DataTrack(start=start(cov)[sel], end=end(cov)[sel], data=runValue(cov)[sel],
                                          name=names(GdObject), genome=genome(GdObject), chromosome=chromosome(GdObject)) else
            DataTrack(name=names(GdObject), genome=genome(GdObject), chromosome=chromosome(GdObject))
            displayPars(dtr) <- displayPars(GdObject)
            displayPars(dtr) <- list(ylim=if(s=="+") ylim else rev(ylim))
            drawGD(dtr, minBase, maxBase, prepare=prepare, ...)
            popViewport(1)
        }
        panel.abline(h=0.5, col="lightgray", lwd=2)
        grid.circle(xx, c(0.25, 0.75), unit(rad, "native"), gp=gpar(fill="lightgray", col="lightgray"),
                    default.units="native")
        grid.segments(c(rep(xx-radv[1]+(radv[1]/2),2), xx), c(0.25, 0.75, 0.75-(radv[2]/2)),
                      c(rep(xx+radv[1]-(radv[1]/2),2), xx), c(0.25, 0.75, 0.75+radv[2]/2),
                      gp=gpar(col="white", lwd=2, lineend="square"), default.units="native")
        return(invisible(GdObject))
    }
    if(detail=="reads")
    {
        if(GdObject@coverageOnly){
            pushViewport(viewport())
            grid.text("Coverage information only for this object.\nUnable to plot read details.",
                      gp=gpar(col="darkgray"))
            panel.abline(h=0.5, col="lightgray", lwd=2)
            recMid <- c(0.25, 0.75)
        } else {
            gdSplit <- split(GdObject, factor(strand(GdObject), levels=c("+", "-")))
            st <- lapply(gdSplit, function(x) if(length(x)) stacks(setStacks(x))-1 else 0)
            omax <- sum(sapply(st, max))
            space <- 0.1
            ratios <- sapply(st[c("-", "+")], function(x) if(omax==0) 0.5 else max(x)/omax)+(space/(2:1))
            y <- if(ratios[["-"]]<ratios[["+"]]) c(0, ratios[["-"]]+space/2) else c(0, ratios[["-"]]-space/2)
            h <- if(ratios[["-"]]<ratios[["+"]]) c(max(space/2, ratios["-"]-space/2),
                                                   max(space/2, ratios["+"]-ratios["-"]-space/2)) else
            c(max(space/2, ratios["-"]-space), max(space/2, ratios["+"]-ratios["-"]-space/2))
            names(y) <- names(h) <- names(ratios)
            pushViewport(viewport(height=0.95, just="center", yscale=c(0, max(ratios))))
            for(s in c("+", "-"))
            {
                gdSub <- gdSplit[[s]]
                if(length(gdSub))
                {	
                    ylim <- c(0, if(max(st[[s]])==0) 1 else max(st[[s]]))
                    pushViewport(viewport(height=h[s], y=y[s], just=c("center", "bottom"), 
                                          yscale=if(s=="+") ylim else rev(ylim), xscale=c(minBase, maxBase),
                                          default.units="native"))
                    ## We need to handle the grid here individually
                    if(.dpOrDefault(GdObject, "grid", FALSE))
                    {
                        pushViewport(dataViewport(xData=c(minBase, maxBase), extension=c(0, 0), yData=0:1, clip=TRUE))
                        panel.grid(h=0, v=.dpOrDefault(GdObject, "v", -1),
                                   col=.dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty=.dpOrDefault(GdObject, "lty.grid", 1),
                                   lwd=.dpOrDefault(GdObject, "lwd.grid", 1))
                        popViewport(1)
                    }
                    grid.segments(start(gdSub), st[[s]], end(gdSub), st[[s]], default.units="native")
                    popViewport(1)
                }
            }
            al <- h["-"]+space/4
            panel.abline(h=al, col="lightgray", lwd=2)
            recMid <- c(al/2, al+(max(ratios)-al)/2)
        }
        grid.circle(xx, recMid, rad, gp=gpar(fill="lightgray", col="lightgray"), default.units="native")
        grid.segments(c(rep(xx-radv[1]+(radv[1]/2),2), xx), c(0.25, 0.75, 0.75-(radv[2]/2)),
                      c(rep(xx+radv[1]-(radv[1]/2),2), xx), c(0.25, 0.75, 0.75+radv[2]/2),
                      gp=gpar(col="white", lwd=3, lineend="square"), default.units="native")
        
        grid.segments(c(rep(xx-radv[1]+(radv[1]/2),2), xx), c(recMid, recMid[2]-radv[2]/2),
                      c(rep(xx+radv[1]-(radv[1]/2),2), xx), c(recMid, recMid[2]+radv[2]/2),
                      gp=gpar(col="white", lwd=2, lineend="square"), default.units="native")
        popViewport(1)
    }
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------		
		



##----------------------------------------------------------------------------------------------------------------------------
## Draw an ideogram track
##----------------------------------------------------------------------------------------------------------------------------
## Helper function to compute coordinates for a rounded ideogram cap
.roundedCap <- function(bl, tr, side=c("left", "right"), bevel=0.4, n=100)
{
    side <- match.arg(side)
    bevel <- max(1/n, min(bevel, 0.5))
    coords <- if(bevel<=0) cbind(c(0,1,1,0), c(1.5,1.5,-1.5,-1.5)) else
    cbind(c(0, sin(seq(0,pi/2,len=max(2, n*bevel)))+bevel*4, sin(seq(pi/2, pi, len=max(2,  n*bevel)))+bevel*4, 0)/(1+4*bevel),
          (c(1+0.5-bevel, cos(seq(0,pi/2,len=max(2, n*bevel)))+(0.5-bevel),
             cos(seq(pi/2, pi, len=max(2, n*bevel)))-(0.5-bevel), cos(pi)-(0.5-bevel))+1+(0.5-bevel))/(2+2*(0.5-bevel)))
    if(side=="right"){
        coords[,1] <- coords[,1]*abs(diff(c(bl[1], tr[1])))+bl[1]
        coords[,2] <- coords[,2]*abs(diff(c(bl[2], tr[2])))+bl[2]
    }else{
        coords[,1] <- (1-coords[,1])*abs(diff(c(bl[1], tr[1])))+bl[1]
        coords[,2] <- (1-coords[,2])*abs(diff(c(bl[2], tr[2])))+bl[2]
    }
    return(coords)
}

## The actual drawing method
setMethod("drawGD", signature("IdeogramTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...) {
    imageMap(GdObject) <- NULL
    chrnam <- paste("Chromosome", gsub("chr", "", chromosome(GdObject)))
    cex <- .dpOrDefault(GdObject, "cex", 1)
    ## Nothing to do if there are no ranges in the object
    if(!length(GdObject))
        return(invisible(GdObject))
    ## In prepare mode we just want to figure out the optimal size
    if(prepare)
    {
        pres <- .pxResolution()
        nsp <-  if(.dpOrDefault(GdObject, "showId", TRUE))
        {
            (as.numeric(convertHeight(stringHeight(chrnam),"native"))*cex) +
                 as.numeric(convertHeight(unit(20, "points"), "native"))
        } else as.numeric(convertHeight(unit(25, "points"), "native"))
        nsp <- nsp/pres["y"]
        displayPars(GdObject) <- list("neededVerticalSpace"=nsp)
        return(invisible(GdObject))
    }
    ## Do we need some space for the chromosome name?
    if(.dpOrDefault(GdObject, "showId", TRUE))
    {
      
        fontface <- .dpOrDefault(GdObject, "fontface", 1)
        width <- vpLocation()$isize["width"]
        width <- as.numeric(convertWidth(stringWidth(chrnam), "inches"))*cex+0.2
        wfac <- vpLocation()$isize["width"]
        nspace <- min(width/wfac, 0.75)
        if((width/wfac)>0.75)
            cex <- cex*(0.75/(width/wfac))
        pushViewport(viewport(x=0, width=nspace, just=0))
        grid.text(chrnam, 0, gp=gpar(cex=cex, fontface=fontface),
                  default.units="native", just=c("left", "center"))
        popViewport(1)
    } else nspace <- 0
    pushViewport(viewport(x=nspace, width=1-nspace, just=0))
    ## A box indicating the current range on the chromosome
    len <- end(range(range(GdObject)))
    fill <- .dpOrDefault(GdObject, "fill", "#FFE3E6")
    if(!missing(minBase) && !missing(maxBase))
        grid.rect(minBase/len, 0.1, width=min(1,(maxBase-minBase)/len), height=0.8, just=c("left","bottom"),
                  gp=gpar(col="transparent", fill=fill))
    ## Color mapping for the bands, grayscale for regular and red for centrosome/repetetive regions
    gpcols <- unique(grep("gpos", values(GdObject)$type, value=TRUE))
    posCol <- if(length(gpcols)==1) 50 else sort(as.integer(gsub("gpos", "", gpcols)))
    cols <- c("white", colorRampPalette(c("white","black"))(100)[posCol], "black", rep("darkred", 2))
    names(cols) <- c("gneg", paste("gpos", if(length(gpcols)==1) "" else posCol, sep=""), "gvar", "acen", "stalk")
    vals <- data.frame(values(GdObject), col=cols[as.character(values(GdObject)$type)], stringsAsFactors=FALSE)
    ## For the rounded caps we need  to figure out the overlap with existing bands for proper coloring
    bevel <- 0.02
    ol <- queryHits(findOverlaps(range(GdObject), IRanges(start=c(bevel, 1-bevel)*len, width=1)))
    st <- start(range(GdObject))/len
    stExt <- c(st[1:ol[1]], bevel, st[(ol[1]+1):ol[2]], 1-bevel)
    valsExt <- rbind(vals[1:ol[1],], vals[ol[1],], vals[(ol[1]+1):ol[2],], vals[ol[2],])
    if(ol[2]<length(st)){
        stExt <- c(stExt, st[(ol[2]+1):length(st)])
        valsExt <- rbind(valsExt, vals[(ol[2]+1):length(st),])
    }
    wd <- diff(c(stExt,1))
    ls <- ol[1]+1
    rs <- ol[2]+1
    ## The centrosome is treated separately
    cent <- grep("acen", valsExt$type)
    if(length(cent))
    {
        bef <- ls:(min(cent)-1)
        aft <- (max(cent)+1):rs
    }else{
        bef <- ls:rs
        aft <- NULL
    }
    margin <- 0.3
    ## First the normal bands
    grid.rect(stExt[bef], margin, width=wd[bef], height=1-margin*2, gp=gpar(col=valsExt[bef,"col"], fill=valsExt[bef,"col"]),
              just=c("left","bottom"))
    if(!is.null(aft))
        grid.rect(stExt[aft], margin, width=wd[aft], height=1-margin*2, gp=gpar(col=valsExt[aft,"col"],
                                                                                fill=valsExt[aft,"col"]),
                  just=c("left","bottom"))
    ## Now the centrosome, if there is any
    if(length(cent))
    {
        grid.polygon(c(stExt[min(cent)], stExt[min(cent)]+wd[min(cent)], rep(stExt[min(cent)],2)),
                     c(margin, 0.5, (1-margin), margin),
                     gp=gpar(col=cols["acen"], fill=cols["acen"]))
        grid.polygon(c(stExt[max(cent)], rep(stExt[max(cent)]+wd[max(cent)], 2), stExt[max(cent)]),
                     c(0.5, margin, (1-margin), 0.5),
                     gp=gpar(col=cols["acen"], fill=cols["acen"]))
    }
    ## Now the caps
    lc <- .roundedCap(c(stExt[1], margin), c(stExt[ls], 1-margin), side="left", bevel=.dpOrDefault(GdObject, "bevel", 0.45))
    grid.polygon(lc[,1], lc[,2], gp=gpar(col=valsExt[1,"col"], fill=valsExt[1,"col"]))
    rc <- .roundedCap(c(tail(stExt,1), margin), c(1, 1-margin), side="right", bevel=.dpOrDefault(GdObject, "bevel", 0.45))
    grid.polygon(rc[,1], rc[,2], gp=gpar(col=tail(valsExt[,"col"],1), fill=tail(valsExt[,"col"],1)))
    lcol <- "black"; lwd <- 1; lty <- 1
    ## And finally some outlines
    grid.lines(lc[,1], lc[,2], gp=gpar(col=lcol, lwd=lwd, lty=lty))
    grid.lines(rc[,1], rc[,2],gp= gpar(col=lcol, lwd=lwd, lty=lty))
    if(length(cent))
    {
        x0 <- c(rep(stExt[2], 2), rep(stExt[max(cent)+1],2), rep(stExt[min(cent)],2), rep(stExt[max(cent)],2))
        y0 <-  c(rep(c(margin, (1-margin)), 3), 0.5, 0.5)  
        x1 <- c(rep(stExt[min(cent)],2), rep(tail(stExt, 1),2), rep(stExt[min(cent)]+wd[min(cent)],2),
                rep(stExt[max(cent)]+wd[max(cent)],2))
        y1 <-  c(rep(c(margin, (1-margin)), 2), 0.5, 0.5, margin, (1-margin))
    } else {
        x0 <- rep(stExt[2], 2)
        y0 <- c(margin, (1-margin))
        x1 <- rep(max(stExt), 2)
        y1 <- y0
    }
    grid.segments(x0, y0, x1, y1, gp=gpar(col=lcol, lwd=lwd, lty=lty))
    ## Finally the outlines of the box
    if(!missing(minBase) && !missing(maxBase))
    {
        col <- .dpOrDefault(GdObject, "col", "red")
        lwd <- .dpOrDefault(GdObject, "lwd", 1)
        lty <- .dpOrDefault(GdObject, "lty", "solid")
        grid.rect(minBase/len, 0.1, width=min(1,(maxBase-minBase)/len), height=0.8, just=c("left","bottom"),
                  gp=gpar(col=col, fill="transparent", lwd=lwd, lty=lty))
    }
    popViewport(1)
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------------
## Object coercion
##----------------------------------------------------------------------------------------------------------------------------
setAs("AnnotationTrack", "UCSCData",
          function(from, to){
              ranges <- range(from)
              dcolor <-  as.integer(col2rgb(getPar(from, "col")))
              line <- new("BasicTrackLine", name=names(from),
                          description=names(from),
                          visibility=stacking(from), color=dcolor, itemRgb=TRUE)
              vals <- values(from)
              color <- .getBiotypeColor(from)
              strand <- as.character(strand(from))
              strand[strand=="*"] <- "+"
              new("UCSCData", GenomicData(ranges, chrom=chromosome(from),
                                          id=gsub(" ", "_", vals$id),
                                          name=gsub(" ", "_", as.character(vals$id)), itemRgb=color,
                                          genome=genome(from), strand=strand),
                  trackLine = line)   
          })

setAs("GeneRegionTrack", "UCSCData",
          function(from, to){
              ranges <- cbind(as(as(ranges(from), "DataFrame"), "data.frame"), start=start(from),
                              end=end(from), color=.getBiotypeColor(from), strand=strand(from))
              ranges <- ranges[order(start(from)),]
              ranges <- split(ranges, ranges[,"X.transcript"])
              start <- sapply(ranges, function(x) min(x$start))
              end <- sapply(ranges, function(x) max(x$end))
              name <- as.character(sapply(ranges, function(x) unique(x$X.symbol)))
              color <- as.character(sapply(ranges, function(x) unique(x$color)))
              strand <- as.character(sapply(ranges, function(x) unique(x$strand)))
              strand[strand=="*"] <- "+"
              id <- names(ranges)
              blocks <- sapply(ranges, nrow)
              bsizes <- sapply(ranges, function(x) paste(x$end-x$start+1, collapse=","))
              bstarts <- sapply(ranges, function(x) paste(x$start - min(x$start), collapse=","))
              dcolor <-  as.integer(col2rgb(getPar(from, "col")))
              line <- new("BasicTrackLine", name=names(from),
                          description=names(from),
                          visibility=stacking(from), color=dcolor, itemRgb=TRUE)
              new("UCSCData", GenomicData(IRanges(start, end), chrom=chromosome(from),
                                          id=id, name=name, itemRgb=color, blockCount=blocks,
                                          blockSizes=bsizes, blockStarts=bstarts,
                                          genome=genome(from), strand=strand),
                  trackLine = line)   
          })
          
setAs("RangeTrack", "data.frame",
          function(from, to) as(as(ranges(from), "DataFrame"), "data.frame"))

setAs("InferredDisplayPars", "list",
          function(from, to) {ll <- from@.Data; names(ll) <- names(from); ll})

setAs("DisplayPars", "list", function(from, to) as.list(from@pars))

setMethod("as.list", "DisplayPars", function(x) as(x, "list"))

setMethod("as.list", "InferredDisplayPars", function(x) as(x, "list"))

setMethod("head", "InferredDisplayPars", function(x, n=10, ...){
    sel <- 1:min(n, length(x))
    return(new("InferredDisplayPars", name=x@name, inheritance=x@inheritance[sel],
               structure(x@.Data[sel], names=names(x)[sel])))})

setMethod("tail", "InferredDisplayPars", function(x, n=10, ...){
    sel <- max(1, length(x)-n):length(x)
    return(new("InferredDisplayPars", name=x@name, inheritance=x@inheritance[sel],
               structure(x@.Data[sel], names=names(x)[sel])))})

##---------------------------------------------------------------------------------



##---------------------------------------------------------------------------------
## Helper method to build a GRanges object from the input arguments.
##---------------------------------------------------------------------------------
## Coordinates for grouped elements may be passed in as comma-separated values (e.g. "1,5,9"), in which case
## we need to split and convert to numerics. This also implies that the additional arguments (like feature, group, etc.)
## have to be replicated accordingly. We handle this by passing along the repeat vector 'by' to the numeric method below.
setMethod(".buildRange", signature("NULLOrMissing", "FactorOrCharacterOrNULL", "FactorOrCharacterOrNULL",
                                   "FactorOrCharacterOrNULL"),
          function(range, start, end, width, asIRanges=FALSE, ...){
              ## The inputs coordinates are all empty
              if(!length(start) && !length(end))
                  return(if(asIRanges) IRanges() else GRanges())
              delim <- ","
              coords <- c("start", "end", "width")
              lengths <-sapply(coords, function(x) length(get(x)))
              items <- structure(as.list(rep(0, 3)), names=coords)
              by <- NULL
              for(i in coords)
              {
                  val <- get(i)
                  if(!is.null(val))
                  {
                      val <- strsplit(as.character(val), delim)
                      lengths[i] <- length(val)
                      items[[i]] <- listLen(val)
                      val <- unlist(val)
                      assign(i, as.numeric(val))
                  }
              }
              len <- max(lengths)
              if(!all(sapply(items, function(x) length(x)==1 && x==0)))
                  by <- rowMax(matrix(sapply(items, function(x)
                                             if(length(x)==1) rep(x,len) else if(length(x)) x else rep(0,len)), nrow=len))
              return(.buildRange(start=start, end=end, width=width, asIRanges=asIRanges, by=by, len=len, ...))})

## Helper function to handle settings and defaults
.fillWithDefaults <- function(range=as.data.frame(matrix(ncol=0, nrow=len)), defaults, args, len, by=NULL, ignore=NULL)
{
    for(a in setdiff(names(defaults), ignore))
    {
        range[[a]] <- if(is.null(args[[a]])){
            if(is.null(defaults[[a]]))
                stop("The mandatory argument '", a, "' is missing with no default")
            val <- defaults[[a]]
            if(length(val)==1)
                val <- rep(val, len)
            if(!length(val) %in% c(len, sum(by)))
                stop("Number of elements in '", a, "' is invalid")
            if(!is.null(by) && length(val)!=sum(by)) rep(val, by) else val
        } else{
            val <- args[[a]]
            if(length(val)==1)
                val <- rep(val, len)
            if(!length(val) %in% c(len, sum(by)))
                stop("Number of elements in '", a, "' is invalid")
            if(!is.null(by) && length(val)!=sum(by)) rep(val, by) else val
        }
    }
    return(range)
}

## For numeric vectors we can immediately create a data frame after some sanity checking and pass that on to the next method.
setMethod(".buildRange", signature("NULLOrMissing", "NumericOrNULL", "NumericOrNULL", "NumericOrNULL"),
          function(range, start, end, width, asIRanges=FALSE, by=NULL, len, args=list(), defaults=list(), ...){
              range <- data.frame()
              ## Some of the arguments are mutually exclusive and we want to catch this here.
              if(is.null(width))
              {
                  if (is.null(start) || is.null(end))
                      stop("Must specify either start and end or width")
              } else {
                  if(is.null(end) && !is.null(start))
                      end <- as.integer(start)+as.integer(width)
                  else if(is.null(start) && !is.null(end))
                      start <- as.integer(end)-as.integer(width)
                  else stop("Can't pass all three of 'start', 'end' and 'width'")
              }
              if(length(start)!=1 && length(end)!=1 && length(start) != length(end))
                  stop("Start and end must be vectors of the same length")
              if(missing(len))
                  len <- length(start)
              if(length(start) > 0)
              {
                  if(asIRanges)
                      return(IRanges(start=as.integer(start), end=as.integer(end)))
                  range <- .fillWithDefaults(data.frame(start=as.integer(start), end=as.integer(end)), defaults, args, len, by)
              }
              return(.buildRange(range=range, asIRanges=asIRanges, args=list(), defaults=list(), ...))})

## For data.frames we need to check for additional arguments (like feature, group, etc.), the chromosome information
## and create the final GRanges object
setMethod(".buildRange", signature("data.frame"),
          function(range, asIRanges=FALSE, args=list(), defaults=list(), chromosome, ...){
              if(asIRanges){
                  range <- .fillWithDefaults(range, defaults, args, len=nrow(range), ignore=setdiff(names(defaults), c("start", "end")))
                  return(IRanges(start=as.integer(range$start), end=as.integer(range$end)))
              }
              mandArgs <- c("start", "end", names(defaults))
              missing <- setdiff(mandArgs, colnames(range))
              range <- .fillWithDefaults(range, defaults[missing], args[missing], len=nrow(range))
              snames <- if(missing(chromosome)){ if(is.null(range$id)) "dummy" else range$id } else chromosome
              grange <- GRanges(ranges=IRanges(start=range$start, end=range$end), strand=range$strand, seqnames=snames)
              elementMetadata(grange) <- range[,setdiff(colnames(range), c("start", "end", "strand"))]
              return(grange)})


## For GRanges we just need to check for the existence of additional arguments (like feature, group, etc.)
setMethod(".buildRange", signature("GRanges"),
          function(range, asIRanges=FALSE, args=list(), defaults=list(), ...){
              if(asIRanges)
                  return(ranges(range))
              if(length(range))
              {
                  mandArgs <- setdiff(names(defaults), "strand")
                  vals <- values(range)
                  for(a in mandArgs)
                  {
                      if(is.null(vals[[a]]))
                          vals[[a]] <-  if(is.null(defaults[[a]])) stop("The mandatory column '", a, "' is missing with no default") else {
                              defaults[[a]]}
                  }
                  elementMetadata(range) <- vals
              }
              return(range)})

## For IRanges we need to deal with additional arguments (like feature, group, etc.) and create the final GRanges object
setMethod(".buildRange", signature("IRanges"),
          function(range, asIRanges=FALSE, args=list(), defaults=list(), chromosome, strand, ...){
              if(asIRanges)
                  return(range)
              range <- GRanges(seqnames=if(missing(chromosome)) "dummy" else chromosome, range=range,
                               strand=if(!is.null(args$strand)) args$strand else "*")
              if(length(range))
              {
                  vals <- .fillWithDefaults(defaults=defaults, args=args, len=(length(range)), by=NULL, ignore="strand")
                  elementMetadata(range) <- vals
              }
              return(range)})
##---------------------------------------------------------------------------------



##---------------------------------------------------------------------------------
## Interact with ImageMap objects
##---------------------------------------------------------------------------------
setMethod("coords", "NULL", function(ImageMap) NULL)
setMethod("coords", "ImageMap", function(ImageMap) ImageMap@coords)
setMethod("coords", "GdObject", function(ImageMap) coords(imageMap(ImageMap)))
setMethod("tags", "NULL", function(ImageMap) NULL)
setMethod("tags", "ImageMap", function(ImageMap) ImageMap@tags)
setMethod("tags", "GdObject", function(ImageMap) tags(imageMap(ImageMap)))
##---------------------------------------------------------------------------------



##---------------------------------------------------------------------------------
## Show methods for the various classes
##---------------------------------------------------------------------------------
setMethod("show",signature(object="DataTrack"),
          function(object){
              cat(sprintf(paste("Data track '%s' at %i position%s containing %i sample%s mapping",
                                "to chromosome %s on the %s strand of the %s genome:\n"),
                          names(object), length(object),
                          ifelse(length(object)==1, "", "s"),
                          nrow(values(object)),
                          ifelse(nrow(values(object))==1, "", "s"),
                          gsub("^chr", "", chromosome(object)),
                          strand(object)[1], genome(object)), "\n")
              print(ranges(object)[chromosome(object) == seqnames(object)])
          })

setMethod("show",signature(object="AnnotationTrack"),
          function(object){
              cat(sprintf(paste("Annotation track '%s' containing %i item%s and mapping",
                                "to chromosome %s of the %s genome:\n"),
                          names(object), length(object),
                          ifelse(length(object)==1, "", "s"),
                          gsub("^chr", "", chromosome(object)),
                          genome(object)), "\n")
              print(ranges(object)[chromosome(object) == seqnames(object)])
          })

setMethod("show",signature(object="GeneRegionTrack"),
          function(object) {
              cat(sprintf(paste("Gene region '%s' ranging from bp %i to bp %i ",
                                "of chromosome %s of the %s genome containing %i object%s:\n"),
                          names(object), object@start, object@end,
                          gsub("^chr", "", chromosome(object)),
                          genome(object), length(object),
                          ifelse(length(object)==1, "", "s")))
              print(ranges(object)[chromosome(object) == seqnames(object)])
          })

setMethod("show",signature(object="GenomeAxisTrack"),
          function(object) {
              cat(sprintf("Genome axis '%s'\n", names(object)))
              if(.dpOrDefault(object, "add53", FALSE))
                  cat("5->3 label is set\n")
              if(.dpOrDefault(object, "add35", FALSE))
                  cat("3->5 label is set\n")
              if(.dpOrDefault(object, "littleTicks", FALSE))
                  cat("littleTicks label is set\n")
              if(length(object))
              {
                  cat("There are annotated axis regions:\n")
                  print(ranges(object))
              }
          })

setMethod("show",signature(object="IdeogramTrack"),
		  function(object){
			  cat(sprintf(paste("Ideogram track '%s' for chromosome %s of the %s genome"),
							  names(object),
							  gsub("^chr", "", chromosome(object)),
							  genome(object)), "\n")
                      })
  
setMethod("show",signature(object="AlignedReadTrack"),
		  function(object){
			  cat(sprintf(paste("AlignedRead track '%s' containing %i read%s all mapping",
									  "to chromosome %s of the %s genome:\n"),
							  names(object), length(object),
							  ifelse(length(object)==1, "", "s"),
							  gsub("^chr", "", chromosome(object)),
							  genome(object)), "\n")
			  print(ranges(object))
		  })
  
setMethod("show", "DisplayPars", function(object) {
    cat("Display parameters:\n")
    for(i in ls(object@pars))
    {
        cat(i, " = ", sep="")
        o <- try(as.character(object@pars[[i]]), silent=TRUE)
        if(is(o, "try-error")) print(object@pars[[i]]) else cat(o, "\n")
    }
})

setMethod("show", "InferredDisplayPars", function(object){
    cat("\nThe following display parameters are available for '", object@name, "' objects:\n",
        "(see ? ", object@name, " for details on their usage)\n\n", sep="")
    for(i in names(object))
    {
        cat(i, ifelse(object@inheritance[i]==object@name, "",
                      paste(" (inherited from class '", object@inheritance[i], "')", sep="")),
            ": ", sep="")
        if(is.null(object[[i]]))
            cat("NULL\n") else {
                 o <- try(as.character(object[[i]]), silent=TRUE)
                 if(is(o, "try-error")) print(object[[i]]) else cat(o, "\n")
             }
    }
})
##---------------------------------------------------------------------------------
