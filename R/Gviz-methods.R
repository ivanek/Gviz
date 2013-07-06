##----------------------------------------------------------------------------------------------------------------------------
## Some rather general accessors to extract information from all kinds of GdObjects
##----------------------------------------------------------------------------------------------------------------------------
## Extract the full GRanges object from the range slot of an object inheriting from RangeTrack
setMethod("ranges", "RangeTrack", function(x) x@range)
setReplaceMethod("ranges", "RangeTrack", function(x, value) {
    x@range <- value
    return(x)})
setMethod("ranges", "GenomeAxisTrack", function(x) x@range)
setReplaceMethod("ranges", "GenomeAxisTrack", function(x, value) {
    x@range <- value
    return(x)})


## Extract the IRanges part of the GRanges object from the range slot of an object inheriting from RangeTrack
setMethod("range", "RangeTrack", function(x) ranges(x@range))
setMethod("range", "GenomeAxisTrack", function(x) ranges(x@range))

## seqnames, levels and infofrom the range track
setMethod("seqnames", "RangeTrack", function(x) as.character(seqnames(ranges(x))))
setMethod("seqnames", "SequenceDNAStringSetTrack", function(x) as.character(names(x@sequence)))
setMethod("seqnames", "SequenceBSgenomeTrack", function(x) as.character(seqnames(x@sequence)))
setMethod("seqlevels", "RangeTrack", function(x) unique(seqnames(x)))
setMethod("seqlevels", "SequenceDNAStringSetTrack", function(x) seqnames(x)[width(x@sequence)>0])
setMethod("seqlevels", "SequenceBSgenomeTrack", function(x) seqnames)
setMethod("seqinfo", "RangeTrack", function(x) table(seqnames(x)))

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
setMethod("start", "SequenceTrack", function(x) NULL)
setMethod("end", "IdeogramTrack", function(x) NULL)
setMethod("end", "SequenceTrack", function(x) NULL)
setMethod("width", "IdeogramTrack", function(x) NULL)
setMethod("width", "SequenceTrack", function(x) NULL)

## Return the number of individual annotation items (independent of any grouping) in a RangeTrack
setMethod("length", "RangeTrack", function(x) sum(seqnames(x) == chromosome(x)))
setMethod("length", "GenomeAxisTrack", function(x) length(ranges(x)))
setMethod("length", "IdeogramTrack", function(x) length(ranges(x)))
setMethod("length", "SequenceTrack", function(x)
          if(chromosome(x) %in% seqnames(x)) length(x@sequence[[chromosome(x)]]) else 0)
## setMethod("length", "ReferenceAnnotationTrack", function(x) 0)
## setMethod("length", "ReferenceGeneRegionTrack", function(x) 0)
## setMethod("length", "ReferenceDataTrack", function(x) 0)

## Extract the elementMetadata slot from the GRanges object of an object inheriting from RangeTrack as a data.frame.
## For a DataTrack object these values are stored as a numeric matrix in the data slot, and we return this instead.
setMethod("values", "RangeTrack", function(x) as.data.frame(values(ranges(x)), stringsAsFactors=FALSE))
setMethod("values", "GenomeAxisTrack", function(x) as.data.frame(values(ranges(x)), stringsAsFactors=FALSE))
setMethod("values", "DataTrack", function(x, all=FALSE){
    if(sum(dim(x@data))==0) x@data else{
        sel <- if(all) rep(TRUE, ncol(x@data)) else seqnames(x) == chromosome(x)
        x@data[,sel, drop=FALSE]
    }
})
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


## Extract a subsequence from a SequenceTrack. For performance reasons we restrict this to a maximum
## of one million nucleotides (which is already more than plenty...)
setMethod("subseq", "SequenceTrack", function(x, start=NA, end=NA, width=NA){
    padding <- "-"
    if(!is.na(start[1]+end[1]+width[1])){
        warning("All 'start', 'stop' and 'width' are provided, ignoring 'width'")
        width <- NA
    }
    ## We want start and end to be set if width is provided
    if(!is.na(width[1])){
        if(is.na(start) && is.na(end))
            stop("Two out of the three in 'start', 'end' and 'width' have to be provided")
        if(is.na(start))
            start <- end-width[1]+1
        if(is.na(end))
            end <- start+width[1]-1
    }
    w <- length(x)
    if(is.na(start))
        start <- 1
    if(w>0){
        if(is.na(end))
            end <- w
        rstart <- max(1, start[1], na.rm=TRUE)
        rend <- max(rstart, min(end[1], w, na.rm=TRUE))
    }else{
        if(is.na(end))
            end <- start
        rend <- end
        rstart <- start
    }
    if(rend<rstart)
        stop("'end' has to be bigger than 'start'")
    if((rend-rstart+1)>10e6)
        stop("Sequence is too big! Unable to extract")
    finalSeq <- rep(DNAString(padding), end-start+1)
    if(chromosome(x) %in% seqnames(x) && rend>rstart){
        chrSeq <- x@sequence[[chromosome(x)]]
        seq <- subseq(chrSeq, start=rstart, end=rend)
        if(is(x, "SequenceBSgenomeTrack")) seq <- unmasked(seq)
        subseq(finalSeq, ifelse(start<1, abs(start)+2, 1), width=rend-rstart+1) <- seq
    }
    if(is(x, "SequenceBSgenomeTrack") && chromosome(x) %in% seqnames(x))
        x@pointerCache[[chromosome(x)]] <- x@sequence[[chromosome(x)]]
    if(.dpOrDefault(x, "complement", FALSE))
        finalSeq <- complement(finalSeq)
    return(finalSeq)
})

setMethod("subseq", "ReferenceSequenceTrack", function(x, start=NA, end=NA, width=NA){
    ## We want start and end to be set if width is provided
    if(!is.na(width[1])){
        if(is.na(start) && is.na(end))
            stop("Two out of the three in 'start', 'end' and 'width' have to be provided")
        if(is.na(start))
            start <- end-width[1]+1
        if(is.na(end))
            end <- start+width[1]-1
    }
    x@sequence <- x@stream(file=x@reference, selection=GRanges(chromosome(x), ranges=IRanges(start, end)))
    return(callNextMethod(x=x, start=start, end=end, width=width))
})



## Set or extract the chromosome from a RangeTrack object
setMethod("chromosome", "GdObject", function(GdObject) return(NULL))
setMethod("chromosome", "RangeTrack", function(GdObject) GdObject@chromosome)
setMethod("chromosome", "SequenceTrack", function(GdObject) GdObject@chromosome)
setReplaceMethod("chromosome", "GdObject", function(GdObject, value){
    return(GdObject)
})
setReplaceMethod("chromosome", "RangeTrack", function(GdObject, value){
    GdObject@chromosome <- .chrName(value[1])
    return(GdObject)
})
setReplaceMethod("chromosome", "SequenceTrack", function(GdObject, value){
    GdObject@chromosome <- .chrName(value[1])
    return(GdObject)
})
setReplaceMethod("chromosome", "IdeogramTrack", function(GdObject, value){
    ## We have changed the class definition to include the bands for all chromosomes, but still want the old objects to work
    chromosome <- .chrName(value[1])
    if(.hasSlot(GdObject, "bandTable") && chromosome %in% as.character(GdObject@bandTable$chrom))
    {
        ranges <- GdObject@bandTable[GdObject@bandTable$chrom==chromosome,]
        ranges <- GRanges(seqnames=as.vector(ranges$name), ranges=IRanges(start=ranges$chromStart, end=ranges$chromEnd),
                          name=as.vector(ranges$name), type=ranges$gieStain)
        GdObject@range <- ranges
        GdObject@chromosome <- chromosome
        return(GdObject)
    }
    message("Updating chromosome band information")
    tmp <- IdeogramTrack(genome=genome(GdObject), chromosome=.chrName(value[1]), name=names(GdObject))
    displayPars(tmp) <- displayPars(GdObject)
    return(tmp)
})
setMethod("isActiveSeq", "RangeTrack", function(x) chromosome(x))
setReplaceMethod("isActiveSeq", "GdObject", function(x, value){
    chromosome(x) <- value
    return(x)
})

## Set or extract the genome from a RangeTrack object
setMethod("genome", "RangeTrack", function(x) x@genome)
setMethod("genome", "SequenceTrack", function(x) x@genome)
setReplaceMethod("genome", "GdObject", function(x, value){
    return(x)
})
setReplaceMethod("genome", "RangeTrack", function(x, value){
    x@genome <- value[1]
    genome(ranges(x)) <- as.vector(value[1])
    return(x)
})
setReplaceMethod("genome", "IdeogramTrack", function(x, value){
    if(genome(x)!=value)
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
setMethod("[", signature(x="StackedTrack"), function(x, i) {
    x <- callNextMethod(x,i)
    x@stacks <- x@stacks[i]
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
        if(ncol(x@data)>0)
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
## annotation is returned, i.e., exon ids  for GeneRegionTracks and id for AnnotationTracks. When using the identifiers as group
## labels we want to add some white space to separate the label from the last item in the group. This can be done by setting
## add.space=TRUE (is ignored if lowest==TRUE)
setMethod("identifier", "AnnotationTrack", function(GdObject, lowest=FALSE, add.space=FALSE){
    id <- if(lowest) .getAnn(GdObject, "id") else group(GdObject)
    if(!lowest && add.space && .dpOrDefault(GdObject, ".__hasAnno", TRUE))
        id[id!=""] <- paste(id[id!=""], "  ", sep="")
    return(id)
})
setMethod("identifier", "GeneRegionTrack", function(GdObject, lowest=FALSE, add.space=FALSE){
    id <- if(lowest) exon(GdObject) else if(.dpOrDefault(GdObject, "geneSymbols", TRUE)){
        symbol(GdObject)} else gene(GdObject)
    id[is.na(id)] <- "NA"
    if(!lowest && add.space && .dpOrDefault(GdObject, ".__hasAnno", TRUE))
        id[id!=""] <- paste(id[id!=""], "  ", sep="")
    return(id)
})        
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
                              paste(pt@stackingValues, collapse=", "), "\n")
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
##   o AnnotationTrack and GeneRegionTrack: features can be grouped, in which case we have to avoid overlapping of the whole group region,
##      i.e, from the start of the first group item to the end of the last. In addition to grouping we have to factor in additional space
##      before each group if group/transcript annotation is enabled (gpar showId==TRUE). To do so we need to figure out the current
##      fontsize on the device, which means that a device already has to be open and the appropriate viewport has been
##      pushed to the stack. Hence we have to call setStacks immediately before the actual plotting.
## 'stacks' should return a vector of stacks, where each factor level of the vector indicates membeship
## to a particular stacking level. 'setStacks' returns the updated GdObject.
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
    if(!.needsStacking(GdObject) || length(GdObject)==0)
    {
        bins <- rep(1, length(GdObject))
    } else {
        uid <- if(is(GdObject, "GeneRegionTrack") && .dpOrDefault(GdObject, "collapseTranscripts", FALSE))
            sprintf("uid%i", seq_along(identifier(GdObject))) else make.unique(identifier(GdObject, lowest=TRUE))
        gp <- group(GdObject)
        needsGrp <- any(duplicated(gp))
        gRanges <- if(!length(gp)){
            IRanges()
        }else{
            if(needsGrp){
                groups <- split(range(GdObject), gp)
                uidSplit <- split(uid, gp)
                unlist(range(groups))
            }else{
                range(GdObject)
            }
        }
        if(.dpOrDefault(GdObject, ".__hasAnno", FALSE))
        {
            cex <- .dpOrDefault(GdObject, "cex", 1) * .dpOrDefault(GdObject, "cex.symbol", 0.7)
            fontfamily <- .dpOrDefault(GdObject, "fontfamily", 1)
            fontsize <- .dpOrDefault(GdObject, "fontsize", 12)
            fontface <- .dpOrDefault(GdObject, "fontface.symbol", 2)
            pushViewport(dataViewport(xData=c(from, to), extension=0, yscale=c(0, 40), clip=TRUE,
                                      gp=gpar(cex=cex, fontfamily=fontfamily, fonface=fontface, fontsize=fontsize)))
            if(needsGrp)
            {
                ids <- sapply(split(identifier(GdObject, add.space=TRUE), gp), head, 1)
                start(gRanges) <- start(gRanges)-(as.numeric(convertWidth(stringWidth(ids),"native"))*1.3)
            } else {
                start(gRanges) <- start(gRanges)-(as.numeric(convertWidth(stringWidth(identifier(GdObject, add.space=TRUE)),"native"))*1.3)
            }
            popViewport(1)
        } 
        bins <- if(needsGrp) rep(disjointBins(gRanges), sapply(groups, length)) else disjointBins(gRanges)
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
setMethod("position", signature("RangeTrack"), definition=function(GdObject, from=NULL, to=NULL, sort=FALSE, ...)
      {
          if(!is.null(from) && !is.null(to))
              GdObject <- subset(GdObject, from=from, to=to, sort=sort, ...)
          pos <- if(length(GdObject)) rowMeans(cbind(start(GdObject), end(GdObject))) else numeric()
          return(pos)
      })
setMethod("position", signature("IdeogramTrack"), definition=function(GdObject, ...) NULL)
     
## The numeric values of data tracks
setMethod("score", signature("DataTrack"), function(x, from=NULL, to=NULL, sort=FALSE, transformation=TRUE, ...)
      {
          if(!is.null(from) && !is.null(to))
              x <- subset(x, from=from, to=to, sort=sort, ...)
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



##----------------------------------------------------------------------------------------------------------------------------
## Before starting of the plotting operation there are a bunch of housekeeping task that should be performed on each
## track, and the mileage may vary between track types, hence we add a layer of abstraction here by using a method.
## Available arguments are:
##    o GdObject: the input track object
##    o chromosome: the currently active chromosome which may have to set for a RangeTrack object
##    o ...: additional arguments that are considered to be display parameters
##----------------------------------------------------------------------------------------------------------------------------
## For all track types we want to update the display parameters
setMethod("consolidateTrack", signature(GdObject="GdObject"), function(GdObject, ...) {
    pars <- list(...)
    pars <- pars[names(pars)!=""]
    displayPars(GdObject) <- pars
    return(GdObject)
})
## For RangeTracks and SequenceTracks we want to set the chromosome
setMethod("consolidateTrack", signature(GdObject="RangeTrack"), function(GdObject, chromosome, ...) {
    if(!is.null(chromosome))
        chromosome(GdObject) <- chromosome
    GdObject <- callNextMethod()
    return(GdObject)
})
setMethod("consolidateTrack", signature(GdObject="SequenceTrack"), function(GdObject, chromosome, ...) {
    if(!is.null(chromosome))
        chromosome(GdObject) <- chromosome
    GdObject <- callNextMethod()
    return(GdObject)
})
## For StackedTracks we want to set the stacking (which could have been passed in as a display parameter)
setMethod("consolidateTrack", signature(GdObject="StackedTrack"), function(GdObject, ...) {
    GdObject <- callNextMethod()
    st <- displayPars(GdObject, "stacking")
    if(!is.null(st))
        stacking(GdObject) <- st
    return(GdObject)
})
## For AnnotationTracks we need to determine whether there is group label annotation or not
setMethod("consolidateTrack", signature(GdObject="AnnotationTrack"), function(GdObject, ...) {
    GdObject <- callNextMethod()
    ids <- identifier(GdObject)
    hasAnno <- .dpOrDefault(GdObject, "showId", FALSE) & !all(ids=="")
    displayPars(GdObject) <- list(".__hasAnno"=hasAnno)
    return(GdObject)
})
##----------------------------------------------------------------------------------------------------------------------------


              
##----------------------------------------------------------------------------------------------------------------------------
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
##    o GdObject: the input AnnotationTrack or GeneRegionTrack track object
##    o min.width: the minimum width in pixel, everything that's smaller will be expanded to this size.
##    o min.distance: the minimum distance between two features in pixels below which to start collapsing track items.
##    o collapse: logical, collapse overlapping items into a single meta-item.
##    o diff: the equivalent of 1 pixel in the native coordinate system.
##    o xrange: the data range on the x axis. Can be used for some preliminary subsetting to speed things up
##----------------------------------------------------------------------------------------------------------------------------
## A slightly quicker function to compute overlaps between two GRanges objects
.myFindOverlaps <- function(gr1, gr2)
{
    gr1 <- sort(gr1)
    gr2 <- sort(gr2)
    gr1 <- split(gr1, seqnames(gr1))
    gr2 <-split(gr2, seqnames(gr2))
    queryHits(findOverlaps(ranges(gr1), ranges(gr2)))
}
## Find all elements in a GRanges object 'grange' width distance smaller than 'minXDist' and merge them along with their additional
## elementMetadata. 'elements' is a frequency table of items per group, and it is needed to figure out whether all items of a given
## group have been merged. 'GdObject' is the input tack object from which certain information has to be extracted. The output of this
## function is a list with elements
##   o range: the updated GRanges object
##   o needsRestacking: logical flag indicating whether stacks have to be recomputed
##   o split: the original merged annotation (need to workon this, currently not needed)
##   o merged: logical vector indicating which of the elements in 'range' constitute fully merged groups 
.collapseAnnotation <- function(grange, minXDist, elements, GdObject, offset=0)
{
    needsRestacking <- TRUE
    annoSplit <- merged <- NULL
    anno <- as.data.frame(grange)
    for(i in colnames(anno))
        if(is.factor(anno[,i]))
            anno[,i] <- as.character(anno[,i])
    cols <- c("strand", "density", "gdensity", "feature", "id", "start", "end", if(is(GdObject, "GeneRegionTrack"))
              c("gene", "exon", "transcript", "symbol", "rank") else "group")
    missing <- which(!cols %in% colnames(anno))
    for(i in missing)
        anno[,cols[missing]] <- if(cols[i]=="density") 1 else NA
    rRed <- if(length(grange)>1) reduce(grange, min.gapwidth=minXDist) else grange
    if(length(rRed) < length(grange))
    {
        ## Some of the items have to be merged and we need to make sure that the additional annotation data that comes with it
        ## is processed in a sane way.
        needsRestacking <- TRUE
        ##mapping <- .myFindOverlaps(rRed, grange)
        mapping <- queryHits(findOverlaps(rRed, grange))
        ## We start by finding the items that have not been reduced
        identical <- mapping %in% which(table(mapping)==1)
        newVals <- anno[identical,cols]
        ## Here we hijack the seqnames column to indicate whether the whole group has been merged
        if(nrow(newVals)){
            newVals$seqnames <- elements[as.character(anno[identical,"seqnames"])]==1
            newVals$gdensity <- ifelse(elements[as.character(anno[identical,"seqnames"])]==1, 1, NA)
        }
        ## Now find out which original items have been merged
        grange <- grange[!identical]
        rRed <- rRed[-(mapping[identical])]
        index <- mapping[!identical]
        annoSplit <- split(anno[!identical,], index)
        cid <- function(j) sprintf("[Cluster_%i]  ", j+offset)
        ## FIXME: We could speed this up by running it in C
        newVals <- rbind(newVals, as.data.frame(t(sapply(seq_along(annoSplit), function(i){
            x <- annoSplit[[i]]
            if(is(GdObject, "GeneRegionTrack")){
                c(strand=ifelse(length(unique(x[,"strand"]))==1, as.character(x[1,"strand"]), "*"),
                  density=sum(as.integer(x[,"density"])),
                  gdensity=ifelse(is.na(head(x[,"gdensity"], 1)), 1, sum(as.integer(x[,"gdensity"]))),
                  feature=ifelse(length(unique(x[,"feature"]))==1, as.character(x[1,"feature"]), "composite"),
                  id=ifelse(length(unique(x[,"id"]))==1, as.character(x[1,"id"]), cid(i)),
                  start=min(x[,"start"]),
                  end=max(x[,"end"]),
                  gene=ifelse(length(unique(x[,"gene"]))==1, as.character(x[1,"gene"]), cid(i)),
                  exon=ifelse(length(unique(x[,"exon"]))==1, as.character(x[1,"exon"]), cid(i)),
                  transcript=ifelse(length(unique(x[,"transcript"]))==1, as.character(x[1,"transcript"]), cid(i)),
                  symbol=ifelse(length(unique(x[,"symbol"]))==1, as.character(x[1,"symbol"]), cid(i)),
                  rank=min(as.integer(x[,"rank"])), seqnames=as.vector(nrow(x)==elements[x[1,"seqnames"]]))
            }else{
                c(strand=ifelse(length(unique(x[,"strand"]))==1, as.character(x[1,"strand"]), "*"),
                  density=sum(as.integer(x[,"density"])),
                  gdensity=ifelse(is.na(head(x[,"gdensity"], 1)) , 1, sum(as.integer(x[,"gdensity"]))),
                  feature=ifelse(length(unique(x[,"feature"]))==1, as.character(x[1,"feature"]), "composite"),
                  id=ifelse(length(unique(x[,"id"]))==1, as.character(x[1,"id"]), cid(i)),
                  start=min(x[,"start"]),
                  end=max(x[,"end"]),
                  group=ifelse(length(unique(x[,"group"]))==1, as.character(x[1,"group"]), cid(i)),
                  seqnames=as.vector(nrow(x)==elements[x[1,"seqnames"]]))
            }
        })), stringsAsFactors=FALSE))
        merged <- as.logical(newVals$seqnames)
        grange <- GRanges(seqnames=chromosome(GdObject), strand=newVals[, "strand"],
                          ranges=IRanges(start=as.integer(newVals[, "start"]), end=as.integer(newVals[, "end"])))
        cnMatch <- match(c(colnames(values(GdObject)), "gdensity"), colnames(newVals))
        elementMetadata(grange) <-
            if(any(is.na(cnMatch))) newVals[, setdiff(colnames(newVals), c("strand", "start", "end", "seqnames"))] else newVals[, cnMatch]
    }else{
        grange2 <-  GRanges(seqnames=chromosome(GdObject), strand=strand(grange), ranges=ranges(grange))
        elementMetadata(grange2) <- elementMetadata(grange)
        grange <- grange2
    }
    return(list(range=grange, needsRestacking=needsRestacking, split=annoSplit, merged=merged, offset=length(annoSplit)))
}

## For AnnotationTracks we need to collapse the all regions along with the additional annotation.
## For GeneRegionTracks we essentially need to do the same thing as for AnnotationTracks, however the additional annotation columns
## are quite different. We do this in multiple turn with increasing levels of complexity:
##    1.) merge all individual items within a group that can no longer be separated
##    2.) merge overlapping groups with just a single remaining item (optional, if mergeGroups==TRUE)
setMethod("collapseTrack", signature(GdObject="AnnotationTrack"),        
          function(GdObject, diff=.pxResolution(coord="x"), xrange) {
              ## We first add the original unmodified GdObject as a display parameter to be able to reference back if we ever need to
              displayPars(GdObject) <- list(".__OriginalGdObject"=.deepCopyPars(GdObject))
              collapse <- .dpOrDefault(GdObject, "collapse", TRUE)
              min.width <- .dpOrDefault(GdObject, "min.width", 2)
              min.distance <- .dpOrDefault(GdObject, "min.distance", 2)
              minXDist <- max(0, ceiling(min.distance*diff))
              r <- ranges(GdObject)
              ## Compute native coordinate equivalent to 1 pixel and resize
              rNew <- .resize(r, min.width, diff)
              needsRestacking <- any(r!=rNew)
              r <- rNew
              ## Collapse all items within a group to a single meta-item if (if collapseTranscripts==TRUE)
              if(is(GdObject, "GeneRegionTrack") && .dpOrDefault(GdObject, "collapseTranscripts", FALSE)){
                  newVals <- unlist(endoapply(split(values(r), paste(gene(GdObject), strand(GdObject))), head, 1))
                  newVals$exon <- NA
                  newVals$transcript <- newVals$gene
                  r <- unlist(range(split(r, gene(GdObject))))
                  elementMetadata(r) <- newVals
                  GdObject@range <- r
              }
              ## Collapse overlapping ranges (less than minXDist space between them) and process the annotation data
              if(collapse)
              {
                  ## Merge all items in those groups for which no individual items can be separated
                  elements <- table(group(GdObject))
                  rr <- GRanges(seqnames=as.character(group(GdObject)), ranges=IRanges(start=start(r), end=end(r)), strand=strand(r))
                  elementMetadata(rr) <- elementMetadata(r)
                  rr <- sort(unique(rr))
                  mergedAnn <- .collapseAnnotation(rr, minXDist, elements, GdObject)
                  needsRestacking <- needsRestacking || mergedAnn$needsRestacking
                  ## Now we take a look whether there are any groups that could be merged (if mergeGroups is TRUE)
                  if(.dpOrDefault(GdObject, "mergeGroups", FALSE) && any(mergedAnn$merged)){
                      rr <- sort(mergedAnn$range[mergedAnn$merged])
                      strand(rr) <- "*"
                      mergedAnn2 <- .collapseAnnotation(rr, minXDist, elements, GdObject, mergedAnn$offset)
                      needsRestacking <- needsRestacking || mergedAnn2$needsRestacking
                      mergedAnn$range <- c(mergedAnn$range[!mergedAnn$merged], mergedAnn2$range)
                  }
                  r <- mergedAnn$range
              }
              ## Reconstuct the track object and return
              GdObject@range <- r
              ##if(needsRestacking)
              GdObject <- setStacks(GdObject, xrange[1], xrange[2])
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
                                    suppressWarnings(runValue(runmean(Rle(as.numeric(rm)), k=windowSize, endrule="constant", na.rm=TRUE)))[seqSel]}), 
                                       nrow=nrow(sc)-1, byrow=TRUE))
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
    ## If groups need to be averaged we have to do it here
    groups <- .dpOrDefault(GdObject, "groups")
    if(!is.null(groups) && .dpOrDefault(GdObject, "aggregateGroups", FALSE)){
        if(!is.factor(groups))
            groups <- factor(groups)
        agFun <- .aggregator(GdObject)
        dat <- values(GdObject)
        rownames(dat) <- groups
        datNew <- matrix(t(sapply(levels(groups), function(x) agFun(t(dat[groups==x,,drop=FALSE])), USE.NAMES=FALSE)), nrow=nlevels(groups))
        GdObject@data <- datNew
        displayPars(GdObject) <- list(groups=levels(groups))
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

## For normal ranges we clip everything outside of the boundaries (keeping one extra item left and right
## in order to assure continuation)
setMethod("subset", signature(x="RangeTrack"), function(x, from=NULL, to=NULL, sort=FALSE, drop=TRUE, use.defaults=TRUE, ...){
    ## Not needed anymore...
    ## Subset to a single chromosome first
    if(drop){
        csel <- seqnames(x) != chromosome(x)
        if(any(csel))
            x <- x[,!csel]
    }
    if(!length(x))
        return(x)
    ranges <- if(use.defaults) .defaultRange(x, from=from, to=to) else c(from=ifelse(is.null(from), -Inf, from), to=ifelse(is.null(to), Inf, to))
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

## For DataTracks we cut exactly, and also reduce to the current chromosome unless told explicitely not to
setMethod("subset", signature(x="DataTrack"), function(x, from=NULL, to=NULL, sort=FALSE, drop=TRUE, use.defaults=TRUE, ...){
    ## Subset to a single chromosome first
    if(drop){
        csel <- seqnames(x) != chromosome(x)
        if(any(csel))
            x <- x[,!csel]
    }
    if(!length(x))
        return(x)
    ranges <- if(use.defaults) .defaultRange(x, from=from, to=to) else c(from=ifelse(is.null(from), -Inf, from), to=ifelse(is.null(to), Inf, to))
    x <- x[,start(x)>=ranges["from"] & end(x)<=ranges["to"]]
    if(sort)
        x <- x[,order(range(x))]
    return(x)
})

## ReferenceDataTracks need to stream the data from file and then pass the results on to the next method
setMethod("subset", signature(x="ReferenceDataTrack"), function(x, from, to, chromosome, ...){
    ## We only need to reach out into the referenced file once if the range is already contained in the object
    if(missing(from) || is.null(from) || missing(to) || is.null(to))
        stop("Need both start and end location to subset a ReferenceDataTrack")
    if(missing(chromosome) || is.null(chromosome))
        chromosome <- Gviz::chromosome(x)
    subRegion <- GRanges(seqnames=chromosome[1], ranges=IRanges(start=from, end=to))
    if(length(ranges(x))==0 || !all(overlapsAny(ranges(x),subRegion))){
        vals <- x@stream(x@reference, subRegion)
        x@range <- vals
        mcols(x@range) <- NULL
        x@data <- .prepareDtData(if(ncol(values(vals))) as.data.frame(values(vals)) else matrix(nrow=0, ncol=0), length(vals))
        chromosome(x) <- chromosome[1]
    }
    return(callNextMethod(x=x, from=from, to=to, drop=FALSE, ...))
})

## Only recompute the stacks here
setMethod("subset", signature(x="StackedTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE, ...){
    x <- callNextMethod(x=x, from=from, to=to, sort=sort)
    if(stacks)
        x <- setStacks(x)
    return(x)
})

## In order to keep the grouping information for track regions in the clipped areas we have to
## keep all group elements that overlap with the range
setMethod("subset", signature(x="AnnotationTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE, use.defaults=TRUE, ...){
    ## Subset to a single chromosome first
    csel <- seqnames(x) != chromosome(x)
    if(any(csel))
        x <- x[!csel]
    if(length(x))
    {
        ## Nothing to do if everything is within the range
        ranges <- if(use.defaults) .defaultRange(x, from=from, to=to) else c(from=ifelse(is.null(from), -Inf, from), to=ifelse(is.null(to), Inf, to))
        if(!(any(end(x)<ranges["from"] | start(x)> ranges["to"]))){
            if(stacks)
                x <- setStacks(x)
            return(x)
        }
        ## Now remove everything except for the overlapping groups by first subselecting all groups in the range...
        granges <- unlist(range(split(ranges(x), group(x))))
        gsel <- names(granges)[subjectHits(findOverlaps(GRanges(seqnames=chromosome(x), ranges=IRanges(ranges["from"], ranges["to"])), granges))]
        x <- x[group(x) %in% gsel]
        if(sort)
            x <- x[order(range(x)),]
        if(stacks)
            x <- setStacks(x)
    }
    return(x)
})

## ReferenceDataTracks need to stream the data from file and then pass the results on to the next method
setMethod("subset", signature(x="ReferenceAnnotationTrack"), function(x, from, to, chromosome, ...){
    ## We only need to reach out into the referenced file once if the range is already contained in the object
    if(missing(from) || is.null(from) || missing(to) || is.null(to))
        stop("Need both start and end location to subset a ReferenceAnnotationTrack")
    if(missing(chromosome) || is.null(chromosome))
        chromosome <- Gviz::chromosome(x)
    subRegion <- GRanges(seqnames=chromosome[1], ranges=IRanges(start=from, end=to))
    if(length(ranges(x))==0 || all(overlapsAny(ranges(x),subRegion))){
        cMap <- .resolveColMapping(x@stream(x@reference, subRegion), x@args, x@mapping)
        x@range <- .buildRange(cMap$data, args=cMap$args, defaults=x@defaults, trackType="AnnotationTrack")
        chromosome(x) <- chromosome[1]
    }
    return(callNextMethod(x=x, from=from, to=to, drop=FALSE, ...))
})

## FIXME: Still needs to be implemented
setMethod("subset", signature(x="ReferenceGeneRegionTrack"), function(x, ...){
    warning("ReferenceGeneRegionTrack objects are not supported yet.")
    return(callNextMethod())
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
setMethod("subset", signature(x="AlignedReadTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE, ...){
    if(x@coverageOnly)
    {
        if(is.null(from))
            from <- min(unlist(lapply(x@coverage, function(y) if(length(y)) min(start(y)))))
        if(is.null(to))
            to <- max(unlist(lapply(x@coverage, function(y) if(length(y)) max(start(y)))))
        x@coverage <- lapply(x@coverage, function(y){runValue(y)[end(y)<from | start(y)>to] <- 0; y})
        x@coverage <- lapply(x@coverage, function(y){ if (length(y) < to) y <- c(y, Rle(0, to-length(y))); y})
        ## 
        ##from <- min(unlist(lapply(x@coverage, function(y) if (length(y)) head(start(y), 2)[2])))
        if (max(unlist(lapply(x@coverage, function(y) {length(runLength(y)[runValue(y)!=0])}))))
        {
        from <- min(unlist(lapply(x@coverage, function(y) if(length(y)) head(start(y)[runValue(y)!=0],1))))
        to <- max(unlist(lapply(x@coverage, function(y) if(length(y)) tail(end(y),2)[1])))
        }
        x@range <- GRanges(range=IRanges(start=from, end=to), strand=names(x@coverage), seqnames=x@chromosome)
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
    if(as.logical(.dpOrDefault(GdObject, "legend", FALSE)) && !is.null(getPar(GdObject, ".__groupLevels"))){
        pushViewport(viewport(y=1, height=unit(1, "npc") - unit(getPar(GdObject, ".__verticalSpace"), "inches"),
                              just=c(0.5, 1)))
        on.exit(popViewport(1))
    }
    if(.dpOrDefault(GdObject, "showAxis", TRUE)) {
	callNextMethod()
    } else {
	return(NULL)
    }
})

setMethod("drawAxis", signature(GdObject="NumericTrack"), function(GdObject, from, to, ...) {
    type <- match.arg(.dpOrDefault(GdObject, "type", "p"), c("p", "l", "b", "a", "s", "g", "r", "S", "smooth", "polygon",
                                                             "histogram", "mountain", "h", "boxplot", "gradient", "heatmap"),
                      several.ok=TRUE)
    yvals <- values(GdObject)
    ylim <- .dpOrDefault(GdObject, "ylim", if(!is.null(yvals) && length(yvals)) 
                         range(yvals, na.rm=TRUE, finite=TRUE) else c(-1,1))
    if(diff(ylim)==0)
        ylim <- ylim+c(-1,1) 
    hSpaceAvail <- vpLocation()$isize["width"]/6
    yscale <- extendrange(r=ylim, f=0.05)
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
    nlevs <- max(1, nlevels(factor(getPar(GdObject, "groups"))))
    if(type=="heatmap" && .dpOrDefault(GdObject, "showSampleNames", FALSE)){
        groups <- .dpOrDefault(GdObject, "groups")
        sn <- if(is.null(groups)) rownames(values(GdObject)) else rev(unlist(split(rownames(values(GdObject)), factor(groups))))
        cex.sn <- .dpOrDefault(GdObject, "cex.sampleNames", acex)
        col.cn <- .dpOrDefault(GdObject, "col.sampleNames", "white")
        wd <- max(as.numeric(convertWidth(stringWidth(sn) + unit(10, "points"), "npc"))) * cex.sn
        samNames <- viewport(x=1, width=wd, just=1, yscale=c(-0.05, 1.05))
        pushViewport(samNames)
        nr <- nrow(values(GdObject))
        yy <- head(seq(0.05, 0.95, len=nr+1), -1)
        yy <- yy + diff(yy)[[1]]/2
        grid.text(x=rep(0.5, nr), y=yy, label=rev(sn), just=0.5, gp=gpar(cex=cex.sn, col=col.cn))
        popViewport(1)
        samAxis <- viewport(x=1-wd, width=1-wd, just=1)
        pushViewport(samAxis)
        on.exit(popViewport(1))
    }
    ## if any of the types are gradient or heatmap we want the gradient scale
    if(any(type %in% c("gradient", "heatmap")) && .dpOrDefault(GdObject, "showColorBar", TRUE)){
        ## viewport to hold the color strip
        shift <- ifelse(all(type %in% c("gradient", "heatmap")), 1, 0)
        pcols <- .getPlottingFeatures(GdObject)
        ncolor <- .dpOrDefault(GdObject, "ncolor", 100)
        vpAxisCont <- viewport(x=unit(1, "npc")-unit(2-shift, "points"), width=unit(1, "npc")-unit(2-shift, "points"), just=1)
        pushViewport(vpAxisCont)
        for(i in seq_len(nlevs)){
            ## create color palette
            cr <- c("white", pcols$col[i])
	    if(nlevs<2)
		cr <- .dpOrDefault(GdObject, "gradient", cr)
            palette <- colorRampPalette(cr)(ncolor+5)[-(1:5)]
            pshift <- ifelse(i==nlevs, 1-shift, 0)
            vpTitleAxis <- viewport(x=unit(1, "npc")-unit(4*(i-1), "points"), width=unit(4+pshift, "points"),
                                    yscale=yscale, just=1)
            pushViewport(vpTitleAxis)
            ## draw a rectangle for each color
            if(all(type %in% c("gradient", "heatmap"))){
                if(i==nlevs)
                    suppressWarnings(grid.yaxis(gp=gpar(col=acol, cex=acex), at=at))
                grid.rect(y=unit(seq(ylim[1],ylim[2],length.out=ncolor+1),"native")[-(ncolor+1)], x=unit(0, "npc")-unit(1, "points"),
                          width=1, height=1/ncolor, gp=gpar(fill=palette, lty=0), just=c("left","bottom"))
            }else{
                grid.rect(y=unit(seq(ylim[1],ylim[2],length.out=ncolor+1),"native")[-(ncolor+1)], x=0,
                          width=1, height=1/ncolor, gp=gpar(fill=palette, lty=0), just=c("left","bottom"))
                if(i==nlevs){
                    suppressWarnings(grid.yaxis(gp=gpar(col=acol, cex=acex), at=at))
                    grid.lines(x=c(0,0), y=ylim, gp=gpar(col=acol), default.units="native")
                }
            }
            popViewport(1)
        }
	popViewport(1)
    } else {
	vpTitleAxis <- viewport(x=0.95, width=0.2, yscale=yscale, just=0)
	pushViewport(vpTitleAxis)
	suppressWarnings(grid.yaxis(gp=gpar(col=acol, cex=acex), at=at))
	grid.lines(x=c(0,0), y=ylim, gp=gpar(col=acol), default.units="native")
	popViewport(1)
    }
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
                            ## We have to figure out the data range, taking transformation into account
                            ylim <- .dpOrDefault(GdObject, "ylim")
                            if (is.null(ylim)) {
                                maxs <- sapply(c("+", "-"), function(s) {
                                    cvr <- coverage(GdObject, strand=s)
                                    if (length(cvr)) max(cvr, na.rm=TRUE, finite=TRUE) else 0L
                                })
                                y.max <- max(maxs, na.rm=TRUE, finite=TRUE)
                                ylim <- c(0, if (y.max == 0) 1 else y.max)
                                trans <- displayPars(GdObject, "transformation")[[1]]
                                if (!is.null(trans))
                                    ylim <- c(0, trans(ylim[2]))
                            }
                            for(s in c("+", "-")) {
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
##    In plotting mode: the object is plotted. Return value is the object with optional HTML image map information
##    added to the imageMap slot
## Since subsetting can be potentially expensive when the data are large we want to minimize this operation. Essentially it
## should be done only once before any other plotting or computation starts, hence we expect the GdObject in the drawGD
## methods to already be trimmed to the correct size
##----------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------
## The default method for all StackedTrack types which always should be called (this has to be done explicitely using
## callNextMethod)
##----------------------------------------------------------------------------------------------------------------------------
## Although the stacking type is not stored as a displayParameter we still want to check whether it is
## included there and set the actual stacking of the object accordingly
setMethod("drawGD", signature("StackedTrack"), function(GdObject, ...){
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if(debug || debug=="prepare")
        browser()
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

## Compute the coordinates and colors for all track items (e.g. exons in a GeneRegionTrack)
.boxes <- function(GdObject, offsets)
{
    ylim <- c(0, 1)
    h <- diff(ylim)
    middle <- mean(ylim)
    sh <- max(0, min(h, .dpOrDefault(GdObject, "stackHeight", 0.75)))
    space <- (h-(h*sh))/2
    if (inherits(GdObject, "GeneRegionTrack")) {
        thinBox <- .dpOrDefault(GdObject, "thinBoxFeature", c("utr", "ncRNA", "utr3", "utr5", "miRNA", "lincRNA"))
        space <- ifelse(feature(GdObject) %in% thinBox, space + ((middle -
                                                                  space) / 2), space)
    }
    shape <- .dpOrDefault(GdObject, "shape", "arrow")  
    color <- .getBiotypeColor(GdObject)
    id <- identifier(GdObject, lowest=TRUE)
    sel <- grepl("\\[Cluster_[0-9]*\\]", id)
    id[sel] <- sprintf("%i merged\n%s", as.integer(.getAnn(GdObject, "density")[sel]),
                       ifelse(class(GdObject) %in% c("AnnotationTrack", "DetailsAnnotationTrack"), "features", "exons"))
    boxes <- data.frame(cx1=start(GdObject), cy1=ylim[1]+space+offsets, cx2=start(GdObject)+width(GdObject), cy2=ylim[2]-space+offsets,
                        fill=color, strand=strand(GdObject), text=id, textX=start(GdObject)+(width(GdObject)/2), textY=middle+offsets,
                        .getImageMap(cbind(start(GdObject), ylim[1]+space+offsets, end(GdObject), ylim[2]-space+offsets)),
                        start=start(GdObject), end=end(GdObject), values(GdObject), stringsAsFactors=FALSE)
    rownames(boxes) <- if(is(GdObject, "GeneRegionTrack") && .dpOrDefault(GdObject, "collapseTranscripts", FALSE))
            sprintf("uid%i", seq_along(identifier(GdObject))) else make.unique(identifier(GdObject, lowest=TRUE))
    return(boxes)
}

## Compute the coordinates for the bars connecting grouped items and the group labels
.barsAndLabels <- function(GdObject)
{
    bins <- stacks(GdObject)
    stacks <- max(bins)
    res <- .pxResolution(coord="x")
    gp <- group(GdObject)
    grpSplit <- split(range(GdObject), gp)
    grpRanges <- unlist(range(grpSplit))
    needBar <- sapply(grpSplit, length)>1 & width(grpRanges) > res
    ## If we draw the bar from start to end of the range we sometimes see little overlaps that extend beyond the first or last item.
    ## In order to fix this, we just substract the equivalent of min.width pixels from both ends of each group range
    min.swidth <- res*.dpOrDefault(GdObject, "min.width", 2)
    nstart <- start(grpRanges[needBar])+min.swidth
    nend <- end(grpRanges[needBar])-min.swidth
    sel <- (nend-nstart)>0
    start(grpRanges[needBar][sel]) <- nstart[sel]
    end(grpRanges[needBar][sel]) <- nend[sel]
    strand <- sapply(split(strand(GdObject), gp), function(x){
        tmp <- unique(x)
        if(length(tmp)>1) "*" else tmp
    })
    yloc <- sapply(split((stacks-bins)+1, gp), function(x) unique(x))+0.5
    color <- if(length(grep("__desatCol", values(GdObject)$feature[1])))
        .dpOrDefault(GdObject, "fill", .DEFAULT_FILL_COL) else sapply(split(.getBiotypeColor(GdObject), gp), head, 1)
    bars <- data.frame(sx1=start(grpRanges)[needBar], sx2=end(grpRanges)[needBar], y=yloc[needBar], strand=strand[needBar],
                      col=color[needBar], stringsAsFactors=FALSE)
    labs <- sapply(split(identifier(GdObject, add.space=TRUE), gp), head, 1)
    lsel <- grepl("\\[Cluster_[0-9]*\\]", labs)
    if(any(lsel)){
        gdens <- as.integer(sapply(split(.getAnn(GdObject, "gdensity"), gp), head, 1))
        labs[lsel] <- sprintf("%i merged %s  ", gdens[lsel],
                              ifelse(class(GdObject) %in% c("AnnotationTrack", "DetailsAnnotationTrack"), "groups", "gene models"))
    }
    offs <- rep(min.swidth, length(grpRanges))
    offs[sapply(grpSplit, length)<=1] <- 0
    labels <- data.frame(txt=labs, x=start(grpRanges)-offs, y=yloc, stringsAsFactors=FALSE)
    return(list(bars=bars, labels=labels))
}

## The actual drawing method
setMethod("drawGD", signature("AnnotationTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...){
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if(debug || debug=="prepare")
        browser()
    imageMap(GdObject) <- NULL
    if(!length(GdObject))
        return(invisible(GdObject))
    ## In prepare mode we need to make sure that the stacking information is updated from the optional display parameter (by calling
    ## the StackedTrack drawGD method) and also perform the collapsing of track items which could potentially lead to re-stacking.
    if(prepare){
        GdObject <- callNextMethod()
        bins <- stacks(GdObject)
        stacks <- max(bins)
        ## We need to collapse the track object based on the current screen resolution (note that this may trigger re-stacking)
        pushViewport(dataViewport(xData=c(minBase, maxBase), extension=0, yscale=c(1, stacks+1), clip=TRUE))
        GdObject <- collapseTrack(GdObject, diff=.pxResolution(coord="x"), xrange=c(minBase, maxBase))
        popViewport(1)
        return(invisible(GdObject))
    }
    if(debug || debug=="draw")
        browser()
    ## If there are too many stacks for the available device resolution we cast an error
    bins <- stacks(GdObject)
    stacks <- max(bins)
    yscale <- if(!.dpOrDefault(GdObject, "reverseStacking", FALSE)) c(1, stacks+1) else c(stacks+1, 1)
    pushViewport(dataViewport(xData=c(minBase, maxBase), extension=0, yscale=yscale, clip=TRUE))
    res <- .pxResolution(coord="x")
    curVp <- vpLocation()
    if(curVp$size["height"]/stacks < .dpOrDefault(GdObject, "min.height", 3))
        stop("Too many stacks to draw. Either increase the device size or limit the drawing to a smaller region.")
    ## We adjust the color saturation to indicate overplotting if necessary
    if(.dpOrDefault(GdObject, "showOverplotting", FALSE))
    {
        dens <- as.numeric(values(GdObject)$density)
        if(length(unique(dens))!=1)
        {
            minSat <- max(0.25, 1/max(dens))
            minDens <- min(dens)
            rDens <- diff(range(dens))
            saturation <- minSat+((dens-minDens)/rDens/(1/(1-minSat)))
            bc <- unique(.getBiotypeColor(GdObject))
            baseCol <- rgb2hsv(col2rgb(bc))
            desatCols <- unlist(lapply(saturation, function(x) hsv(baseCol[1,], x, baseCol[3,])))
            names(desatCols) <- paste(unique(feature(GdObject)), rep(dens, each=length(bc)), sep="_")
            feature(GdObject) <- paste(feature(GdObject), dens, sep="_")
            desatCols <- desatCols[unique(names(desatCols))]
            displayPars(GdObject) <- as.list(desatCols)
        }
    }
    ## Now we can pre-compute all the coordinates and settings for the elements to be drawn...
    box <- .boxes(GdObject, (stacks-bins)+1)
    barsAndLab <- .barsAndLabels(GdObject)
    bar <- barsAndLab$bars
    bartext <- barsAndLab$labels
    ## ... and then draw whatever is needed
    shape <- .dpOrDefault(GdObject, "shape", "arrow")
    border <- .dpOrDefault(GdObject, "col")[1]
    col.line <- .dpOrDefault(GdObject, "col.line")[1]
    if(is.null(border))
        border <- ifelse(is(GdObject, "GeneRegionTrack"), NA, "transparent")
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
    fontcolor.group <- .dpOrDefault(GdObject, "fontcolor.group", .DEFAULT_SHADED_COL)[1]
    cex.group <- .dpOrDefault(GdObject, "cex", 1) * .dpOrDefault(GdObject, "cex.group", 0.6)
    fontsize.group <- .dpOrDefault(GdObject, "fontsize.group", fontsize)
    fontface.group  <- .dpOrDefault(GdObject, "fontface.group", fontface)
    fontfamily.group <- .dpOrDefault(GdObject, "fontfamily.group", fontfamily)
    if(nrow(box)>0){
        if(nrow(bar)>0)
            .arrowBar(bar$sx1, bar$sx2, y=bar$y, bar$strand, box[,1:4, drop=FALSE],
                      col=if(is.null(col.line)) bar$col else rep(col.line, length(bar$col)), lwd=lwd, lty=lty,
                      alpha=alpha, barOnly=(!"smallArrow" %in% .dpOrDefault(GdObject, "shape", "box") || stacking(GdObject)=="dense"),
                      diff=res, min.height=.dpOrDefault(GdObject, "min.height", 3))
        if("box" %in% shape || ("smallArrow" %in% shape && !"arrow" %in% shape))
            grid.rect(box$cx2, box$cy1, width=box$cx2-box$cx1, height=box$cy2-box$cy1,
                      gp=gpar(col=if(is.na(border)) box$fill else border, fill=box$fill, lwd=lwd, lty=lty, alpha=alpha),
                      default.units="native", just=c("right", "bottom"))
       
        if("ellipse" %in% shape){
            ellCoords <- .box2Ellipse(box)
            grid.polygon(x=ellCoords$x1, y=ellCoords$y1, id=ellCoords$id,
                         gp=gpar(col=if(is.na(border)) box$fill else border, fill=box$fill, lwd=lwd, lty=lty, alpha=alpha),
                         default.units="native")
        }
        if("arrow" %in% shape && !"box" %in% shape){
            .filledArrow(box[,1:4], col=border, fill=box$fill, lwd=lwd, lty=lty, alpha=alpha, strand=box$strand, min.width=6*res)
        }
        if(.dpOrDefault(GdObject, "showFeatureId", FALSE))
            grid.text(box$text, box$textX, box$textY, rot=rotation,
                      gp=gpar(col=fontcolor, cex=cex, fontsize=fontsize, fontface=fontface, lineheight=lineheight,
                              fontfamily=fontfamily), default.units="native", just=c("center", "center"))
      
        if(.dpOrDefault(GdObject, "showId", FALSE) && nrow(bartext)>0 && stacking(GdObject)!="dense")
            grid.text(bartext$txt, bartext$x, bartext$y, gp=gpar(col=fontcolor.group, cex=cex.group, fontsize=fontsize.group,
                                                                 fontface=fontface.group, fontfamily=fontfamily.group),
                      default.units="native", just=c("right", "center"))
    }
    popViewport(1)
    ## Finaly we set up the image map
    ## FIXME: we may want to record the merging information here
    im <- if(!is.null(box)) {
        coords <- as.matrix(box[,c("x1", "y1", "x2", "y2"),drop=FALSE])
        restCols <- setdiff(colnames(box), c("x1", "x2", "y1", "y2", "cx1", "cx2", "cy1", "cy2", "textX", "textY"))
        tags <- sapply(restCols, function(x){
            tmp <- as.character(box[,x])
            names(tmp) <- rownames(coords)
            tmp}, simplify=FALSE)
        tags$title <- identifier(GdObject)
        ImageMap(coords=coords, tags=tags) } else NULL
    imageMap(GdObject) <- im
    return(invisible(GdObject))
})

## For a GeneRegionTrack we just set the showExonId alias and then call the AnnotationTrack method
setMethod("drawGD", signature("GeneRegionTrack"), function(GdObject,  ...){
    displayPars(GdObject) <- list(showFeatureId=as.vector(displayPars(GdObject, "showExonId")))
    GdObject <- callNextMethod(GdObject, ...)
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------
## Draw a genome axis
##----------------------------------------------------------------------------------------------------------------------------
.expLabel <- function(GdObject, tckText, prune=FALSE){
    tck <- tckText
    exponent <- if(is.null(.dpOrDefault(GdObject, "exponent", NULL))){
        exp <- 0
        while(all(tck[tck>0]/10^exp >= 1))
            exp <- exp+3
        exp-3
    } else  max(0, getPar(GdObject, "exponent"))
    if(exponent > 0){
        tckText <- tckText/(10^exponent)
    }
    if(prune){
        tmp <- as.character(tckText)
        count <- max(nchar(gsub("*.\\.", "", tmp)))
        while(count>1 && !any(duplicated(round(tckText, count)))){
            count <- count-1
        }
        tckText <- round(tckText, count+1)
    }
    return(switch(as.character(exponent),
                  "0"=sprintf("%i", as.integer(tckText)),
                  "3"=sprintf("%s kb", tckText),
                  "6"=sprintf("%s mb", tckText),
                  "9"=sprintf("%s gb", tckText),
                  sapply(tckText, function(x) bquote(paste(.(x), " ",10^.(exponent))))))
 }

setMethod("drawGD", signature("GenomeAxisTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if((is.logical(debug) && debug) || debug=="prepare")
        browser()
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
    ids <- as.character(values(GdObject)$id)
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
        nsp <- if(is.null(.dpOrDefault(GdObject, "scale", NULL))){
            (sum(tickHeight, pyOff*2, textYOff*2 + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex)*2*1.3)/pres["y"]
        } else {
            labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
            if(labelPos %in% c("above", "below")){
                (sum(tickHeight, pyOff*2 + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex)*2)/pres["y"]
            } else {
                (sum(tickHeight, pyOff*2 + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex))/pres["y"]
            }
        }
        displayPars(GdObject) <- list("neededVerticalSpace"=nsp)
        popViewport(1)
        return(invisible(GdObject))
    }
    if((is.logical(debug) && debug) || debug=="draw")
        browser()
    ## Plot range if there is any
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
	
    ## in "scale" mode we just plot a simple scale and return ...
    scaleLen <- .dpOrDefault(GdObject, "scale", NULL)
    if(!is.null(scaleLen))
    {
        len <- (maxBase-minBase + 1)
        if(scaleLen>len)
		{
            warning(paste("scale (", scaleLen,
                          ") cannot be larger than plotted region",
                          len, " - setting to ~5%\n", sep=""))
            scaleLen = 0.05
        }
        xoff <- len * 0.03 + minBase
        labelPos <- match.arg(labelPos, c("alternating", "revAlternating", "above", "below", "beside"))
		if(scaleLen<=1 && scaleLen>0) { # calculate and round the scale 
			scaleLen <- len * scaleLen
			ex <- .dpOrDefault(GdObject, "exponent", floor(log10(scaleLen)))
			v <- round(scaleLen, -ex)
			if(v==0)  v <- scaleLen
		} else { # if the scale is an absolute value don't round
			ex <- .dpOrDefault(GdObject, "exponent", floor(log10(scaleLen)))
			v <- scaleLen
		}
	
        ## work out exponent/unit
        label <- .expLabel(GdObject, v)
        grid.lines(x=c(xoff, v+xoff), y=c(0,0), default.units="native", gp=gpar(col=color, lwd=lwd, alpha=alpha))
        grid.segments(x0=c(xoff, v+xoff), y0=c(0-tickHeight, 0-tickHeight),
                      x1=c(xoff, v+xoff), y1=c(tickHeight,tickHeight, tickHeight),
                      default.units="native", gp=gpar(col=color, lwd=lwd, alpha=alpha))
        z <- len * 0.01
        if(labelPos=="below"){
            grid.text(label=if(is.character(label)) label else label[[1]], x=xoff+v/2, y=0-(tickHeight/1.5*dfact),
                      just=c("center", "top"), gp=gpar(alpha=alpha, col=color, cex=cex, fontface=fontface), default.units="native")
        } else if(labelPos=="above"){
            grid.text(label=if(is.character(label)) label else label[[1]], x=xoff+v/2, y=tickHeight/1.5*dfact, just=c("center", "bottom"),
                      gp=gpar(alpha=alpha, col=color, cex=cex, fontface=fontface), default.units="native")
        } else {
            grid.text(label=if(is.character(label)) label else label[[1]], x=v+xoff+z, y=0, just=c("left", "center"),
                      gp=gpar(alpha=alpha, col=color, cex=cex, fontface=fontface), default.units="native")
        }
        popViewport(1)
        return(invisible(GdObject))
    }

    GdObject <- GdObject[end(GdObject) > axRange[1] & start(GdObject) < axRange[2]]
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
        if(is.null(ids) || length(ids)==0)
            ids <- as.character(seq_len(nrow(map)))
        rownames(map) <- make.unique(as.character(ids))
        tags <- lapply(list(title=ids, start=as.character(start(GdObject)), end=as.character(end(GdObject))),
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
    label <- .expLabel(GdObject, tck)
    ylabs <- y1t + (ifelse(y1t>0, 1, -1) * (textYOff + (as.numeric(convertHeight(stringHeight("1"),"native"))/2)*cex))
    ttck <- if(min(diff(tck))==1) tck+0.5 else tck
    if(is.character(label)){
        grid.text(label=label, x=ttck, y=ylabs, just=c("centre", "centre"),
                  gp=gpar(cex=cex, fontface=fontface), default.units="native")
    }else{
        for(i in seq_along(label))
            grid.text(label=label[[i]], x=ttck[i], y=ylabs[i], just=c("centre", "centre"),
                      gp=gpar(cex=cex, fontface=fontface), default.units="native")
    }
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
            llabel <- .expLabel(GdObject, ltck[sel], prune=TRUE)
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
                  just=c("left", "top"), gp=gpar(cex=cex*0.75, fontface=fontface),
                  default.units="native")
    }
    popViewport()
    return(invisible(GdObject))})
##----------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------
## Draw DetailsAnnotationTrack
##----------------------------------------------------------------------------------------------------------------------------
## Create a data.frame with the distinct details function arguments (like start, end, ...)
.buildArgsDf <- function(GdObject)
{
      groupDetails <- .dpOrDefault(GdObject, "groupDetails", FALSE)
      rr <- if(groupDetails) unlist(range(split(ranges(GdObject), group(GdObject)))) else ranges(GdObject)
      args <- data.frame(start=as.integer(start(rr)), end=as.integer(end(rr)), strand=as.character(strand(rr)),
                         chromosome=as.character(seqnames(rr)),
                         identifier=as.character(if(groupDetails) names(rr) else identifier(GdObject, lowest=TRUE)),
                         stringsAsFactors=FALSE)
      return(args)
}


setMethod("drawGD", signature("DetailsAnnotationTrack"),
          function(GdObject,  minBase, maxBase, prepare=FALSE, ...){
              debug <- .dpOrDefault(GdObject, "debug", FALSE)
              if((is.logical(debug) && debug) || debug=="prepare")
                  browser()
              adf <- .buildArgsDf(GdObject)
              args <- .dpOrDefault(GdObject, "detailsFunArgs", fromPrototype=TRUE)
              groupDetails <- .dpOrDefault(GdObject, "groupDetails", FALSE)
              if(prepare){
                  GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
                  GdObject <- GdObject[order(start(GdObject))]
                  indices <- if(groupDetails) seq_len(length(unique(group(GdObject)))) else seq_len(length(GdObject))
                  pushViewport(viewport(xscale=c(minBase, maxBase)))
                  hasWarned <- FALSE
                  select <- sapply(indices, function(i){
                      iargs <- as.list(adf[i,])
                      iargs$index <- i
                      iargs$GdObject <- GdObject
                      iargs$GdObject.original <- .dpOrDefault(GdObject, ".__OriginalGdObject", GdObject)
                      args <- c(args[setdiff(names(args), names(iargs))], iargs)
                      res <- do.call(GdObject@selectFun, args)
                      if(length(res)!=1 || !is.logical(res) || is.na(res)){
                          if(!hasWarned)
                              warning("The result of function 'selectFun' has to be a single logical value. Forcing the value to 'TRUE'")
                          hasWarned <<- TRUE
                          res <- TRUE
                      }
                      res
                  })
                  popViewport(1)
                  displayPars(GdObject) <- list(".__select"=select)
                  return(invisible(GdObject))
              }
              if((is.logical(debug) && debug) || debug=="draw")
                  browser()
              n <- length(GdObject)
              col <- rep(.dpOrDefault(GdObject, "detailsConnector.col", fromPrototype=TRUE), n)[1:n]
              lty <- rep(.dpOrDefault(GdObject, "detailsConnector.lty", fromPrototype=TRUE), n)[1:n]
              lwd <- rep(.dpOrDefault(GdObject, "detailsConnector.lwd", fromPrototype=TRUE), n)[1:n]
              pch <- rep(.dpOrDefault(GdObject, "detailsConnector.pch", fromPrototype=TRUE), n)[1:n]
              cex <- rep(.dpOrDefault(GdObject, "detailsConnector.cex", fromPrototype=TRUE), n)[1:n]
              border.lty <- rep(.dpOrDefault(GdObject, "detailsBorder.lty", fromPrototype=TRUE), n)[1:n]
              border.lwd <- rep(.dpOrDefault(GdObject, "detailsBorder.lwd", fromPrototype=TRUE), n)[1:n]
              border.col <- rep(.dpOrDefault(GdObject, "detailsBorder.col", fromPrototype=TRUE), n)[1:n]
              border.fill <-rep(.dpOrDefault(GdObject, "detailsBorder.fill", fromPrototype=TRUE), n)[1:n]
              minwidth <- .dpOrDefault(GdObject, "details.minWidth", fromPrototype=TRUE)
              size <- .dpOrDefault(GdObject, "details.size", fromPrototype=TRUE)
              xyratio <- .dpOrDefault(GdObject, "details.ratio", fromPrototype=TRUE)
              if ( 0 >= size || size > 1 ) {
                  warning("details.size must be >0 and <1 - reset to 0.5")
                  size = 0.5
              }
              selection <- .dpOrDefault(GdObject, ".__select", rep(TRUE, length(GdObject)))
              len <- sum(selection)
              bins <- if(!groupDetails) stacks(GdObject) else sapply(split(stacks(GdObject) , group(GdObject)), unique)
              stacks <- max(bins)
              if(len>0){
                  if( ((maxBase-minBase)/len)/.pxResolution(coord="x") < minwidth ) {
                      warning("too much detail for available space (plot fewer annotation or increase details.minWidth)!")
                      popViewport(1)
                      GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
                      return(GdObject)
                  }
                  rr <- if(groupDetails) unlist(range(split(ranges(GdObject), group(GdObject)))) else ranges(GdObject)
                  xloc1 <- (end(rr) - start(rr))/2+start(rr)
                  yloc1 <- (stacks - (bins - 0.5)+1)
                  xloc2 <- ((1/len*seq_len(len))-1/len + (1/len*0.5))
                  yloc2 <- rep(1, len) 
                  ## draw details plots (via user supplied function 'fun')
                  pushViewport(viewport(height=size, y=1-size, just=c(0.5, 0)))
                  w <- 1
                  v <- 0
                  vpl <- vpLocation()
                  r <- vpl$size["width"]/len/vpl$size["height"]
                  if ( r > xyratio ) {
                      w <- xyratio/r
                      v <- ((1/len) - (1/len*w))/2
                  }
                  indices <- if(groupDetails) seq_len(length(unique(group(GdObject)))) else seq_len(length(GdObject))
                  j <- 1
                  pres <- list()
                  hasError <- FALSE
                  for(i in indices[selection]) {
                      pushViewport(viewport(width=1/len*w, x=((1/len*j)-1/len)+(v), just=c(0, 0.5)))
                      grid.rect(gp=gpar(col=border.col[i], lwd=border.lwd[i], lty=border.lty[i], fill=border.fill[i]))
                      iargs <- as.list(adf[i,])
                      iargs$index <- i
                      iargs$GdObject <- GdObject
                      iargs$GdObject.original <- .dpOrDefault(GdObject, ".__OriginalGdObject", GdObject)
                      args <- c(args[setdiff(names(args), names(iargs))], iargs)
                      pres[[as.character(j)]] <- try(do.call(GdObject@fun, args), silent=TRUE)
                      if(!is.null(pres) && is(pres[[as.character(j)]], "try-error")){
                          hasError <- TRUE
                          grid.segments(x0=c(0.1,0.1), x1=c(0.9,0.9), y0=c(0.9,0.1), y1=c(0.1,0.9), gp=gpar(col="red", lwd=3))
                      }
                      popViewport(1)
                      j <- j+1
                  }
                  if(hasError)
                      warning("There have been errors in the detail plotting function:\n", paste(pres, collapse="\n"))
                  popViewport(1)
                  ## plot AnnotationTrack and connectors to details
                  pushViewport(viewport(xscale=c(minBase, maxBase),
                                        yscale=c(1, stacks+1), clip=FALSE,
                                        height=1-size, y=0, just=c(.5, 0)))
                  GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...)
                  grid.segments(x0=unit(xloc1[selection], "native"), x1=xloc2, y0=unit(yloc1[selection], "native"),
                                y1=yloc2, gp=gpar(col=col, lwd=lwd, lty=lty, cex=cex))
                  grid.points(x=unit(xloc2, "npc"), y=unit(yloc2, "npc"), gp=gpar(col=col, cex=cex), pch=pch)
                  grid.points(x=unit(xloc1[selection], "native"), y=unit(yloc1[selection], "native"), gp=gpar(col=col, cex=cex), pch=pch)
                  popViewport(1)
              }else{
                  GdObject <- callNextMethod(GdObject,  minBase, maxBase, prepare=prepare, ...) 
              }
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
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if((is.logical(debug) && debug) || debug=="prepare")
        browser()
    imageMap(GdObject) <- NULL
    type <- .dpOrDefault(GdObject, "type", "p")
    type <- match.arg(type, c("p", "l", "b", "a", "s", "g", "r", "S", "smooth",
                              "histogram", "mountain", "h", "boxplot", "gradient", "heatmap", "polygon"),
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
                valsS <- if(ncol(vals)) matrix(sapply(split(vals, groups),
                                                      function(x) agFun(t(matrix(x, ncol=ncol(vals))))), nrow=ncol(vals)) else{
                    matrix(nrow=nlevels(groups), ncol=0, dimnames=list(levels(groups)))}
                displayPars(GdObject) <- list(".__valsS"=valsS)
                if(stacked==TRUE && is.null(ylim))
                {
                    ylim <- suppressWarnings(range(unlist(apply(valsS, 1, function(x){
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
                        tmp}))))
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
        if(!is.factor(grps))
            grps <- factor(grps)
        if(is.null(grps) || length(grps)==1 || length(setdiff(type, c("gradient", "mountain", "grid"))) == 0)
            displayPars(GdObject) <- list(legend=FALSE)
        if(as.logical(as.logical(.dpOrDefault(GdObject, "legend", FALSE))) && nlevels(grps)>1){
            cex <- .dpOrDefault(GdObject, "cex.legend", 0.8)
            fontsize <- .dpOrDefault(GdObject, "fontsize.legend", 12)
            fontface <- .dpOrDefault(GdObject, "fontface.legend", 1)
            lineheight <- .dpOrDefault(GdObject, "lineheight.legend", 1)
            fontfamily <- .dpOrDefault(GdObject, "fontfamily.legend", 1)
            pushViewport(viewport(width=unit(1, "npc")-unit(0.2,"inches"),
                                  gp=gpar(cex=cex, fontsize=fontsize, fontface=fontface,
                                          lineheight=lineheight)))
            grps <- levels(grps)
            legInfo <- .legendInfo()[type,, drop=FALSE]
            for(i in colnames(legInfo))
                legInfo[,i] <- any(legInfo[,i]) && !any(duplicated(pcols[[i]][1:length(grps)]))
            legFactors <- sort(names(which(apply(legInfo, 2, any))))
            boxSize <-  if(length(setdiff(legFactors, c("col", "cex")))==0) 0.1 else 0.3
            spacing <- 0.1
            hspacing <- 0.02
            lengths <- as.numeric(convertUnit(stringWidth(grps),"inches"))
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
        return(invisible(GdObject))
    }
    if((is.logical(debug) && debug) || debug=="draw")
        browser()
    ## We only proceed if there is something to draw within the ranges, but still may have to add the grid and the legend.
    ## Legend drawing causes another viewport for all the other graphics to be opened and will be called after all other
    ## drawing has finished, hence we call it in on.exit
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
    ## The optional legend is plotted below the data
    grpLevels <- .dpOrDefault(GdObject, ".__groupLevels")
    if(as.logical(.dpOrDefault(GdObject, "legend", FALSE)) && !is.null(grpLevels)){
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
                 fontcolor <- .dpOrDefault(GdObject, "fontcolor.legend", .DEFAULT_SHADED_COL)
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
                                   gp=gpar(fill=pcols$col[i], col=.DEFAULT_SHADED_COL))
                     } else {
                         if(any(c("pch", "col.symbol") %in% legFactors))
                             panel.points(unit(boxSize/2, "inches"), 0.5, pch=pcols$pch[i], cex=pcols$cex[i], col=pcols$col.symbol[i])
                         if(any(c("lwd", "lty", "col.lines") %in% legFactors))
                             ##panel.lines(unit(c(0,boxSize), "inches"), c(0.5, 0.5), col=pcols$col.line[i], lwd=pcols$lwd[i], lty=pcols$lty[i])
                             grid.lines(unit(c(0,boxSize), "inches"), c(0.5, 0.5), gp=gpar(col=pcols$col.line[i], lwd=pcols$lwd[i], lty=pcols$lty[i]))
                     }
                     grid.text(x=unit(boxSize+spacing, "inches"), y=0.5, just=c(0, 0.5), label=grpLevels[i], gp=gpar(col=fontcolor))
                     popViewport(1)
                 }
                 popViewport(2)
             })
    }
    if(!length(GdObject))
    {
        if ("g" %in% type)
            panel.grid(h=.dpOrDefault(GdObject, "h", -1), v=.dpOrDefault(GdObject, "v", -1),
                       col=.dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty=.dpOrDefault(GdObject, "lty.grid", 1),
                       lwd=.dpOrDefault(GdObject, "lwd.grid", 1), alpha=alpha)
        return(invisible(GdObject))
    }
    vals <- values(GdObject)
    ylim <- suppressWarnings(.dpOrDefault(GdObject, "ylim", range(vals, na.rm=TRUE, finite=TRUE)))
    if(diff(ylim)==0)
        ylim <- ylim+c(-1,1)
    if(all(is.infinite(ylim)))
        ylim <- c(0,1)
    ylimExt <- extendrange(r=ylim, f=0.05)
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
    ## The special type 'polygon' has to be handled separately
    if("polygon" %in% type) {
        mbaseline <- if(is.null(baseline)) 0 else baseline[1]
        fill.mountain <- .dpOrDefault(GdObject, "fill.mountain", superpose.symbol$fill)[1:2]
        col.mountain <- .dpOrDefault(GdObject, "col.mountain", pcols$col)[1]
        col.baseline <- .dpOrDefault(GdObject, "col.baseline", col.mountain)[1]
        lwd.mountain <- .dpOrDefault(GdObject, "lwd.mountain", pcols$lwd)[1]
        lty.mountain <- .dpOrDefault(GdObject, "lty.mountain", pcols$lty)[1]
        .panel.polygon(x, y, col=col.mountain, fill=fill.mountain, lwd=lwd.mountain,
                        lty=lty.mountain, col.line=col.mountain, alpha=alpha,
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
            grid.rect(rep(start(GdObject), each=nr), rev(yy), width=rep(width(GdObject), each=nr),
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
    if(!any(c("mountain","polygon") %in% type) && !is.null(baseline) && !is.na(baseline))
        panel.abline(h=baseline, col=pcols$col.baseline, lwd=lwd.baseline, lty=lty.baseline, alpha=alpha)
    popViewport(1)
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------------
## Draw a AlignedRead track
##----------------------------------------------------------------------------------------------------------------------------
setMethod("drawGD", signature("AlignedReadTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if((is.logical(debug) && debug) || debug=="prepare")
        browser()
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
    if((is.logical(debug) && debug) || debug=="draw")
        browser()
    ## We only proceed if there is something to draw within the ranges, but still may have to add the grid and the legend.
    ## Legend drawing causes another viewport for all the other graphics to be opened and will be called after all other
    ## drawing has finished, hence we call it in on.exit
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    alpha <- .dpOrDefault(GdObject, "alpha", 1)
    ## The optional legend is plotted below the data
    grpLevels <- .dpOrDefault(GdObject, ".__groupLevels")
    if(as.logical(.dpOrDefault(GdObject, "legend", FALSE)) && !is.null(grpLevels)){
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
                 fontcolor <- .dpOrDefault(GdObject, "fontcolor.legend", .DEFAULT_SHADED_COL)
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
                                   gp=gpar(fill=pcols$col[i], col=.DEFAULT_SHADED_COL))
                     } else {
                         if(any(c("pch", "col.symbol") %in% legFactors))
                             panel.points(unit(boxSize/2, "inches"), 0.5, pch=pcols$pch[i], cex=pcols$cex[i], col=pcols$col.symbol[i])
                         if(any(c("lwd", "lty", "col.lines") %in% legFactors))
                             ##panel.lines(unit(c(0,boxSize), "inches"), c(0.5, 0.5), col=pcols$col.line[i], lwd=pcols$lwd[i], lty=pcols$lty[i])
                             grid.lines(unit(c(0,boxSize), "inches"), c(0.5, 0.5), gp=gpar(col=pcols$col.line[i], lwd=pcols$lwd[i], lty=pcols$lty[i]))
                     }
                     grid.text(x=unit(boxSize+spacing, "inches"), y=0.5, just=c(0, 0.5), label=grpLevels[i], gp=gpar(col=fontcolor))
                     popViewport(1)
                 }
                 popViewport(2)
             })
    }
    if(!length(GdObject))
    {
        if ("g" %in% type)
            panel.grid(h=.dpOrDefault(GdObject, "h", -1), v=.dpOrDefault(GdObject, "v", -1),
                       col=.dpOrDefault(GdObject, "col.grid", "#e6e6e6"), lty=.dpOrDefault(GdObject, "lty.grid", 1),
                       lwd=.dpOrDefault(GdObject, "lwd.grid", 1), alpha=alpha)
        return(invisible(GdObject))
    }
    vals <- values(GdObject)
    ylim <- suppressWarnings(.dpOrDefault(GdObject, "ylim", range(vals, na.rm=TRUE, finite=TRUE)))
    if(diff(ylim)==0)
        ylim <- ylim+c(-1,1)
    if(all(is.infinite(ylim)))
        ylim <- c(0,1)
    ylimExt <- extendrange(r=ylim, f=0.05)
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
    ## The special type 'polygon' has to be handled separately
    if("polygon" %in% type) {
        mbaseline <- if(is.null(baseline)) 0 else baseline[1]
        fill.mountain <- .dpOrDefault(GdObject, "fill.mountain", superpose.symbol$fill)[1:2]
        col.mountain <- .dpOrDefault(GdObject, "col.mountain", pcols$col)[1]
        col.baseline <- .dpOrDefault(GdObject, "col.baseline", col.mountain)[1]
        lwd.mountain <- .dpOrDefault(GdObject, "lwd.mountain", pcols$lwd)[1]
        lty.mountain <- .dpOrDefault(GdObject, "lty.mountain", pcols$lty)[1]
        .panel.polygon(x, y, col=col.mountain, fill=fill.mountain, lwd=lwd.mountain,
                        lty=lty.mountain, col.line=col.mountain, alpha=alpha,
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
    if(!any(c("mountain","polygon") %in% type) && !is.null(baseline) && !is.na(baseline))
        panel.abline(h=baseline, col=pcols$col.baseline, lwd=lwd.baseline, lty=lty.baseline, alpha=alpha)
    popViewport(1)
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------------
## Draw a AlignedRead track
##----------------------------------------------------------------------------------------------------------------------------
setMethod("drawGD", signature("AlignedReadTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, subset=TRUE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if((is.logical(debug) && debug) || debug=="prepare")
        browser()
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
    if((is.logical(debug) && debug) || debug=="draw")
        browser()
    ## In plotting mode we either show all the reads (time-consuming), or the coverage only
    rad <- 0.015
    xx <- -0.01
    loc <- vpLocation()$size
    diff <- .pxResolution(coord="x")
    radv <- rad / if(loc["width"] < loc["height"]) c(1,loc[2]/loc[1]) else c(loc[1]/loc[2],1)
    if(subset)
        GdObject <- subset(GdObject, from=minBase, to=maxBase)
    ## If type is 'coverage' all we need to do is compute a coverage vector, create dummy DataTracks and pass everything on
    if(detail=="coverage")
    {
      if (!any(unlist(lapply(GdObject@coverage, function(y) runValue(y)!=0))))
        {
        ## Nothing there, but we still need the strand separator
        panel.abline(h=0.5, col="lightgray", lwd=2)
        grid.circle(xx, c(0.25, 0.75), rad, gp=gpar(fill="lightgray", col="lightgray"))
        grid.segments(c(rep(xx-radv[1]+(radv[1]/2),2), xx), c(0.25, 0.75, 0.75-(radv[2]/2)),
                      c(rep(xx+radv[1]-(radv[1]/2),2), xx), c(0.25, 0.75, 0.75+radv[2]/2),
                      gp=gpar(col="white", lwd=2, lineend="square"), default.units="native")
        return(invisible(GdObject))
      } else {
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
            sel <- suppressWarnings(runValue(cov)!=0) #changed from >
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
    }
    if(detail=="reads")
    {
        if(!length(GdObject))
        {
            ## No reads, but we still need the strand separator
            panel.abline(h=0.5, col="lightgray", lwd=2)
            grid.circle(xx, c(0.25, 0.75), rad, gp=gpar(fill="lightgray", col="lightgray"))
            grid.segments(c(rep(xx-radv[1]+(radv[1]/2),2), xx), c(0.25, 0.75, 0.75-(radv[2]/2)),
                          c(rep(xx+radv[1]-(radv[1]/2),2), xx), c(0.25, 0.75, 0.75+radv[2]/2),
                          gp=gpar(col="white", lwd=2, lineend="square"), default.units="native")
            return(invisible(GdObject))
        } else {  
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
    }
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------		
		



##----------------------------------------------------------------------------------------------------------------------------
## Draw an ideogram track
##----------------------------------------------------------------------------------------------------------------------------
## Helper function to compute coordinates for a rounded ideogram cap
.roundedCap <- function(bl, tr, st, vals, side=c("left", "right"), bevel=0.4, n=100)
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
    lcS <- split(as.data.frame(coords), cut(coords[,1], st, right=TRUE, include.lowest=TRUE, labels=FALSE), drop=TRUE)
    first <- TRUE
    shift <- ifelse(side=="left", 0, 1)
    for(j in names(lcS)){
        xx <- lcS[[j]][,1]
        yy <- lcS[[j]][,2]
        if(!first){
            prev <- lcS[[as.character(as.numeric(j)-1)]]
            xx <- c(tail(prev[,1], 1), xx, tail(prev[,1], 1))
            yy <- c(1-tail(prev[,2], 1), yy, tail(prev[,2], 1))
        }
        grid.polygon(xx, yy, gp=gpar(col=vals[as.numeric(j)+shift,"col"], fill=vals[as.numeric(j)+shift,"col"]))
        first <- FALSE
    }
    return(coords)
}

## A more generic method to come up with colors for chromosome bands that still relies a bit on biovizBase
.getBioColorIdeo <- function(type){
    type <- as.character(type)
    ocols <- getBioColor("CYTOBAND")
    cols <- c(ocols[c("gneg", "stalk", "acen")], gpos=unname(ocols["gpos100"]), gvar=unname(ocols["gpos100"]))
    gpcols <- unique(grep("gpos", type, value=TRUE))
    crmp <- colorRampPalette(c(cols["gneg"], cols["gpos"]))(100)
    posCols <- setNames(crmp[as.integer(gsub("gpos", "", gpcols))], gpcols)
    return(c(cols, posCols))
}

## The actual drawing method
setMethod("drawGD", signature("IdeogramTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if((is.logical(debug) && debug) || debug=="prepare")
        browser()
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
    if((is.logical(debug) && debug) || debug=="draw")
        browser()
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
    ## Color mapping for the bands taken from the biovizBase package
    cols <- .getBioColorIdeo(values(GdObject)$type)
    vals <- data.frame(values(GdObject), col=cols[as.character(values(GdObject)$type)], stringsAsFactors=FALSE)
    ## For the rounded caps we need  to figure out the overlap with existing bands for proper coloring
    bevel <- 0.02
    ol <- queryHits(findOverlaps(range(GdObject), IRanges(start=c(bevel, 1-bevel)*len, width=1)))
    st <- start(range(GdObject))/len
    ed <- end(range(GdObject))/len
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
    lc <- .roundedCap(c(stExt[1], margin), c(stExt[ls], 1-margin), st, vals, side="left", bevel=.dpOrDefault(GdObject, "bevel", 0.45))
    rc <- .roundedCap(c(tail(stExt,1), margin), c(1, 1-margin), ed, vals, side="right", bevel=.dpOrDefault(GdObject, "bevel", 0.45))
    lcol <- "black"; lwd <- 1; lty <- 1
    ## Now some outlines
    grid.lines(lc[,1], lc[,2], gp=gpar(col=lcol, lwd=lwd, lty=lty))
    grid.lines(rc[,1], rc[,2],gp= gpar(col=lcol, lwd=lwd, lty=lty))
    if(length(cent))
    {
        x0 <- c(rep(max(lc[,1]), 2), rep(stExt[max(cent)+1],2), rep(stExt[min(cent)],2), rep(stExt[max(cent)],2))
        y0 <-  c(rep(c(margin, (1-margin)), 3), 0.5, 0.5)  
        x1 <- c(rep(stExt[min(cent)],2), rep(tail(stExt, 1),2), rep(stExt[min(cent)]+wd[min(cent)],2),
                rep(stExt[max(cent)]+wd[max(cent)],2))
        y1 <-  c(rep(c(margin, (1-margin)), 2), 0.5, 0.5, margin, (1-margin))
    } else {
        x0 <- rep(max(lc[,1]), 2)
        y0 <- c(margin, (1-margin))
        x1 <- rep(max(stExt), 2)
        y1 <- y0
    }
    grid.segments(x0, y0, x1, y1, gp=gpar(col=lcol, lwd=lwd, lty=lty))
    ## The outlines of the box
    if(!missing(minBase) && !missing(maxBase))
    {
        col <- .dpOrDefault(GdObject, "col", "red")
        lwd <- .dpOrDefault(GdObject, "lwd", 1)
        lty <- .dpOrDefault(GdObject, "lty", "solid")
        grid.rect(minBase/len, 0.1, width=min(1,(maxBase-minBase)/len), height=0.8, just=c("left","bottom"),
                  gp=gpar(col=col, fill="transparent", lwd=lwd, lty=lty))
    }
    ## Finally the band annotation if we need it
    if(.dpOrDefault(GdObject, "showBandId", FALSE)){
        bn <- as.character(values(GdObject)$name)
        cval <- rgb2hsv(col2rgb(cols[as.character(values(GdObject)$type)]))["v", ]
        tcol <- ifelse(cval>0.9, "black", "white")
        bwidth <- (c(st[-1], 1)-st)/2
        cex.bands <- .dpOrDefault(GdObject, "cex.bands", 0.7)
        sspace <- as.numeric(convertUnit(unit(0.01, "inches"), "native"))
        swidth <- as.numeric(convertWidth(stringWidth(bn), "native"))*cex.bands+sspace
        sel <- swidth<bwidth
        if(any(sel))
            grid.text(x=(st+bwidth)[sel], y=0.5, label=bn[sel], hjust=0.5, gp=gpar(col=tcol[sel], cex=cex.bands))
    }
    popViewport(1)
    return(invisible(GdObject))
})
##----------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------
## Draw a SequenceTrack
##----------------------------------------------------------------------------------------------------------------------------
setMethod("drawGD", signature("SequenceTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...) {
    debug <- .dpOrDefault(GdObject, "debug", FALSE)
    if((is.logical(debug) && debug) || debug=="prepare")
        browser()
    fcol <- .dpOrDefault(GdObject, "fontcolor", getBioColor("DNA_BASES_N"))
    cex <- max(0.3, .dpOrDefault(GdObject, "cex", 1))
    pushViewport(viewport(xscale=c(minBase, maxBase), clip=TRUE,
                          gp=gpar(alpha=.dpOrDefault(GdObject, "alpha", 1),
                                  fontsize=.dpOrDefault(GdObject, "fontsize", 12),
                                  fontface=.dpOrDefault(GdObject, "fontface", 2),
                                  lineheight=.dpOrDefault(GdObject, "lineheight", 1),
                                  fontfamily=.dpOrDefault(GdObject, "fontfamily", 1),
                                  cex=cex)))
    if(prepare){
        pres <- .pxResolution()
        nsp <-  max(as.numeric(convertHeight(stringHeight(stringWidth(DNA_ALPHABET)),"native")))
        nsp <- nsp/pres["y"]*2
        displayPars(GdObject) <- list("neededVerticalSpace"=nsp)
        popViewport(1)
        return(invisible(GdObject))
    }
    if((is.logical(debug) && debug) || debug=="draw")
        browser()
    imageMap(GdObject) <- NULL
    delta <- maxBase-minBase
    if(delta==0)
        return(invisible(GdObject))
    lwidth <- max(as.numeric(convertUnit(stringWidth(DNA_ALPHABET),"inches")))
    perLetter <- vpLocation()$isize["width"]/(maxBase-minBase+1)
    diff <- .pxResolution(.dpOrDefault(GdObject, "min.width", 2), coord="x")
    ## FIXME: Need to deal with sequences that are too long.
    if(diff>1 || (maxBase-minBase+1)>=10e6){
        grid.lines(x=unit(c(minBase, maxBase), "native"), y=0.5,
                   gp=gpar(col=.dpOrDefault(GdObject, "col", "darkgray"),
                               lwd=.dpOrDefault(GdObject, "lwd", 2)))
    }else{
       
        sequence <- as.character(as(subseq(GdObject, start=minBase, end=maxBase-1), "Rle"))
        at <- seq((minBase+0.5), maxBase - 1 + 0.5, by=1)
        sequence[sequence=="-"] <- ""
        if(perLetter<0.5 && .dpOrDefault(GdObject, "add53", FALSE))
            sequence[c(1, length(sequence))] <- ""
        col <- fcol[toupper(sequence)]
        if(lwidth<perLetter && !.dpOrDefault(GdObject, "noLetters", FALSE)){
            grid.text(x=unit(at, "native"), y=0.5, label=sequence, rot=.dpOrDefault(GdObject, "rotation", 0),
                      gp=gpar(col=col))
        } else {
            grid.rect(x=unit(at, "native"), y=0.05, width=unit(1, "native"), height=0.9,
                      gp=gpar(fill=col, col="white"), just=c(0.5, 0))
        }
    }
    ## The direction indicators
    if(.dpOrDefault(GdObject, "add53", FALSE))
    {
        if(.dpOrDefault(GdObject, "complement", FALSE)){
            grid.text(label=expression("3'"), x=unit(minBase+0.1, "native"), just=c(0, 0.5), gp=gpar(col="#808080", cex=0.8))
            grid.text(label=expression("5'"), x=unit(maxBase-0.1, "native"), just=c(1, 0.5), gp=gpar(col="#808080", cex=0.8))
        }else{
            grid.text(label=expression("5'"), x=unit(minBase+0.1, "native"), just=c(0, 0.5), gp=gpar(col="#808080", cex=0.8))
            grid.text(label=expression("3'"), x=unit(maxBase-0.1, "native"), just=c(1, 0.5), gp=gpar(col="#808080", cex=0.8))
        }
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
                                          genome=genome(from), strand=strand,
                                          asRangedData=TRUE),
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
                                          genome=genome(from), strand=strand,
                                          asRangedData=TRUE),
                  trackLine = line)   
          })
          
setAs("RangeTrack", "data.frame",
          function(from, to) as(as(ranges(from), "DataFrame"), "data.frame"))

setAs("DataTrack", "data.frame", function(from, to){
    tmp <- as.data.frame(ranges(from))
    colnames(tmp)[1] <- "chromosome"
    tmp <- cbind(tmp, as.data.frame(t(values(from, all=TRUE))))
    return(tmp)
})

setAs("GRanges", "DataTrack", function(from, to) DataTrack(range=from))

setAs("GRanges", "AnnotationTrack", function(from, to) AnnotationTrack(range=from))
setAs("GRangesList", "AnnotationTrack", function(from, to) AnnotationTrack(range=from))

setAs("GRanges", "GeneRegionTrack", function(from, to) GeneRegionTrack(range=from))
setAs("GRangesList", "GeneRegionTrack", function(from, to) GeneRegionTrack(range=from))
setAs("TranscriptDb", "GeneRegionTrack", function(from, to) GeneRegionTrack(range=from))

setAs("DNAString", "Rle", function(from, to) Rle(strsplit(as.character(from), "")[[1]]))

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
## Helper methods to build a GRanges object from the input arguments.
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
                stop("Number of elements in argument '", a, "' is invalid")
            if(!is.null(by) && length(val)!=sum(by)) rep(val, by) else val
        }
    }
    return(range)
}

## For numeric vectors we can immediately create a data frame after some sanity checking and pass that on to the next method.
setMethod(".buildRange", signature("NULLOrMissing", "NumericOrNULL", "NumericOrNULL", "NumericOrNULL"),
          function(range, start, end, width, asIRanges=FALSE, by=NULL, len, args, defaults, ...){
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
              return(.buildRange(range=range, asIRanges=asIRanges, args=args["genome"], defaults=defaults, ...))})

## For data.frames we need to check for additional arguments (like feature, group, etc.), the chromosome information
## and create the final GRanges object
setMethod(".buildRange", signature("data.frame"),
          function(range, asIRanges=FALSE, args=list(), defaults=list(), chromosome=NULL, trackType, ...){
              if(asIRanges){
                  range <- .fillWithDefaults(range, defaults, args, len=nrow(range), ignore=setdiff(names(defaults), c("start", "end", "genome")))
                  return(IRanges(start=as.integer(range$start), end=as.integer(range$end)))
              }
              mandArgs <- c("start", "end", "genome", names(defaults))
              ## Not quite sure how whether exisiting chromosome information in a GRanges object should generally have precedence over the
              ## chromosome constructor, but probably that should be the case
              if("chromosome" %in% colnames(range))
                  args$chromosome <- NULL
              missing <- setdiff(union(setdiff(mandArgs, c(colnames(range))), names(which(!sapply(args, is.null)))), "genome")
              range <- .fillWithDefaults(range, defaults[missing], args[missing], len=nrow(range))
              range$chromosome <- .chrName(as.character(range$chromosome))
              grange <- GRanges(ranges=IRanges(start=range$start, end=range$end), strand=range$strand, seqnames=range$chromosome)
              mcols(grange) <- range[,setdiff(colnames(range), c("start", "end", "strand", "width", "chromosome", "genome", "seqnames",
                                                                           "ranges", "seqlevels", "seqlengths", "isCircular", "element"))]
              if(trackType != "DataTrack")
                  mcols(grange) <- mcols(grange)[, intersect(names(defaults), colnames(mcols(grange)))]
              suppressWarnings(genome(grange) <- unname(if(is.null(args[["genome"]])) defaults[["genome"]] else as.character(args[["genome"]])[[1]]))
              return(grange)})


## For GRanges we just need to check for the existence of additional arguments (like feature, group, etc.)
setMethod(".buildRange", signature("GRanges"),
          function(range, asIRanges=FALSE, args=list(), defaults=list(), trackType=NULL, ...){
              if(asIRanges)
                  return(ranges(range))
              if(length(range))
              {
                  mandArgs <- names(defaults)
                  ## Not quite sure how whether exisiting chromosome information in a GRanges object should generally have precedence over the
                  ## chromosome constructor, but probably that should be the case
                  args$chromosome <- NULL
                  range <- renameSeqlevels(range, setNames(.chrName(seqlevels(range)), seqlevels(range)))
                  missing <- setdiff(union(setdiff(mandArgs, c("chromosome", "strand", colnames(mcols(range)))), names(which(!sapply(args, is.null)))), "genome")
                  newVars <- .fillWithDefaults(DataFrame(chromosome=as.character(seqnames(range)), strand=as.character(strand(range)), mcols(range), check.names=FALSE),
                                               defaults[missing], args[missing], len=length(range))
                  if(any(c("start", "end", "strand", "chromosome") %in% colnames(newVars))){
                      gen <- genome(range)
                      range <- GRanges(seqnames=if(is.null(newVars[["chromosome"]])) seqnames(range) else (newVars[["chromosome"]]),
                                       strand=if(is.null(newVars[["strand"]])) strand(range) else (newVars[["strand"]]),
                                       ranges=IRanges(start=if(is.null(newVars[["start"]])) start(range) else (newVars[["start"]]),
                                                      end=if(is.null(newVars[["end"]])) end(range) else (newVars[["end"]])))
                      if(length(unique(gen)) != 1)
                          warning("Tracks can only be defined for a single genome. Forcing all reads to belong to genome '", gen[1], "'")
                      defaults[["genome"]] <- as.character(gen)[1]
                  }
                  mcols(range) <- newVars[, setdiff(colnames(newVars),  c("start", "end", "strand", "width", "chromosome", "genome", "seqnames",
                                                                          "ranges", "seqlevels", "seqlengths", "isCircular", "element")), drop=FALSE]
              }
              if(trackType != "DataTrack")
                  mcols(range) <- mcols(range)[, intersect(names(defaults), colnames(mcols(range)))]
              ## The genome information may or may not be encoded in the GRanges object at this time but we want it in there for sure
              genome <- if(!is.null(args[["genome"]])) args[["genome"]] else .getGenomeFromGRange(range, defaults[["genome"]])
              suppressWarnings(genome(range) <- unname(genome))[1]
              return(range)})

## For IRanges we need to deal with additional arguments (like feature, group, etc.) and create the final GRanges object
setMethod(".buildRange", signature("IRanges"),
          function(range, asIRanges=FALSE, args=list(), defaults=list(), chromosome=NULL, strand, ...){
              if(asIRanges)
                  return(range)
              if(missing(chromosome) || is.null(chromosome))
                  stop("Unable to find chromosome information in any of the arguments")
              range <- GRanges(seqnames=.chrName(chromosome), range=range, strand=if(!is.null(args$strand)) args$strand else "*")
              if(length(range))
              {
                  vals <- .fillWithDefaults(defaults=defaults, args=args, len=(length(range)), by=NULL, ignore="strand")
                  elementMetadata(range) <- vals
              }
              return(range)})

## For GRangesLists we capture the grouping information from the list structure, unlist and use the GRanges method
setMethod(".buildRange", signature("GRangesList"),
          function(range, groupId="group", ...){
              grps <- rep(names(range), elementLengths(range))
              range <- unlist(range)
              names(range) <- NULL
              mcols(range)[[groupId]] <- grps
              return(.buildRange(range=range, ...))})

## For TranscriptDb objects we extract the grouping information and use the GRanges method
setMethod(".buildRange", signature("TranscriptDb"),
          function(range, groupId="transcript", tstart, tend, chromosome, args, ...){
              ## If chromosome (and optional start and end) information is present we only extract parts of the annotation data
              noSubset <- is.null(tstart) && is.null(tend)
              if(!is.null(chromosome)){
                  chromosome <- .chrName(chromosome)
                  ## Seems like TranscriptDb objects use pass by reference for the active chromosomes, so we have to
                  ## restore the old values after we are done
                  oldAct <- isActiveSeq(range)
                  oldRange <- range
                  on.exit(isActiveSeq(oldRange)[seqlevels(oldRange)] <- oldAct)
                  isActiveSeq(range)[seqlevels(range)] <- FALSE
                  isActiveSeq(range) <- structure(rep(TRUE, length(unique(chromosome))), names=unique(chromosome))
                  sl <- seqlengths(range)
                  if(is.null(tstart))
                      tstart <- rep(1, length(chromosome))
                  if(is.null(tend)){
                      tend <- sl[chromosome]+1
                      tend[is.na(tend)] <- tstart[is.na(tend)]+1
                  }
                  sRange <- GRanges(seqnames=chromosome, ranges=IRanges(start=tstart, end=tend))
              }
              ## First the mapping of internal transcript ID to transcript name
              txs <- as.data.frame(values(transcripts(range, columns = c("tx_id", "tx_name"))))
              rownames(txs) <- txs[, "tx_id"]
              ## Now the CDS ranges
              t2c <- cdsBy(range, "tx")
              names(t2c) <- txs[names(t2c), 2]
              tids <- rep(names(t2c), elementLengths(t2c))
              t2c <- unlist(t2c)
              if(length(t2c)){
                  t2c$tx_id <- tids
                  t2c$feature_type <- "CDS"
              }
              ## And the 5'UTRS
              t2f <- fiveUTRsByTranscript(range)
              names(t2f) <- txs[names(t2f), 2]
              tids <- rep(names(t2f), elementLengths(t2f))
              t2f <- unlist(t2f)
              if(length(t2f)){
                  t2f$tx_id <- tids
                  t2f$feature_type <- "utr5"
              }
              ## And the 3'UTRS
              t2t <- threeUTRsByTranscript(range)
              names(t2t) <- txs[names(t2t), 2]
              tids <- rep(names(t2t), elementLengths(t2t))
              t2t <- unlist(t2t)
              if(length(t2t)){
                  t2t$tx_id <- tids
                  t2t$feature_type <- "utr3"
              }
              ## And finally all the non-coding transcripts
              nt2e <- exonsBy(range, "tx")
              names(nt2e) <- txs[names(nt2e), 2]
              nt2e <- nt2e[!names(nt2e) %in% c(values(t2c)$tx_id, values(t2f)$tx_id, values(t2t)$tx_id)]
              tids <- rep(names(nt2e), elementLengths(nt2e))
              nt2e <- unlist(nt2e)
              if(length(nt2e)){
                  nt2e$tx_id <- tids
                  nt2e$feature_type <- "ncRNA"
              }
              ## Now we can merge the three back together (we need to change the column names of t2c to make them all the same)
              colnames(values(t2c))[1:2] <- c("exon_id", "exon_name")
              ## t2e <- c(t2c, t2f, t2t, nt2e) ## This is super-slow, much more efficient if we build the GRanges object from the individual bits and pieces
              vals <- DataFrame(exon_id=c(values(t2c)$exon_id, values(t2f)$exon_id, values(t2t)$exon_id, values(nt2e)$exon_id),
                                exon_name=c(values(t2c)$exon_name, values(t2f)$exon_name, values(t2t)$exon_name, values(nt2e)$exon_name),
                                exon_rank=c(values(t2c)$exon_rank, values(t2f)$exon_rank, values(t2t)$exon_rank, values(nt2e)$exon_rank),
                                tx_id=c(values(t2c)$tx_id, values(t2f)$tx_id, values(t2t)$tx_id, values(nt2e)$tx_id),
                                feature_type=c(values(t2c)$feature_type, values(t2f)$feature_type, values(t2t)$feature_type, values(nt2e)$feature_type))
              t2e <- GRanges(seqnames=c(seqnames(t2c), seqnames(t2f), seqnames(t2t), seqnames(nt2e)),
                             ranges=IRanges(start=c(start(t2c), start(t2f), start(t2t), start(nt2e)),
                                            end=c(end(t2c), end(t2f), end(t2t), end(nt2e))),
                             strand=c(strand(t2c), strand(t2f), strand(t2t), strand(nt2e)))
              values(t2e) <- vals
              if(length(t2e)==0)
                  return(GRanges())
              ## Add the gene level annotation
              g2t <- transcriptsBy(range, "gene")
              gids <- rep(names(g2t), elementLengths(g2t))
              g2t <- unlist(g2t)
              values(g2t)[["gene_id"]] <- gids
              values(t2e)$gene_id <- gids[match(values(t2e)$tx_id, as.character(txs[as.character(values(g2t)$tx_id),2]))]
              vals <- values(t2e)[c("tx_id", "exon_id", "exon_rank", "feature_type", "tx_id", "gene_id")]
              colnames(vals) <- c("transcript", "exon", "rank", "feature", "symbol", "gene")
              ## Add the genome information
              genome(t2e) <- unique(genome(range))
              ## Finally we re-assign, subset if necessary, and sort
              range <- t2e
              values(range) <- vals
              if(!noSubset && !is.null(chromosome)){
	      	  ## We have to keep all exons for all the overlapping transcripts
	          txSel <- unique(subsetByOverlaps(g2t, sRange)$tx_name)
                  range <- range[range$transcript %in% txSel]
 	      }		  
              args <- list(genome=genome(range)[1])
              return(.buildRange(range=sort(range), chromosome=chromosome, args=args, ...))})


## For character scalars the data need to be extracted from a file and we have to deal with parser functions
## and column assignments here. 
setMethod(".buildRange", signature("character"),
          function(range, importFun=NULL, trackType, stream=FALSE, args, defaults, ...){
              .checkClass(range, "character", 1)
              .checkClass(importFun, c("NULL", "function"), mandatory=FALSE)
              .checkClass(stream, "logical", 1)
              ## We first check for the default column mapping and whether this is a streaming file
              defMap <- .defaultVarMap(.fileExtension(range), trackType, stream, !is.null(importFun))
              isStream <- !is.null(defMap[[".stream"]]) && defMap[[".stream"]]
              defMap[[".stream"]] <- NULL
              if(!isStream){
                  data <- if(is.null(importFun)) .registerImportFun(range) else{
                      if(!"file" %in% names(formals(importFun)))
                          stop("The user-defined import function needs to define a 'file' argument")
                      importFun(range)
                  }
                  if(!is(data, "GRanges"))
                      stop("The import function did not provide a valid GRanges object. Unable to build track from file '",
                           range, "'")
                  if(trackType=="DataTrack"){
                      ## For data tracks we take all numeric data columns regardless of any mapping
                      mc <- .prepareDtData(as.data.frame(mcols(data)), length(data))
                      mcols(data) <- t(mc)
                  } else {
                      ## For the rest we use the mapping as provided by the constructor
                      ## are available
                      cmap <- .resolveColMapping(data, args, defMap)
                      args <- cmap$args
                      data <- cmap$data
                  }
                  args[["chromosome"]] <- as.character(seqnames(data))
                  args[["strand"]] <- as.character(strand(data))
                  return(.buildRange(range=data, args=args, defaults=defaults, trackType=trackType, ...))
              }else{
                  if(trackType!="DataTrack"){
                      for(i in names(defMap)){
                          if(is.character(args[[i]]) && length(args[[i]])==1){
                              defMap[[i]] <- args[[i]]
                          }
                      }
                  }
                  return(list(reference=path.expand(range), mapping=defMap,
                              stream=if(is.null(importFun)) .registerImportFun(range) else importFun))
              }
          })
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

## A helper function to plot information regarding additional features on other chromosomes
.addFeatInfo <- function(object, addfeat){
    freqs <- table(seqnames(object))
    freqs <- freqs[setdiff(names(freqs), chromosome(object))]
    nrChr <- length(freqs)
    msg <- sprintf("There %s %s additional annotation feature%s on %s further chromosome%s%s",
                   ifelse(addfeat>1, "are", "is"),
                   addfeat,
                   ifelse(addfeat>1, "s", ""),
                   nrChr,
                   ifelse(nrChr>1, "s", ""),
                   ifelse(nrChr==1, sprintf(" (%s)", names(freqs)), ""))
    if(nrChr>1){
        msg <- if(nrChr>10){
            c(msg, paste("  ", head(names(freqs), 5), ": ", head(freqs, 5), sep="", collapse="\n"),
              "  ...", paste("  ", tail(names(freqs), 5), ": ", tail(freqs, 5), sep="", collapse="\n"))
        }else{
            c(msg, paste("  ", names(freqs), ": ", freqs, " features", sep="", collapse="\n"))
        }
        msg <- c(msg, paste("Call seqlevels(obj) to list all available chromosomes",
                            "or seqinfo(obj) for more detailed output"))
    }
    return(msg)
}

setMethod("show",signature(object="DataTrack"),
          function(object){
              msg <- sprintf(paste("DataTrack '%s'\n| genome: %s\n| active chromosome: %s\n",
                                   "| positions: %s\n| samples:%s\n| strand: %s", sep=""),
                             names(object),
                             genome(object),
                             chromosome(object),
                             length(object),
                             nrow(values(object)),
                             strand(object)[1])
              addfeat <- ncol(object@data)-length(object)
              if(addfeat>0)
                  msg <- c(msg, .addFeatInfo(object, addfeat), "Call chromosome(obj) <- 'chrId' to change the active chromosome")
              cat(paste(msg, collapse="\n"), "\n")
          })


## A helper function to plot general information about an AnnotationTrack
.annotationTrackInfo <- function(object){
    msg <- sprintf(paste("| genome: %s\n| active chromosome: %s\n",
                         "| annotation features: %s", sep=""),
                   genome(object),
                   chromosome(object),
                   length(object))
    addfeat <- length(object@range)-length(object)
    if(addfeat>0)
        msg <- c(msg, .addFeatInfo(object, addfeat), "Call chromosome(obj) <- 'chrId' to change the active chromosome")
    return(paste(msg, collapse="\n"))
}

## We have to show the name, genome and currently active chromosome, and, if more ranges are available on additional
## chromosomes some information about that
setMethod("show", signature(object="AnnotationTrack"), function(object)
          cat(sprintf("AnnotationTrack '%s'\n%s\n", names(object), .annotationTrackInfo(object))))

setMethod("show", signature(object="GeneRegionTrack"), function(object)
          cat(sprintf("GeneRegionTrack '%s'\n%s\n", names(object), .annotationTrackInfo(object))))

## A helper function to plot general information about a ReferenceTrack
.referenceTrackInfo <- function(object, type){
    cat(sprintf("%s '%s'\n| genome: %s\n| active chromosome: %s\n| referenced file: %s\n",
                type,
                names(object),
                genome(object),
                chromosome(object),
                object@reference))
    if(length(object@mapping) && type != "ReferenceDataTrack")
        cat(sprintf( "| mapping: %s\n",  paste(names(object@mapping), as.character(object@mapping), sep="=", collapse=", ")))
}

setMethod("show",  signature(object="ReferenceAnnotationTrack"), function(object)
    .referenceTrackInfo(object, "ReferenceAnnotationTrack"))

setMethod("show",  signature(object="ReferenceDataTrack"), function(object)
    .referenceTrackInfo(object, "ReferenceDataTrack"))

setMethod("show",  signature(object="ReferenceGeneRegionTrack"), function(object)
    .referenceTrackInfo(object, "ReferenceGeneRegionTrack"))

setMethod("show",  signature(object="ReferenceSequenceTrack"), function(object)
    .referenceTrackInfo(object, "ReferenceSequenceTrack"))
  
setMethod("show", signature(object="GenomeAxisTrack"),
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

## A helper function to print general information about SequenceTracks
.sequenceTrackInfo <- function(object){
    msg <- sprintf(paste("Sequence track '%s':\n",
                         "| genome: %s\n",
                         "| chromosomes: %s\n",
                         "| active chromosome: %s (%s nulceotides)\n", sep=""),
                   names(object),
                   genome(object),
                   length(seqnames(object)),
                   chromosome(object),
                   length(object))
    if(length(seqnames(object))>1)
        msg <- paste(msg, "Call seqnames() to list all available chromosomes\n",
                     "Call chromosome()<- to change the active chromosome\n", sep="")
    return(msg)
}

## We need to show the name, genome, information about the source BSgenome object as well as the currently active chromosome
setMethod("show",signature(object="SequenceBSgenomeTrack"),
		  function(object){
                      cat(.sequenceTrackInfo(object),
                          sprintf(paste("Parent BSgenome object:\n",
                                        "| organism: %s (%s)\n",
                                        "| provider: %s\n",
                                        "| provider version: %s\n",
                                        "| release date: %s\n",
                                        "| release name: %s\n",
                                        "| package name: %s\n", sep=""),
                                  organism(object@sequence),
                                  object@sequence@species,
                                  provider(object@sequence),
                                  providerVersion(object@sequence),
                                  releaseDate(object@sequence),
                                  releaseName(object@sequence),
                                  object@sequence@seqs_pkgname), sep="")
                  })

## Here we only need the name, genome and currently active chromosome information
setMethod("show", signature(object="SequenceDNAStringSetTrack"), function(object) cat(.sequenceTrackInfo(object)))
  
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
