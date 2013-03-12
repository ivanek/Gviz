## All class definitions and constructors/initializers for the
## package. Note that some methods have to be defined in here as well
## because they are used as part of the initializers.  Constructors
## and definitions that are exported in the name space are marked with
## (N)



##----------------------------------------------------------------------------------------------------------------------
## Some usefull class unions to allow for NULL values in certain slots
##----------------------------------------------------------------------------------------------------------------------
setClassUnion("DfOrNULL", c("data.frame", "NULL"))
setClassUnion("MartOrNULL", c("Mart", "NULL"))
setClassUnion("FactorOrCharacterOrNULL", c("factor", "character", "NULL"))
setClassUnion("NumericOrNULL", c("numeric", "NULL"))
setClassUnion("ListOrEnv", c("list", "environment"))
setClassUnion("GRangesOrIRanges", c("GRanges", "IRanges"))
setClassUnion("NULLOrMissing", c("NULL", "missing"))
setClassUnion("BSgenomeOrNULL", c("BSgenome", "NULL"))
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## ImageMap:
##
## A class to hold HTML image map information
## Slots :
##    o coords: a numeric matrix of rectangular image map cordinates, in the order x bl, y bl, x tr, y tr. Rownames are
##       mandatory for the matrix and have to be unique
##    o tags: a list of tags that are to be added to the <map> tag, where the name of the list item is used as the
##      tagname. The value of each list item has to be a named character vector, where the names must match back into
##      the rownames of the 'coords' matrix
##----------------------------------------------------------------------------------------------------------------------
setClass("ImageMap", representation(coords="matrix", tags="list"),
         prototype=prototype(coords=matrix(1, ncol=4, nrow=0), tags=list()))

## Constructor
ImageMap <- function(coords, tags)
{
    if(!(is.matrix(coords) && is.numeric(coords) && ncol(coords)==4))
        stop("'coords' must be a numeric matrix with 4 columns")
    rn <- rownames(coords)
    if(is.null(rn))
        stop("Rownames must be set for the matrix in 'coords'")
    if(!is.list(tags) || is.null(names(tags)) || any(names(tags)=="") || !all(sapply(tags, is.character)))
        stop("'tags' must be a named list with character vector items.")
    n <- unique(unlist(sapply(tags, names)))
    if(is.null(n) || any(n==""))
        stop("All items in the 'tags' list must be named character vectors.")
    m <- n %in% rn
    if(!all(m))
        stop("The following values in the 'tags' list could not be mapped to the 'coords' matrix:\n",
             paste(n[!m], sep="", collapse=", "))
    new("ImageMap", coords=coords, tags=tags)
}

## Allow for NULL value slots
setClassUnion("ImageMapOrNULL", c("ImageMap", "NULL"))
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## DisplayPars:
##
## A class to control the plotting parameters for GdObjects
## Slots :
##    o pars: an environment or a list containing parameter key value pairs
##       The initial idea was for the class to uses pass by reference semantic, allowing to update the object
##       without reassigning to a symbol. However this turned out to be problematic whenever a GdObject was
##       copied and to avoid confusion with the users about unwanted side-effects we decided to deprecate this
##       feature.
##----------------------------------------------------------------------------------------------------------------------
setClass("DisplayPars", representation(pars="ListOrEnv"))

## The initializer needs to create the environment, we can't take this
## directly from the class prototype since this would result in using
## the same environment all the time
setMethod("initialize", "DisplayPars", function(.Object, ...) {
  e = new.env(hash=TRUE)
  args <- list(...)
  n <- names(args)
  if(any(n == ""))
    stop("All supplied arguments must be named.")
  .Object@pars <- args
  return(.Object)
})

## Constructor, all supplied arguments are added to the environment.
DisplayPars <- function(...)
{
  return(new("DisplayPars", ...))
}

## Update function and deprecation message to show that an old environment-based
## DisplayParameter object has been updated to a list-based object.
.updateDp <- function(x, interactive=TRUE)
{
    if(interactive)
        message("Note that the behaviour of the 'setPar' method has changed. You need to reassign the result to an ",
                "object for the side effects to happen. Pass-by-reference semantic is no longer supported.")
    if(class(x@pars)=="environment")
    {
        x@pars <- as.list(x@pars) 
        message("The DisplayPars object has been updated to a list-based representation.")
    }
    return(x)
}

## The accessor methods for the DisplayPars class. (they need to be in
## here because they are being called in the prototypes). For the
## setter method, the input can either be a named list, or a single
## keyword/value pair. The DisplayPars class was first implemented as
## an environment, so we essentially had pass by reference semantic here,
## and the object could be modified without the need to reassign to a
## symbol. However this somewhat broke the the R paradigm, and 'setPar'
## has been deprecated in favour of the more standard 'DisplayPars<-'
## replacement method. The getter methods either return the list
## of all parameters, or a subset of parameters if their names are
## provided as a character vector.  Please note that for convenience
## the result is unlisted if only a single parameter is queried.
setMethod("setPar", signature("DisplayPars", "list"), 
          function(x, value, interactive=TRUE) {
              x <- .updateDp(x, interactive)
              x@pars[names(value)] <- value
              return(x)
          })

setMethod("setPar", signature("DisplayPars", "character"), 
          function(x, name, value, interactive=TRUE) {
              x <- .updateDp(x, interactive)
              x@pars[[name]] <- value
              return(x)
          })

setReplaceMethod("displayPars", signature("DisplayPars", "list"), function(x, value) {
  x <- setPar(x, value, interactive=FALSE)
  return(x)
})

setMethod("getPar", c("DisplayPars", "character"),
          function(x, name){
              if(class(x@pars)=="environment")
              {
                  name <- intersect(name, ls(x@pars, all.names=TRUE))
                  tmp <- mget(name, x@pars)
              }else{
                  name <- intersect(name, names(x@pars))
                  tmp <- x@pars[name]
              }
            tmp <- if(is.list(tmp) && length(tmp)==1) tmp[[1]] else tmp
            return(if(length(tmp)) tmp else NULL)
          })

setMethod("getPar", c("DisplayPars", "missing"), function(x) as.list(x@pars))

setMethod("displayPars", c("DisplayPars", "missing"), function(x) getPar(x))

setMethod("displayPars", c("DisplayPars", "character"), function(x, name) getPar(x, name))

setMethod("as.list", "DisplayPars", function(x) as(x, "list"))

setAs("DisplayPars", "list", function(from, to) if(!is.null(from)) as.list(from@pars) else list())

##----------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------
## InferredDisplayPars:
##
## A class to allow for querrying of available display parameters. Essentially this is a normal list with
## a bit of a fancyfied show method.
## Slots :
##    o name: the name of the class
##    o inheritance: a character vector indicating the inheritance structure
##----------------------------------------------------------------------------------------------------------------------
setClass("InferredDisplayPars", representation(name="character", inheritance="character"), contains="list")

setMethod("as.list", "InferredDisplayPars", function(x) as(x, "list"))

setAs("InferredDisplayPars", "list",
          function(from, to) {ll <- from@.Data; names(ll) <- names(from); ll})
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
##  GdObject:
##
## The GdObject is the parent of all Gviz objects in the system. 
## This class is virtual.
## Slots:
##    o dp: object of class DisplayPars to fine-tune the drawing.
##    o name: the object name, used in the title panel if needed.
##    o imageMap: an ImageMap object holding relevant image map information, or NULL
## A bunch of DisplayPars are set during object instantiation:
##    o fontface, fontcolor, fontsize, fontfamily, lineheight, cex: default settings for
##       all text unless specified someplace else.
##    o col, lwd, lty, fill: default settings for all plotting elements unless specified
##       somewhere else.
##    o col.line, col.symbol: default colors for plot symbols and plot lines. The default is
##       to take the value of the global col parameter.
##    o col.frame: the color of the panel frame, if frame==TRUE
##    o col.grid, lwd.grid, lty.grid, v and h: the default parameters for the grid plotting,
##       both when type=="g" in DataTracks and when grid==TRUE
##    o alpha: the transparancy for all track items.
##    o background.title: the fill color for the title panel. Defaults to lightgray.
##    o col.title: the font color for the title panel. Defaults to white.
##    o cex.title: the expansion factor for the title panel. This effects the fontsize
##       of both the title and the axis, if any. Defaults to NULL, which means that the
##       text size is automatically adjusted to the available space.
##    o fontfamily.title, fontface.title: the font family and font face for the title panel.
##       Defaults to bold sans serif.
##    o col.axis: the font and line color for the y axis, if any. Defaults to white.
##    o cex.axis: the expansion factor for the axis annotation. Defaults to NULL, in which case
##       it is computed based on the available space.
##    o background.panel: the background color of the content panel. Defaults to transparent.
##    o showTitle: boolean, controlling whether to plot a title panel. Although this can be 
##       set individually for each track, in multi-track plots there will still
##       be an empty placeholder in case any of the other tracks include a title. The
##       same holds true for axes. Note that the background color could be set to 
##       transparent in order to completely hide the title/axis panel.
##    o showAxis: boolean, controlling whether to plot a y axis (only track types where axes 
##	 are implemented). 
##    o grid: boolean, switching on/off the plotting of a grid.
##    o collapse: collapse the content of the track to accomodate the minimum current 
##	 device resolution
##    o min.width, min.height: the minimum width and height in pixels to display. All ranges are expanded
##       to this size in order to avoid rendering issues.
##    o min.distance: the minimum pixel distance before collapsing range items, only if collapse==TRUE
##    o frame: draw a frame around the track
##    o size: the relative size of the track
##    o ...: additional DisplayPars are allowed. Unless specified in one of the subclasses,
##       those should take the value of a valid R color descriptors. The parameter names will 
##	 later be matched to optional track item types as defined in the 'feature' range
##	 attribute, and all tracks of the matched types are colored accordingly. See the
##	 documentation of the 'GeneRegion' and 'AnnotationTrack' classes for details.
##----------------------------------------------------------------------------------------------------------------------
setClass("GdObject",
         representation=representation("VIRTUAL",
                                       dp="DisplayPars",
                                       name="character",
                                       imageMap="ImageMapOrNULL"), 
         prototype=prototype(dp=DisplayPars(fontsize=12,
                                            fontcolor="black",
                                            cex=1,
                                            fontfamily="sans",
                                            fontface=1,
                                            lineheight=1,
                                            col=Gviz:::.DEFAULT_SYMBOL_COL,
                                            col.line=NULL,
                                            col.symbol=NULL,
                                            col.grid=Gviz:::.DEFAULT_SHADED_COL,
                                            col.frame="lightgray",
                                            lwd.grid=1,
                                            lty.grid="solid",
                                            v=-1,
                                            h=-1,
                                            alpha=1,                          
                                            fill=Gviz:::.DEFAULT_FILL_COL,
                                            lwd=1,
                                            lty="solid",
                                            background.title="lightgray",
                                            col.title="white",
                                            cex.title=NULL,
                                            fontfamily.title="sans",
                                            fontface.title=2,
                                            col.axis="white",
                                            cex.axis=NULL,
                                            background.panel="transparent",
                                            showTitle=TRUE,
                                            showAxis=TRUE,
                                            grid=FALSE,
                                            collapse=TRUE,
                                            min.width=1,
                                            min.height=3,
                                            min.distance=1,
                                            frame=FALSE,
                                            size=1),
                             name="GdObject",
                             imageMap=NULL))

## We need to set and query DisplayPars in the the initializer, hence
## the appropriate methods have to be defined here first.
setMethod("setPar", signature("GdObject", "character"), function(x, name, value, interactive=TRUE) {
    newDp <- setPar(x@dp, name, value, interactive=interactive)
    x@dp <- newDp
    return(x)
})

setMethod("setPar", signature("GdObject", "list"), function(x, value, interactive=TRUE) {
    newDp <- setPar(x@dp, value, interactive=interactive)
    x@dp <- newDp
    return(x)
})

setReplaceMethod("displayPars", signature("GdObject", "list"), function(x, value) {
    x <- setPar(x, value, interactive=FALSE)
    return(x)
})

setMethod("getPar", c("GdObject", "character"), function(x, name) getPar(x@dp, name))		 

setMethod("getPar", c("GdObject", "missing"), function(x) getPar(x@dp))

setMethod("displayPars", c("GdObject", "character"), function(x, name) getPar(x, name))		 

setMethod("displayPars", c("GdObject", "missing"), function(x) getPar(x))	

## We add everything that hasn't been clobbered up so far as
## additional DisplayParameters.  Also, the dp slot must be
## re-initiated here in order to get a fresh environment for each
## instance of the class.
setMethod("initialize", "GdObject", function(.Object, name, ...) {
    ## update the default parameters first
    .makeParMapping()
    .Object <- .updatePars(.Object, "GdObject")
    ## now rebuild the slot to get a new environment
    pars <- getPar(.Object)
    .Object@dp <- DisplayPars()
    .Object <- setPar(.Object, pars, interactive=FALSE)
    if(!missing(name))
    .Object@name <- if(is.null(name)) "" else name
    ## Finally clobber up everything that's left
    .Object <- setPar(.Object, list(...), interactive=FALSE)
    return(.Object)
})
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## RangeTrack:
##
## Parent class for all region-like annotation tracks, storing information
## about start, end, strand, chromosome and the associated genome. 
## Slots:
##    o range: object of class GRangesOrIRanges containing all the necessary information
##       for plotting. The content of the elementMetadata may vary between 
##	 subclasses. The strand information may be provided in the form '+' for
##	 the Watson strand, '-' for the Crick strand or '*' for any of the two.
##    o chromosome: a character vector giving the active chromosome for which the 
##	 track is defined. Valid chromosome names are: 
##          - a single numeric character
##	    - a string, starting with 'chr', followed by any additional characters
##    o genome: character giving the reference genome for which the track is defined.
##----------------------------------------------------------------------------------------------------------------------
setClass("RangeTrack",
         representation=representation("VIRTUAL",
                                       range="GRangesOrIRanges",
                                       chromosome="character",
                                       genome="character"),
         contains="GdObject",
         prototype=prototype(range=GRanges(),
                             dp=DisplayPars(),
                             genome="ANY",
                             chromosome="chr1",
                             name="RangeTrack"))

## Coercing all input to the appropriate form			 
setMethod("initialize", "RangeTrack", function(.Object, range, chromosome, genome, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "RangeTrack")
    if(!missing(chromosome) && !is.null(chromosome))
    {
        .Object@chromosome <- .chrName(chromosome)[1]
        .Object@genome <- genome
        .Object@range <- range
    }
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## ReferenceTrack:
##
## Parent class for all tracks that provide a reference to data somewhere on the file system. This class is virtual
## and only exists for the purpose of dispatching
## Slots:
##   o stream: the import function to stream data of the disk. Needs to be able to handle the two mandatory arguments
##      'file' (a character containing a valid file path) and 'selection' (a GRanges object with the genomic region to plot)
##   o reference: the path to the file containing the data
##   o mapping: a default mapping between elementMetadata columns of the returned GRanges object from the import function
##      and the elemenMetadata columns that make up the final track object
##   o args: a list with the passed in constructor arguments during object instantiation. Those will be needed when
##     fetching the data in order to fill all necessary slots
##   o defaults: a list with the relevant default values to be used when neither 'mapping' nor 'args' provides the
##     necessary information
setClass("ReferenceTrack",  representation=representation("VIRTUAL",
                                                          stream="function",
                                                          reference="character",
                                                          mapping="list",
                                                          args="list",
                                                          defaults="list"),
         prototype=prototype(stream=function(x, selection){},
                             reference="~",
                             mapping=list()),
         validity=function(object){
             msg <- NULL
             if(!all(c("file", "selection") %in% names(formals(object@stream))))
                 msg <- "The streaming function in the 'stream' slot needs to define two arguments, 'file' and 'selection'"
             if(!file.exists(object@reference))
                 msg <- c(msg, sprintf("The referenced file '%s' does not exist", object@reference))
             return(if(is.null(msg)) TRUE else msg)
         })
setMethod("initialize", "ReferenceTrack", function(.Object, stream, reference, mapping=list(),
                                                   args=list(), defaults=list()) {
    .Object@stream <- stream
    .Object@reference <- reference
    .Object@mapping <- mapping
    .Object@args <- args
    .Object@defaults <- defaults
    validObject(.Object)
    return(.Object)
})
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## NumericTrack:
##
## Parent class for all annotation tracks that include numeric values. This is for
## dispatching purpose only.
##----------------------------------------------------------------------------------------------------------------------
setClass("NumericTrack",
         representation=representation("VIRTUAL"),
         prototype=prototype(name="NumericTrack",
                             dp=DisplayPars()),					  
         contains="RangeTrack")
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## StackedTrack
##
## Parent class for all tracks that involve stacking. This is for
## dispatching purpose only.
## Slots:
##    o stacking: character controlling the stacking of overlapping items on the final plot. 
##       One in 'hide', 'dense', 'squish', 'pack' or 'full'. 
##    o stacks: a numeric vector holding the current stacking information. Not usually set by the user. 
## This is part of the prototype only
##    o stackingValues: possible values of the stacking slot 
##	 [c("hide", "dense", "squish", "pack", "full")]
##----------------------------------------------------------------------------------------------------------------------
setClass("StackedTrack",
         representation=representation("VIRTUAL",
                                       stacking="character",
                                       stacks="numeric"),
         prototype=prototype(name="StackedTrack",
                             stacking="squish",
                             stackingValues=c("hide", "dense", "squish", "pack", "full"),
                             dp=DisplayPars(reverseStacking=FALSE,
                                            stackHeight=0.75)),					  
         contains="RangeTrack")

## Need to fill the stacks slot here, don't want to recompute all the time
setMethod("initialize", "StackedTrack", function(.Object, stacking, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "StackedTrack")
    pt <- getClass("StackedTrack")@prototype
    if(!missing(stacking))
    {
        if(!all(stacking %in% pt@stackingValues))
            stop("Problem initializing AnnotationTrack need the following values for 'stacking':",
                 paste(pt@stackingValues, collpase=", "), "\n")
        .Object@stacking <- stacking
        r <- list(...)$range
        ##stacks <- if(length(r)>0) disjointBins(ranges(r)) else 0
        ##.Object@stacks <- stacks
        .Object@stacks <- numeric()
    }
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## AnnotationTrack:
## 
## A generic Annotation track, storing information about start, end, chromosome, 
## strand and the associated genome. 
## Slots: no additional formal slots are defined, but the following are part of the prototype
##    o columns: column names that are allowed as part of the internal GRanges object in
##       the ranges slot [c("feature", "group")]
##    o featureColumnName: the column in the internal RangedData object identifying the 
##       feature type ["feature"]. If a display parameter of the same name is specified
##	 the software will use its value for the coloring.
## A bunch of DisplayPars are set during object instantiation:
##    o fill: the fill color for untyped items. Defaults to lightblue.
##    o col: the border color for all track items. This is also used to connect grouped items.
##    o lty, lwd: the line type and width for all track items. This is also used to connect grouped items.
##    o lex: the line expansion factor
##    o fontsize, fontfamily, fontface, fontcolor: the face, size, family and color for the 
##       annotation text (i.e., the track item IDs)
##    o cex: the font expansion factor.
##    o size: the relative size of the track, if not explicitely set in the plotTracks function.
##    o lineheight: the text lineheight
##    o showId: boolean controlling whether to plot track item identifiers
##    o showFeatureId: boolean controlling whether to annotate individual exons
##    o col.group, cex.group=0.7, fontface.group=2: the font color, size and face for the
##       group-level annotation
##    o shape: the shape used for the annotation items. Currently only 'box' and 'arrow' are implemented.
##    o rotation: the rotation of the item annotation in degrees.
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("AnnotationTrack",
         contains="StackedTrack",
         prototype=prototype(columns=c("feature", "group", "id"), 
                             stacking="squish",
                             name="AnnotationTrack",
                             dp=DisplayPars(fill="lightblue",
                                            col="transparent",
                                            col.line="darkgray",
                                            lty="solid",
                                            lwd=1,
                                            lex=1,
                                            fontsize=12,
                                            fontcolor="white",
                                            fontfamily="sans",
                                            fontface=1,
                                            cex.group=0.6,
                                            fontface.group=2,
                                            fontcolor.group=Gviz:::.DEFAULT_SHADED_COL,
                                            cex=1,
                                            size=1,
                                            lineheight=1,
                                            showId=FALSE,
                                            showFeatureId=FALSE,
                                            shape="arrow",
                                            rotation=0,
                                            showOverplotting=FALSE,
                                            mergeGroups=FALSE)))
                                           
## Essentially we just check for the correct GRanges columns here
setMethod("initialize", "AnnotationTrack", function(.Object, ...) {
    if(is.null(list(...)$range) && is.null(list(...)$genome) && is.null(list(...)$chromosome))
        return(.Object)
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "AnnotationTrack")
    range <- list(...)$range
    if(!is.null(range) && length(.Object))
    {
        if(!all(.Object@columns %in% colnames(values(range))))
            stop(paste("Problem initializing AnnotationTrack need the following columns:",
                     paste(.Object@columns, collpase = ", ")), "\n")
        grp <- if(is(.Object, "GeneRegionTrack")) values(range)$transcript else values(range)$group
        if(any(sapply(split(as.character(strand(range)), grp), function(x) length(unique(x))) != 1))
            stop("Grouped elments of a RangeTrack can not be on opposing strands")
    }
    .Object <- callNextMethod()
    return(.Object)
})

## The file-based version of the AnnotationTrack class. This will mainly provide a means to dispatch to
## a special 'subset' method which should stream the necessary data from disk. 
setClass("ReferenceAnnotationTrack", contains=c("AnnotationTrack", "ReferenceTrack"))

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritence has some strange features with regards to method selection
setMethod("initialize", "ReferenceAnnotationTrack", function(.Object, stream, reference, mapping=list(),
                                                             args=list(), defaults=list(), ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object=.Object, reference=reference, stream=stream,
                                                            mapping=mapping, args=args, defaults=defaults)
    .Object <- callNextMethod()
    return(.Object)
})

## Constructor. The following arguments are supported:
##    o range: a data.frame or a GRanges object containing the information
##       about the track items. If a data.frame, it needs to be coerceable
##	 to a GRanges object, i.e., it needs at least the mandatory 'start', 'stop' and 
##	 'strand' columns. Additional optional columns are:
##	    - feature: the type of the item. Can be mapped to colors via the DisplayPars.
##	    - group: a grouping factor to connect track items.
##	    - id: a unique identifier for a feature. This will be plotted if showId==TRUE
##             Note that internally we use the value of ID as the seqnames slot in the
##             internal GRanges object. Defaults for all missing columns are generated.
##	Instead of using the 'range' parameter, all these values can also be passed as
##	individual vectors, in which case they need to be of similar length.
##    o start, end, width: numeric vectors of the item start and end coordinates
##    o strand: the strand information may be provided in the form '+' for
##       the Watson strand, '-' for the Crick strand or '*' for any of the two.
##    o feature, group, id: individual vectors of equal length as described above.
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o stacking: character controlling the stacking of overlapping items. One in 'hide', 
##       'dense', 'squish', 'pack' or 'full'.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as further DisplayParameters
## (N)
AnnotationTrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, feature, group, id, strand, chromosome,
                            genome, stacking="squish", name="AnnotationTrack", fun, selectFun, importFunction,
                            stream=FALSE, ...)
{
    ## Some defaults
    covars <- .getCovars(range)
    isStream <- FALSE
    if(!is.character(range)){
        n <- max(c(length(start), length(end), length(width)), nrow(covars))
        if(is.null(covars[["feature"]]) && missing(feature))
            feature <- rep("unknown", n)
        if(is.null(covars[["id"]]) && missing(id))
            id <- make.unique(rep(if(!is.null(feature)) as.character(feature) else covars[["feature"]], n)[1:n])
        if(is.null(covars[["group"]]) && missing(group))
            group <- seq_len(n)
    }
    ## Build a GRanges object from the inputs
    .missingToNull(c("feature", "group", "id", "strand", "chromosome", "importFunction", "genome"))
    args <- list(feature=feature, group=group, id=id, strand=strand, chromosome=chromosome, genome=genome)
    defs <- list(feature="unknown", group="unknown", id="unknown", strand="*", density=1, chromosome="chrNA", genome=NA)
    range <- .buildRange(range=range, groupId="group", start=start, end=end, width=width,
                         args=args, defaults=defs, chromosome=chromosome, trackType="AnnotationTrack",
                         importFun=importFunction, stream=stream)
    if(is.list(range)){
        isStream <- TRUE
        slist <- range
        range <- GRanges()
    }
    ## Pipes have a special meaning for merged groups, so we can't have them in the initial group vector
    mcols(range)[["group"]] <- gsub("|", "", mcols(range)[["group"]], fixed=TRUE)
    ## If no chromosome was explicitely asked for we just take the first one in the GRanges object
    if(missing(chromosome) || is.null(chromosome))
        chromosome <- if(length(range)>0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    ## And finally the object instantiation, we have to distinguish between DetailsAnnotationTracks and normal ones
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    if(missing(fun))
    {
        if(!isStream){
            return(new("AnnotationTrack", chromosome=chromosome[1], range=range,
                       name=name, genome=genome, stacking=stacking, ...))
        }else{
            ## A bit hackish but for some functions we may want to know which track type we need but at the
            ## same time we do not want to enforce this as an additional argument
            e <- new.env()
            e[["._trackType"]] <- "AnnotationTrack"
            environment(slist[["stream"]]) <- e
            return(new("ReferenceAnnotationTrack", chromosome=chromosome[1], range=range,
                       name=name, genome=genome, stacking=stacking, stream=slist[["stream"]], reference=slist[["reference"]],
                       mapping=slist[["mapping"]], args=args, defaults=defs, ...))
        }
    }else{
        if(!is.function(fun))
            stop("'fun' must be a function")
        if(missing(selectFun))
            selectFun <- function(...) return(TRUE)
        if(!is.function(selectFun))
            stop("'selectFun' must be a function")
        return(new("DetailsAnnotationTrack", chromosome=chromosome[1], range=range,
                   name=name, genome=genome, stacking=stacking, fun=fun, selectFun=selectFun, ...))
    }
}
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## DetailsAnnotationTrack is an AnnotationTrack for which each fearture (annotation)
## has a detail plot (e.g. a scatter plot e.g. using lattice). The details plot
## is generated by a user defined funnction that gets called internally with
## arguments start, end and chromosome for the feature. The details plots are placed
## on top of the AnnotationTrack. As details plots can take relatively much space
## compared to features of an AnnotationTrack it only makse sense to prove details
## for an AnnotationTrack with few features. The details plots are distributed evenly
## on top of the annotation and can be connected to their corresponding featurs (they're
## providing details for) with lines for better readability.
##
## The function must have the '...' arguments. All items from the list in
## detailsFunArgs get added as arguments to the function call as well as the
## arguments, start, end, chromosome and identifier.
##
## foo = function(...) {
##    plot(densityplot(rnorm(1000), xlab=NA, ylab=NA), newpage=FALSE, prefix="foo")
## }
##
## Note, use plot with newpage=FALSE and an explicit prefix (otherwise plotting
## get more and more slow the more trellis plots you create!). See '?plot.trellis'
## for details.
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("DetailsAnnotationTrack",
         contains="AnnotationTrack",
         representation=representation(fun="function", selectFun="function"),
         prototype=prototype(fun=function(...){},
                             selectFun=function(...){return(TRUE)},
                             dp=DisplayPars(
                                            details.size=0.5,
                                            details.minWidth=100,
                                            detailsConnector.col="darkgray",
                                            detailsConnector.lty="dashed",         
                                            detailsConnector.lwd=1,
                                            detailsConnector.pch=20,          
                                            detailsConnector.cex=1,
                                            detailsBorder.lty="solid",
                                            detailsBorder.lwd=1,
                                            detailsBorder.col="darkgray",
                                            detailsBorder.fill="transparent",
                                            details.ratio=Inf,
                                            detailsFunArgs=list(),
                                            groupDetails=FALSE)
                             ))

DetailsAnnotationTrack <- function(...) AnnotationTrack(...)

setMethod("initialize", "DetailsAnnotationTrack", function(.Object, fun, selectFun, ...) {
			## the diplay parameter defaults
			.Object <- .updatePars(.Object, "DetailsAnnotationTrack")
			.makeParMapping()
			.Object@fun <- fun
                        .Object@selectFun <- selectFun
			.Object <- callNextMethod()
			return(.Object)
		})
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## GeneRegionTrack:
##
## A track containing all genes in a particular region. The data are usually fetched dynamially
## from an online data store. Of course it would
## also be possible to manully construct an object from local data. Particular data sources
## should be implemented as sub-classes, this is just the commone denominator that is being
## used for plotting later on. There are several levels of data associated to a GeneRegionTrack:
##   - exon level: identifiers are stored in the exon column of the GRanges object. Data may be extracted
##          using the 'exons' method.
##   - transcript level: identifiers are stored in the transcript column of the GRanges object.
##          Data may be extracted using the 'transcripts' method. 
##   - gene level: identifiers are stored in the gene column of the GRanges object, more human-readable
##          versions in the symbol column. Data may be extracted using the 'genes' or the 'symbols' methods.
##   - transcript-type level: information is stored in the feature column of the GRanges object.
##      If a display parameter of the same name is specified the software will use its value for the coloring.
## Slots: no additional formal slots are defined, but the following are part of the prototype
##    o columns: column names that are allowed as part of the internal GRanges object in
##       the ranges slot [c("feature", "transcript", "symbol", "gene")]
## A bunch of DisplayPars are set during object instantiation:
##    o fill: the fill color for untyped items. Defaults to lightblue.
##    o col: the border color for all track items. This is also used to connect grouped items.
##    o alpha: the transparancy for all track items.
##    o lty, lwd: the line type and width for all track items. This is also used to connect grouped items.
##    o lex: the line expansion factor
##    o fontsize, fontfamily, fontface, fontcolor: the face, size, family and color for the 
##       annotation text (i.e., the track item IDs)
##    o cex: the font expansion factor.
##    o lineheight: the text lineheight
##    o showId: boolean controlling whether to plot group identifiers, e.g, gene symbols or transcript IDs.
##    o geneSymbols: boolean controlling whether to use gene symbols or gene ids as gene annotations
##    o shape: the shape used for the annotation items. Currently only 'box' and 'arrow' are implemented.
##    o rotation: the rotation of the exon annotation in degrees.
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("GeneRegionTrack",
         contains="AnnotationTrack",
         representation=representation(start="numeric", end="numeric"),
         prototype=prototype(columns=c("feature", "transcript", "symbol", "gene", "exon"),
                             stacking="squish",
                             stacks=0,
                             start=0,
                             end=0,
                             name="GeneRegionTrack",
                             dp=DisplayPars(fill="orange",
                                            min.distance=0,
                                            col=NULL,
                                            geneSymbols=TRUE,
                                            showExonId=FALSE,
                                            collapseTranscripts=FALSE,
                                            shape=c("smallArrow", "box"),
                                            thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "3UTR", "5UTR", "miRNA", "lincRNA",
                                                             "three_prime_UTR", "five_prime_UTR"))))
                                           

## Making sure all the display parameter defaults are being set
setMethod("initialize", "GeneRegionTrack", function(.Object, start, end, ...){
    if(is.null(list(...)$range) && is.null(list(...)$genome) && is.null(list(...)$chromosome))
        return(.Object)
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "GeneRegionTrack")
    .Object@start <- ifelse(is.null(start), 0 , start)
    .Object@end <- ifelse(is.null(end), 0 , end)
    .Object <- callNextMethod()
    return(.Object)
})

## The file-based version of the GeneRegionTrack class. This will mainly provide a means to dispatch to
## a special 'subset' method which should stream the necessary data from disk. 
setClass("ReferenceGeneRegionTrack", contains=c("GeneRegionTrack", "ReferenceTrack"))

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritence has some strange features with regards to method selection
setMethod("initialize", "ReferenceGeneRegionTrack", function(.Object, stream, reference, mapping=list(),
                                                             args=list(), defaults=list(), ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object=.Object, reference=reference, stream=stream,
                                                            mapping=mapping, args=args, defaults=defaults)
    .Object <- callNextMethod()
    return(.Object)
})

## Constructor. The following arguments are supported:
##    o range: one in a whole number of differerent potential inputs upon which the .buildRanges method will dispatch
##    o start, end: numeric vectors of the track start and end coordinates. 
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o rstarts, rends, rwidths: integer vectors of exon start and end locations or widths, or a character vector of
##      comma-delimited exon locations, one vector element for each transcript
##    o strand, feature, exon, transcript, gene, symbol, chromosome: vectors of equal length containing
##       the exon strand, biotype, exon id, transcript id, gene id, human-readable gene symboland chromosome information
##    o stacking: character controlling the stacking of overlapping items. One in 'hide', 
##       'dense', 'squish', 'pack' or 'full'.
##    o name: the name of the track. This will be used for the title panel.
##    o exonList: boolean, causing the values in starts, rends or rwidths to be interpreted as delim-deparated
##       lists that have to be exploded. All other annotation arguments will be repeated accordingly.
##    o delim: the delimiter if coordinates are in a list
## All additional items in ... are being treated as DisplayParameters
## (N)
GeneRegionTrack <- function(range=NULL, rstarts=NULL, rends=NULL, rwidths=NULL, strand, feature, exon,
                            transcript, gene, symbol, chromosome, genome, stacking="squish",
                            name="GeneRegionTrack", start=NULL, end=NULL, importFunction, stream=FALSE, ...)
{
    ## Some defaults
    covars <- if(is.data.frame(range)) range else if(is(range, "GRanges")) as.data.frame(mcols(range)) else data.frame()
    isStream <- FALSE
    if(!is.character(range)){
        n <- if(is.null(range)) max(c(length(start), length(end), length(width))) else if(is(range, "data.frame")) nrow(range) else length(range)
        if(is.null(covars[["feature"]]) && missing(feature))
            feature <- rep("unknown", n)
        if(is.null(covars[["exon"]]) && missing(exon))
            exon <- make.unique(rep(if(!missing(feature) && !is.null(feature)) as.character(feature) else covars[["feature"]], n)[1:n])
        if(is.null(covars[["transcript"]]) && missing(transcript))
            transcript <- paste("transcript", seq_len(n), sep="_")
        if(is.null(covars[["gene"]]) && missing(gene))
            gene <- paste("gene", seq_len(n), sep="_")
    }
    ## Build a GRanges object from the inputs
    .missingToNull(c("feature", "exon", "transcript", "gene", "symbol", "strand", "chromosome", "importFunction", "genome"))
    args=list(feature=feature, id=exon, exon=exon, transcript=transcript, gene=gene, symbol=symbol, strand=strand,
              chromosome=chromosome, genome=genome)
    defs <- list(feature="unknown", id="unknown", exon="unknown", transcript="unknown", genome=NA,
                 gene="unknown", symbol="unknown", strand="*", density=1, chromosome="chrNA")
    range <- .buildRange(range=range, groupId="transcript", start=rstarts, end=rends, width=rwidths, args=args, defaults=defs,
                         chromosome=chromosome, tstart=start, tend=end, trackType="GeneRegionTrack", importFun=importFunction,
                         genome=genome)
    if(is.list(range)){
        isStream <- TRUE
        slist <- range
        range <- GRanges()
    }
    if(is.null(start))
        start <- if(!length(range)) NULL else min(IRanges::start(range))
    if(is.null(end))
        end <- if(!length(range)) NULL else max(IRanges::end(range))
    if(missing(chromosome) || is.null(chromosome))
        chromosome <- if(length(range)>0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    if(!isStream){
        return(new("GeneRegionTrack", start=start, end=end, chromosome=chromosome[1], range=range,
                   name=name, genome=genome, stacking=stacking, ...))
    }else{
        return(new("ReferenceGeneRegionTrack", start=start, end=end, chromosome=chromosome[1], range=range,
                   name=name, genome=genome, stacking=stacking, stream=slist[["stream"]],
                   reference=slist[["reference"]], mapping=slist[["mapping"]], args=args, defaults=defs, ...))
    }
}
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## BiomartGeneRegionTrack:
##
## A GeneRegionTrack that fetches its information from Biomart
## Slots:
##    o biomart: a biomaRt object providing the connection to the data source.
##        currently we assume that this connects to a esembl_gene source for
##        a particular organism.
##    o filter: a named list of additional filters that are passed on to the Biomart query.
## This is defined only in the prototype:
##    o columns: column names that are allowed as part of the internal GRanges object in
##        the ranges slot [c("feature", "transcript", "symbol", "gene")]
## A bunch of DisplayPars are set during object instantiation, all of which are
## transcript types as returned from Biomart in the 'biotype' field.
## This class mainly exists for dispatching purpose and to keep things as flexible
## as possible, the actual plottable information is all contained in the parent
## class definition.
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("BiomartGeneRegionTrack",
         contains="GeneRegionTrack",
         representation=representation(biomart="MartOrNULL", filter="list"),
         prototype=prototype(biomart=NULL,
                             filters=list(),
                             columns=c("feature", "transcript", "symbol", "gene", "rank"),
                             name="BiomartGeneRegionTrack",
                             dp=DisplayPars(C_segment="burlywood4",
                                            D_segment="lightblue",
                                            J_segment="dodgerblue2",
                                            miRNA="cornflowerblue",
                                            miRNA_pseudogene="cornsilk",
                                            misc_RNA="cornsilk3",
                                            misc_RNA_pseudogene="cornsilk4",
                                            Mt_rRNA="yellow",
                                            Mt_tRNA="darkgoldenrod",
                                            Mt_tRNA_pseudogene="darkgoldenrod1",
                                            protein_coding="orange",
                                            utr5="orange",
                                            utr3="orange",
                                            pseudogene="brown1",         
                                            retrotransposed="blueviolet",
                                            rRNA="darkolivegreen1",
                                            rRNA_pseudogene="darkolivegreen" ,   
                                            scRNA="gold4",
                                            scRNA_pseudogene="darkorange2",
                                            snoRNA="cyan",           
                                            snoRNA_pseudogene="cyan2",
                                            snRNA="coral",
                                            snRNA_pseudogene="coral3",   
                                            tRNA_pseudogene="antiquewhite3",
                                            V_segment="aquamarine")))

## Retrieving information from Biomart.
setMethod("initialize", "BiomartGeneRegionTrack", function(.Object, start, end, biomart, filters=list(), range, genome, chromosome, strand, ...){
    if((missing(range) || is.null(range)) && is.null(genome) && is.null(chromosome))
        return(.Object)
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "BiomartGeneRegionTrack")
    if(missing(start) || missing(end))
    {
        .Object <- setPar(.Object, "size", 0, interactive=FALSE)
        return(.Object)
    }
    ## changing start and end positions to capture genes on the edges.
    sstart <- start - 2000
    send <- end + 2000
    ## fetching data from Biomart
    .Object@biomart <- biomart
    if (!is.null(.Object@biomart))
    {
        attributes <- c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start",
                        "exon_chrom_end", "rank", "strand", "external_gene_id", "gene_biotype", "chromosome_name",
                        "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end", "phase")
        filterNames <- c("chromosome_name", "start", "end", names(filters))
        filterValues <- c(list(gsub("^chr", "", chromosome), sstart, send), as.list(filters))
        strand <- .strandName(strand, extended=TRUE)
        if(strand %in% 0:1)
        {
            filterNames <- c(filterNames, "strand")
            filterValues <- c(filterValues, c(1, -1)[strand+1])
        }
        ens <- getBM(attributes, filters=filterNames,                           
                     values=filterValues,
                     mart=.Object@biomart, uniqueRows=TRUE)
        colnames(ens) <- c("gene_id","transcript_id","exon_id","start",
                           "end", "rank", "strand", "symbol", "biotype",
                           "chromosome_name", "u5s", "u5e", "u3s", "u3e", "phase")
        ## We may have to split exons if they contain UTRs
        hasUtr <- !is.na(ens$u5s) | !is.na(ens$u3s)
        ensUtr <- ens[hasUtr,, drop=FALSE]
        ensUtr$feature <- ifelse(is.na(ensUtr$u5s), "utr3", "utr5")
        ensUtr$us <- ifelse(ensUtr$feature=="utr3", ensUtr$u3s, ensUtr$u5s)
        ensUtr$ue <- ifelse(ensUtr$feature=="utr3", ensUtr$u3e, ensUtr$u5e)
        ensUtr$u5e <- ensUtr$u5s <- ensUtr$u3e <- ensUtr$u3s <- NULL
        allUtr <- ensUtr$us == ensUtr$start & ensUtr$ue == ensUtr$end
        utrFinal <- ensUtr[allUtr,, drop=FALSE]
        ensUtr <- ensUtr[!allUtr,, drop=FALSE]
        ensUtrS <- split(ensUtr, ifelse(ensUtr$start==ensUtr$us, "left", "right"))
        utrFinal <- rbind(utrFinal, do.call(rbind, lapply(names(ensUtrS), function(i){
            y <- ensUtrS[[i]]
            if(nrow(y)==0)
                return(NULL)
            yy <- y[rep(1:nrow(y), each=2),]
            sel <- seq(1, nrow(yy), by=2)
            yy[sel, "end"] <- if(i=="left") yy[sel, "ue"] else yy[sel, "us"]-1
            yy[sel, "feature"] <-  yy[sel, ifelse(i=="left", "feature", "biotype")]
            yy[sel, "phase"] <-  if(i=="left") -1 else 0
            sel <- seq(2, nrow(yy), by=2)
            yy[sel, "start"] <- if(i=="left") yy[sel, "ue"]+1 else yy[sel, "us"]
            yy[sel, "feature"] <-  yy[sel, ifelse(i=="left", "biotype", "feature")]
            yy[sel, "phase"] <- if(i=="left") yy[sel, "phase"] else -1
            yy
        })))
        utrFinal$biotype <- utrFinal$feature
        keep <-  c("gene_id","transcript_id","exon_id","start",
                   "end", "rank", "strand", "symbol", "biotype",
                   "chromosome_name", "phase")
        ens <- rbind(ens[!hasUtr,keep, drop=FALSE], utrFinal[,keep])
        prefix <- ifelse(!grepl("^chr", chromosome) && getOption("ucscChromosomeNames"), "chr", "")
        range <- GRanges(seqnames=paste(prefix, .chrName(ens$chromosome_name), sep="")[seq_len(nrow(ens))],
                         ranges=IRanges(start=ens$start, end=ens$end),
                         strand=ens$strand, feature=as.character(ens$biotype),
                         gene=as.character(ens$gene_id), exon=as.character(ens$exon_id),
                         transcript=as.character(ens$transcript_id), symbol=as.character(ens$symbol),
                         rank=as.numeric(ens$rank), phase=as.integer(ens$phase))
        suppressWarnings(genome(range) <- unname(genome[1]))
        range <- sort(range)
    }
    if(length(range)==0)
        .Object <- setPar(.Object, "size", 0, interactive=FALSE)
    .Object <- callNextMethod(.Object=.Object, range=range, start=start, end=end,
                              genome=genome, chromosome=chromosome, strand=strand, ...)
    return(.Object)
})

## Constructor. The following arguments are supported:
##    o start, end: numeric vectors of the item start and end coordinates
##    o biomart: a biomaRt object used to to query for gene annotations
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o strand: character, search for gene models on the plus strand ("+" or 0), the
##       minus strand ("-" or 1) or both strands ("+-" or "-+" or 2)
##    o stacking: character controlling the stacking of overlapping items. One in 'hide', 
##       'dense', 'squish', 'pack' or 'full'.
##    o filters: list of additional filters for the biomaRt query, where the item names
##     	 are the filter identifiers and the item values are the filter values.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters
## (N)
BiomartGeneRegionTrack <- function(start, end, biomart, chromosome, strand, genome,
                                   stacking="squish", filters=list(), name="BiomartGeneRegionTrack", ...)
{
    ## Some default checking
    .missingToNull(c("genome"))
    if(missing(strand))
        strand <- "*"
    if(missing(start) || missing(end))
        stop("Need to specify a start and end for creating a BiomartGeneRegionTrack.")
    if(missing(chromosome))
        stop("A chromosome must be specified for this annotation track.")
    chromosome <- .chrName(chromosome)[1]
    if(is.null(genome))
        stop("A genome must be specified for this annotation track.")
    if(missing(biomart))
    {
        if(is.null(genome))
            stop("Need either a valid BiomaRt connection object or a UCSC genome identifier as the 'genome' argument.")
        biomart <- .genome2Dataset(genome)
    }
    new("BiomartGeneRegionTrack", start=start, end=end, chromosome=chromosome, strand=strand,
        biomart=biomart, name=name, genome=genome, stacking=stacking, filters=filters, ...)
}


## This filters for genes with refseq IDs only
.refseqFilter <- list("with_refseq_dna"=TRUE)
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## GenomeAxisTrack:
##
## A track for the genome axis
## Slots:
##    o range: an object of class GRanges containing ranges to be highlighted along the axis.
## A bunch of DisplayPars are set during object instantiation:
##    o col: the color for the axis text, lines and tickmarks
##    o fill.range: the fill color for the optional range annotation
##    o col.range: the border color for the optional range annotation
##    o fontsize: the font size for the axis annotation
##    o showTitle: boolean, show or hide the title panel text
##    o background.title: the background color of the title panel
##    o cex: the cex value for the axis text
##    o exponent: the exponent for the axis coordinates. E.g., 3 means mb, 6 means gb, etc.
##    o distFromAxis: numeric, distance of text from the axis
##    o labelPos: character giving the position of axis labels, one in "alternating", "revAlternating",
##       "above" or "below"
##    o add53: boolean, add a 5'->3' indicator
##    o add35: boolean, add a 3'->5' indicator
##    o littleTicks: boolean, add second level of smaller tick marks
##    o size: the relative size of the track
##    o col.id, cex.id: text settings for the optional range annotation
##    o showId: boolean, show range annotation
##    o scale: numeric, if not NULL a small scale is drawn instead of the full axis,
##      if between 0 and 1 it is interpreted as a fraction of the region (otherwise absolute).
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("GenomeAxisTrack",
         contains="GdObject", 
         representation=representation(range="GRanges"),
         prototype(range=GRanges(),
                   name="GenomeAxisTrack",
                   dp=DisplayPars(col="darkgray",
                                  fontcolor="#808080",
                                  fill.range="cornsilk3",
                                  col.range="cornsilk4",
                                  fontsize=10,
                                  showTitle=FALSE,
                                  background.title="transparent",
                                  cex=0.8,
                                  exponent=NULL,
                                  distFromAxis=1,
                                  labelPos="alternating",
                                  add53=FALSE,
                                  add35=FALSE,
                                  col.id="white",
                                  cex.id=0.7,
                                  showId=FALSE,
                                  littleTicks=FALSE,
                                  size=NULL,
                                  scale=NULL,
                                  lwd=2)))

## Only pass on the stuff to the GdObject initializer
setMethod("initialize", "GenomeAxisTrack", function(.Object, range, ids, ...){
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "GenomeAxisTrack")
    if(missing(range) || is.null(range))
        range <- GRanges()
    if(is(range, "IRanges"))
        range <- GRanges(range=range, seqnames="dummy", id=ids)
    .Object@range <- range
    .Object<- callNextMethod()
    return(.Object)
})

## Constructor. The following arguments are supported:
##    o range: an object of class 'GRanges' containing regions to be highlighted on the axis by colored boxes
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters
## (N)
GenomeAxisTrack <- function(range=NULL, name="Axis", id, ...)
{
    if(missing(id))
        id <- names(range)
    new("GenomeAxisTrack", name=name, range=range, id=id, ...)
}
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## DataTrack:
##
## A track for all sorts of quantitative values
## Slots:
##    o data: a numeric matrix containing the data values
## A bunch of DisplayPars are set during object instantiation:
##    o jitter.x, jitter.y, amount, factor: parameters controlling the jittering in xy-type plots.
##    o span, degree, family, evaluation: parameters controlling the loess calculation.
##    o col.mountain, fill.mountain, col.baseline, lwd.baseline, lty.baseline, baseline:
##       parameters controlling the colors and the baseline position in mountain-type plots.
##    o fill.histogram, col.histogram: fill and line colors for the histogram-type plots
##    o box.ratio, box.width, varwidth, notch, notch.frac, levels.fos, stats, coef, do.out:
##       parameters controlling the boxplot appearance
##    o size: the relative size of the track
##    o type: the plot type, one or several in c("p", "l", "b", "a", "s", "g", "r", "S", "smooth",
##       "histogram", "mountain", "h", "boxplot", "gradient", "heatmap", "polygon")
##    o cex: the default pixel size
##    o ncolor, gradient: the number of colors and the base colors for the gradient type
##    o collpase: collapse overlapping ranges
##    o min.distance: the mimimum distance in pixel below which to collapse
##    o window: average the rows of the data matrix to 'window' slices on the data range
##    o separator: number of pixels used to separate individual samples in heatmap-type plots
##    o transformation: a function applied on the data matrix prior to plotting. The function
##         should accept exactly one argument and the return value needs to be a numeric vector
##         which can be coerced back into a data matrix of identical dimensionality as the input
##         data.
##    o aggregation: a function to aggregate values in windows or for collapsed items. Either a
##      function that collapses a numeric vector into a single number, or one of the predefined
##      options "mean", "median", "sum" "min", "max" or "extreme". Defaults to "mean"
##    o stackedBars: logical, draw stacked histograms if groups!=NULL. Else show grouped data as
##         side-by-side bars.
##    o na.rm: remove NA values before plotting
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("DataTrack",
         contains="NumericTrack",
         representation=representation(data="matrix", strand="character"),
         prototype=prototype(columns=c("score"),
                             name="DataTrack",
                             dp=DisplayPars(size=NULL,
                                            type="p",
                                            cex=0.7,
                                            jitter.x=FALSE,
                                            jitter.y=FALSE,
                                            factor=0.5,
                                            amount=NULL,
                                            span=1/5, 
                                            degree=1,
                                            pch=20,
                                            family="symmetric",
                                            evaluation=50,
                                            col=trellis.par.get("superpose.line")[["col"]],
                                            col.mountain=NULL,
                                            lwd.mountain=NULL,
                                            lty.mountain=NULL,
                                            fill.mountain=c("#CCFFFF", "#FFCCFF"),
                                            baseline=NULL,
                                            col.baseline=NULL,
                                            lwd.baseline=NULL,
                                            lty.baseline=NULL,
                                            box.ratio=1,
                                            box.width=NULL,
                                            varwidth=FALSE,
                                            notch=FALSE,
                                            notch.frac=0.5,
                                            levels.fos=NULL,
                                            stats=boxplot.stats,
                                            coef=1.5,
                                            do.out=TRUE,
                                            transformation=NULL,
                                            ncolor=100,
                                            gradient=brewer.pal(9, "Blues"),
                                            min.distance=0,
                                            collapse=FALSE,
                                            window=NULL,
                                            windowSize=NULL,
                                            fill.histogram=NULL,
                                            col.histogram=Gviz:::.DEFAULT_SHADED_COL,
                                            stackedBars=TRUE,
                                            groups=NULL,
                                            separator=0,
                                            aggregation="mean",
                                            aggregateGroups=FALSE,
                                            ylim=NULL,
                                            na.rm=FALSE,
                                            showSampleNames=FALSE,
                                            cex.sampleNames=NULL,
                                            col.sampleNames="white")))

## Only pass on the stuff to the GdObject initializer
setMethod("initialize", "DataTrack", function(.Object, data=matrix(), strand, ...){
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "DataTrack")
    .Object@data <- data
    if(!missing(strand))
        .Object@strand <- unique(strand)
    .Object <- callNextMethod()
    return(.Object)
})

## The file-based version of the DataTrack class. This will mainly provide a means to dispatch to
## a special 'subset' method which should stream the necessary data from disk. 
setClass("ReferenceDataTrack", contains=c("DataTrack", "ReferenceTrack"))

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritence has some strange features with regards to method selection
setMethod("initialize", "ReferenceDataTrack", function(.Object, stream, reference, mapping=list(),
                                                       args=list(), defaults=list(), ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object=.Object, reference=reference, stream=stream,
                                                            mapping=mapping, args=args, defaults=defaults)
    .Object <- callNextMethod()
    return(.Object)
})

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
## (N)
DataTrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, data, chromosome, strand, genome,
                      name="DataTrack", importFunction, stream=FALSE, ...)
{
    ## Build a GRanges object from the inputs
    wasGR <- is(range, "GRanges") || is.character(range)
    fromFile <- is.character(range)
    isStream <- FALSE
    .missingToNull(c("strand", "chromosome", "importFunction", "genome"))
    args <- list(strand=strand, chromosome=chromosome, genome=genome)
    defs <- list(strand="*", chromosome="chrNA")
    range <- .buildRange(range=range, start=start, end=end, width=width, args=args, defaults=defs,
                         asIRanges=FALSE, chromosome=chromosome, genome=NA, trackType="DataTrack",
                         importFun=importFunction, stream=stream)
    if(is.list(range)){
        isStream <- TRUE
        slist <- range
        range <- GRanges()
    }
    ## Some default checking
    if(length(unique(strand(range)))>1)
        stop("The strand has to be unique for all ranges in a DataTrack object.")
    if(!missing(data))
    {
        if(is.character(data))
        {
            if(!wasGR)
                stop("Columns indices for the data section are only allowed when 'range' is of class 'GRanges'")
            mt <- is.na(match(data, colnames(values(range))))
            if(any(mt))
                warning("Unable to match data columns: ", paste(data[mt], collapse=","))
            data <- as.data.frame(values(range)[,data[!mt], drop=FALSE])
        }
        if(is.null(dim(data)))
            dim(data) <- c(1, length(data))
        if(is.matrix(data))
            data <- as.data.frame(t(data))
    } else {
        data <- if(ncol(values(range))) as.data.frame(values(range)) else matrix(nrow=0, ncol=0)
    }
    data <- .prepareDtData(data, len=length(range))
    if(missing(chromosome) || is.null(chromosome))
        chromosome <- if(length(range)>0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    genome <- .getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    values(range) <- NULL
    if(!isStream){
        return(new("DataTrack", chromosome=chromosome, strand=as.character(strand(range)), range=range,
                   name=name, genome=genome, data=data, ...))
    }else{
        ## A bit hackish but for some functions we may want to know which track type we need but at the
        ## same time we do not want to enforce this as an additional argument
        e <- new.env()
        e[["._trackType"]] <- "DataTrack"
        environment(slist[["stream"]]) <- e
        return(new("ReferenceDataTrack", chromosome=chromosome, strand=as.character(strand(range)), range=range,
                   name=name, genome=genome, data=data, stream=slist[["stream"]], reference=slist[["reference"]],
                   mapping=slist[["mapping"]], args=args, defaults=defs, ...))
    }
}
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## IdeogramTrack:
##
## A track for chromosome ideograms
## A bunch of DisplayPars are set during object instantiation:
##    o size: the relative size of the track
##    o col, fill: the border and fill color used for the highlighted currently displayed region on the chromosome
##    o fontcolor, fontface, fontfamily, fonsize, lineheight, cex : the color, family, face, size, lineheight and
##       expansion factor for the chromosome name text
##    o showId: indicate the chromosome name next to the ideogram
##    o showBandId: show the identifiers of the chromosome bands
##    o cex.bands: character expansion factor for the chromosome band information
##    o bevel: the amount of beveling at the ends of the ideogram. A number between 0 and 1.
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("IdeogramTrack", contains = "RangeTrack",
         representation=representation(bandTable="data.frame"),
         prototype=prototype(name="IdeogramTrack",
                             bandTable=data.frame(),
                             dp=DisplayPars(size=NULL,
                                            col="red",
                                            fill="#FFE3E6",
                                            fontcolor="#808080",
                                            cex=0.8,
                                            fontsize=10,
                                            showTitle=FALSE,
                                            background.title="transparent",
                                            showId=TRUE,
                                            showBandId=FALSE,
                                            cex.bands=0.7,
                                            bevel=0.45)))

## Grab the chromosome band and length information from UCSC and fill the ranges slot.
setMethod("initialize", "IdeogramTrack", function(.Object, genome, chromosome, bands, name, ...){
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "IdeogramTrack")
    if(is.null(bands) && (missing(genome) || missing(chromosome)))
        return(callNextMethod(.Object=.Object, range=GRanges(), genome=NULL, chromosome=NULL, ...))
    if(is.null(bands)){
        sessionInfo <- .cacheGenomes(genome=genome)
        .Object@bandTable <- sessionInfo$bands
        bands <- sessionInfo$bands
    }else{
        .checkClass(bands, "data.frame")
        cols <- c("chrom", "chromStart", "chromEnd", "name", "gieStain") 
        miss <- ! cols %in% colnames(bands)
        if(any(miss))
            stop(sprintf("The following column%s missing from the bands table: %s",
                         ifelse(sum(miss)>1, "s are", " is"), paste(cols[miss], collapse=", ")))
        .Object@bandTable <- bands
    }
    chromosome <- if(is.null(chromosome)) as.character(bands[1, "chrom"]) else .chrName(chromosome)[1]
    bands <- bands[bands$chrom==chromosome,]
    if(nrow(bands)==0)
        stop("Chromosome '", chromosome, "' does not exist on UCSC genome '", genome, "'")
    if(is.null(name))
        name <- .chrName(chromosome)[1]
    ranges <- GRanges(seqnames=as.character(bands$name), range=IRanges(start=bands$chromStart, end=bands$chromEnd),
                      name=as.character(bands$name), type=as.character(bands$gieStain))
    .Object <- callNextMethod(.Object=.Object, range=ranges, genome=genome, chromosome=chromosome, name=name, ...)
    return(.Object)
})

## Constructor. The following arguments are supported:
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters
## (N)
IdeogramTrack <- function(chromosome=NULL, genome, name=NULL, bands=NULL, ...){
    if(missing(genome)) stop("Need to specify genome for creating an IdeogramTrack")
    new("IdeogramTrack", chromosome=chromosome, genome=genome, name=name, bands=bands, ...)
}
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## UcscTrack:
##
## Strictly speaking this is not a class, but rather some sort of meta-constructor for several of the previously
## defined track types directly from UCSC data. It will fetch online data from a particular track (or a sub-table of
## a track) and feed it to one of the original constructors with user-provided argument mappings.
##----------------------------------------------------------------------------------------------------------------------
## Constructor. The following arguments are supported:
##    o track: a character of one of the available UCSC tracks
##    o table: a character of one one of the sub-tables of the track, or NULL to fetch all
##    o trackType: a character giving the name of the constructor to pass the data to, one in
##       c("AnnotationTrack", "GeneRegionTrack", "DataTrack", "GenomeAxisTrack")
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o from: the starting and end coordinates of the track data
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters

## A simple caching mechanism for UCSC session information. The overhead for establishing a connection to UCSC is
## quite significant and we can shave off 5 to 10 seconds here by caching sessions and associated information
## for a particular genome and chromosome.
.ucscCache <- new.env()
.ensemblCache <- new.env()
.doCache <- function(token, expression, env, callEnv=environment())
{
     if(!token %in% ls(env))
     {
         res <- eval(expression, envir=callEnv)
         assign(x=token, value=res, envir=env)
         res
     } else env[[token]]
}
.cacheTracks <- function(genome, chromosome, track, env=.ucscCache)
{
    genomes <- .doCache("availableGenomes", expression(ucscGenomes()), env)
    if(!genome %in% as.character(genomes[,"db"]))
        stop("'", genome, "' is not a valid UCSC genome.")
    sessionToken <- paste("session", genome, sep="_")
    tracksToken <- paste("tracks", genome, sep="_")
    tablesToken <- paste("tables", track, genome, sep="_")
    cenv <- environment()
    session <- .doCache(sessionToken,
                        expression({tmp <- browserSession()
                                    genome(tmp) <- genome
                                    tmp}), env, cenv)
    availTracks <- .doCache(tracksToken, expression(trackNames(session)), env, cenv)
    track <- match.arg(track, sort(c(availTracks, names(availTracks))))
    if(!is.na(availTracks[track]))
        track <- availTracks[track]
    availTables <- .doCache(tablesToken, expression({query <- ucscTableQuery(session, track)
                                                     sort(tableNames(query))}), env, cenv)
    chrInfo <- seqlengths(session)
    return(list(session=session, availTracks=availTracks, availTables=availTables, track=track,
                chrInfo=chrInfo))
}
.cacheGenomes <- function(genome=NULL, env=.ucscCache)
{
    availToken <- "availableGenomes"
    genomesToken <- paste("genomeBands", genome, sep="_")
    genomes <- .doCache(availToken, expression(ucscGenomes()), env)
    bands <- NULL
    if(!is.null(genome))
    {
        cenv <- environment()
        bands <- .doCache(genomesToken, expression({
            if(!genome %in% as.character(genomes[,"db"]))
                stop("'", genome, "' is not a valid UCSC genome.")
            sessionToken <- paste("session", genome, sep="_")
            session <- .doCache(sessionToken,
                                expression({tmp <- browserSession()
                                            genome(tmp) <- genome
                                            tmp}), env, cenv)
            query <-  tryCatch(ucscTableQuery(session, "cytoBandIdeo"), error=function(e)
                               stop("There doesn't seem to be any cytoband data available for genome '", genome, "' at UCSC."))
            getTable(query)}), env, cenv)
    }
    return(list(availableGenomes=genomes, bands=bands))
}

## empty the session cache
clearSessionCache <- function(){
    assignInNamespace(".ucscCache", new.env(), ns="Gviz")
    assignInNamespace(".ensemblCache", new.env(), ns="Gviz")
}

    

## (N)
UcscTrack <- function(track, table=NULL, trackType=c("AnnotationTrack", "GeneRegionTrack",
                                           "DataTrack", "GenomeAxisTrack"),
                      genome, chromosome, name=NULL, from, to, ...)
{
    trackType <- match.arg(trackType)
    if(missing(genome) || !IRanges:::isSingleString(genome)) stop("Need to specify genome for creating a UcscTrack")
    if(missing(chromosome)) stop("Need to specify chromosome for creating a UcscTrack")
    chromosome <- .chrName(chromosome)[1]
    sessionInfo <- .cacheTracks(genome=genome, chromosome=chromosome, track=track, env=.ucscCache)
    if(missing(from))
        from <- 1
    if(missing(to))
        to <- sessionInfo$chrInfo[chromosome]
    gr <- GRanges(ranges=IRanges(start=from, end=to), seqnames=chromosome)
    suppressWarnings(genome(gr) <- unname(genome))[1]
    query <- ucscTableQuery(sessionInfo$session, sessionInfo$track, gr)
    if(!is.null(table))
    {
        table <- match.arg(table, sessionInfo$availTables)
        tableName(query) <- table
    }
    if(is.null(name))
      name <- if(is.null(table)) track else paste(sessionInfo$track, table)
    tableDat <- if(trackType=="DataTrack"){
        tmp <- try(track(query, asRangedData=FALSE), silent=TRUE)
        if(is(tmp, "try-error")){
            warning(tmp)
            data.frame()
        } else as.data.frame(tmp)} else {
            tmp <- try(getTable(query), silent=TRUE)
            if(is(tmp, "try-error")){
                warning(tmp)
                data.frame()
            } else tmp}
    if(is(tmp, "try-error") && nrow(tableDat)==0)
        stop("Error fetching data from UCSC")
    args <- lapply(list(...), function(x) if(is.character(x) && length(x)==1)
                   if(!x %in% colnames(tableDat)) x else tableDat[,x] else x)
    if(trackType=="GeneRegionTrack")
    {
        args$start <- from
        args$end <- to
    }
    args <- lapply(args, function(x) if(!length(x))  NULL else x)
    trackObject <- do.call(trackType, args=c(list(chromosome=chromosome, genome=genome, name=name), args))
    return(trackObject)
  }
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## AlignedReadTrack:
## 
## A track to visualize sequence reads, as typically produced by NGS experiments
## Slots: no additional formal slots are defined
## A bunch of DisplayPars are set during object instantiation:
##    o detail: the amount of plotting details to show for the aligned reads, one in c("reads", "coverage")
##----------------------------------------------------------------------------------------------------------------------
## (N)
setClass("AlignedReadTrack",
         representation=representation(coverage="list",
                                       coverageOnly="logical"),
         contains="StackedTrack",
         prototype=prototype(stacking="squish",
                             name="AlignedReadTrack",
                             coverageOnly=FALSE,
                             dp=DisplayPars(detail="coverage",
                                            type="histogram",
                                            fill="#0080ff",
                                            size=NULL,
                                            collapse=FALSE)))

## Recompute coverage on plus and minus strand and combined strands for AlignedRead tracks
## and update the respective slot to hold this information
setMethod("setCoverage", signature("AlignedReadTrack"),
          definition=function(GdObject){
              if(length(GdObject))
              {
                  str <- factor(strand(GdObject), levels=c("+", "-"))
                  gdSplit <- as.list(split(range(GdObject), str))
                  covs <- lapply(gdSplit, coverage)
                  covs[["*"]] <- coverage(range(GdObject))
                  GdObject@coverage <- covs
              }
              return(GdObject)
          })

## Essentially we just update the display parameters here and precompute the coverage
setMethod("initialize", "AlignedReadTrack", function(.Object, coverageOnly=FALSE, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "AlignedReadTrack")
    .Object <- callNextMethod()
    .Object <- setCoverage(.Object)
    if(coverageOnly)
    {
      ## from <- min(unlist(lapply(.Object@coverage, function(y) if(length(y)) min(start(y)))))
      ## to <- max(unlist(lapply(.Object@coverage, function(y) if(length(y)) max(end(y)))))
      from <- min(start(range(.Object)))
      to <- max(end(range(.Object)))
        .Object@range <- GRanges(range=IRanges(start=from, end=to),
                                 strand=names(.Object@coverage), seqnames=.Object@chromosome)
        .Object@coverageOnly <- coverageOnly
    }
    return(.Object)
})

## Constructor. The following arguments are supported:
##    o range: a data.frame or a GRanges object containing the information
##       about the track items. If a data.frame, it needs to be coerceable
##	 to a GRanges object, i.e., it needs at least the mandatory 'start', 'stop' and 
##	 'strand' columns.
##	Instead of using the 'range' parameter, all these values can also be passed as
##	individual vectors, in which case they need to be of similar length.
##    o start, end, width: numeric vectors of the item start and end coordinates, or their widths
##    o strand: the strand information needs to be provided in the form '+' for
##       the Watson strand, '-' for the Crick strand.
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o stacking: character controlling the stacking of overlapping items. One in 'hide', 
##       'dense', 'squish', 'pack' or 'full'.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as further DisplayParameters
## (N)
AlignedReadTrack <- function(range=NULL, start=NULL, end=NULL, width=NULL, chromosome, strand, genome, stacking="squish", 
                             name="AlignedReadTrack", coverageOnly=FALSE, ...)
{
    .missingToNull(c("strand", "chromosome", "genome"))
    ## Build a GRanges object from the inputs
    range <- .buildRange(range=range, start=start, end=end, width=width,
                         args=list(strand=strand, genome=genome, chromosome=chromosome),
                         defaults=list(strand="+", genome=NA, chromosome="chrNA"), chromosome=chromosome,
                         trackType="AlignedReadTrack")
    str <- unique(as.character(GenomicRanges::strand(range)))
    if("*" %in% str)
        stop("Only '+' and '-' strand information is allowed for AlignedReadTrack objects.")
    if(missing(chromosome) || is.null(chromosome))
        chromosome <- if(length(range)>0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    ## And finally the object instantiation
    return(new("AlignedReadTrack", chromosome=chromosome, range=range, name=name, genome=genome(range)[1],
               stacking=stacking, coverageOnly=coverageOnly, ...))
}
##----------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------
## SequenceTrack:
## 
## A generic track to visualize nucleotide sequences. This class is virtual.
## Slots:
##    o chromosome: a character vector giving the active chromosome for which the 
##	 track is defined. Valid chromosome names are: 
##          - a single numeric character
##	    - a string, starting with 'chr', followed by any additional characters
##    o genome: character giving the reference genome for which the track is defined.
## A bunch of DisplayPars are set during object instantiation:
##    o foo: bar
setClass("SequenceTrack",
         representation=representation("VIRTUAL",
                                       chromosome="character",
                                       genome="character"),
         contains="GdObject",
         prototype=prototype(name="Sequence",
                             dp=DisplayPars(size=NULL,
                                            fontcolor=getBioColor("DNA_BASES_N"),
                                            fontsize=10,
                                            fontface=2,
                                            lwd=2,
                                            col="darkgray",
                                            min.width=2,
                                            showTitle=FALSE,
                                            background.title="transparent",
                                            noLetters=FALSE,
                                            complement=FALSE,
                                            add53=FALSE),
                             genome=as.character(NA),
                             chromosome="chrNA"))

## Essentially we just update the display parameters here and set the chromosome and the genome
setMethod("initialize", "SequenceTrack", function(.Object, chromosome, genome, ...) {
    ## the diplay parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "SequenceTrack")
     if(!missing(chromosome) && !is.null(chromosome)){
        .Object@chromosome <- .chrName(chromosome)[1]
    }
    if(missing(genome) || is.null(genome))
        genome <- as.character(NA)
    .Object@genome <- genome
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

## We want the following behaviour in the constructor:
##   a) sequence is missing (NULL) => build SequenceDNAStringSetTrack with chromosome NA and genome as supplied or NA if missing
##   b) sequence is DNAStringSet => build SequenceDNAStringSetTrack where chromosome is names(sequence)[1] or the supplied
##      chromosome if available, and genome as supplied or NA if missing
##   c) sequence is BSgenome => build SequenceBSgenomeTrack where chromosome is seqnames(sequence)[1] or the supplied
##      chromosome if available, and genome is the supplied genome or the one extracted from the BSgenome object
SequenceTrack <- function(sequence, chromosome, genome, name="SequenceTrack", importFunction, stream=FALSE, ...){
    .missingToNull(c("chromosome", "genome", "sequence"))
    if(is.null(sequence)){
        return(new("SequenceDNAStringSetTrack", chromosome=chromosome, genome=genome, name=name, ...))
    }
    if(is(sequence, "BSgenome")){
        if(is.null(genome))
            genome <- providerVersion(sequence)
        if(is.null(chromosome))
            chromosome <- seqnames(sequence)[1]
        obj <- new("SequenceBSgenomeTrack", sequence=sequence, chromosome=chromosome, genome=genome, name=name, ...)
    }else if(is(sequence, "DNAStringSet")){
        if(is.null(names(sequence)))
            stop("The sequences in the DNAStringSet must be named")
        if(any(duplicated(names(sequence))))
            stop("The sequence names in the DNAStringSet must be unique")
        if(is.null(chromosome))
            chromosome <- names(sequence)[1]
        obj <- new("SequenceDNAStringSetTrack", sequence=sequence, chromosome=chromosome, genome=genome, name=name, ...)
    } else if(is.character(sequence)){
        sequence <- sequence[1]
        if(!file.exists(sequence))
            stop(sprintf("'%s' is not a valid file.", sequence))
        ext <- .fileExtension(sequence)
        obj <- if(missing(importFunction) && ext %in% c("fa", "fasta")){
            if(!file.exists(paste(sequence, "fai", sep="."))){
                new("SequenceDNAStringSetTrack", sequence=readDNAStringSet(sequence), chromosome=chromosome,
                    genome=genome, name=name, ...)
            }else{
                new("ReferenceSequenceTrack", chromosome=chromosome, genome=genome, name=name,
                    stream=.import.fasta, reference=path.expand(sequence), ...)
            }
        }else if(missing(importFunction) && ext == "2bit"){
            new("ReferenceSequenceTrack", chromosome=chromosome, genome=genome, name=name,
                stream=.import.2bit, reference=path.expand(sequence), ...)
        }else{
            if(missing(importFunction)){
                stop(sprintf("No predefined import function exists for files with extension '%s'. Please manually provide an import function.",
                             ext))
            }else{
                if(!stream){
                    seq <- importFunction(file=sequence)
                    if(!is(seq, "DNAStringSet"))
                        stop("The import function did not provide a valid DNAStringSet object. Unable to build track from file '",
                             sequence, "'")
                    new("SequenceDNAStringSetTrack", sequence=importFunction(file=sequence), chromosome=chromosome,
                        genome=genome, name=name, ...)
                }else{
                    new("ReferenceSequenceTrack", chromosome=chromosome, genome=genome, name=name,
                        stream=importFunction, reference=path.expand(sequence), ...)
                }
            }
        }
    }else{
        stop("Argument sequence must be of class 'BSgenome', 'DNAStringSet' or 'character'")
    }
    return(obj)  
}
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## SequenceDNAStringSetTrack:
## 
## A track to visualize nucleotide sequences that are stored in a DNSStringSet
## Slots:
##    o sequence: a DNAStringSet object that contains all the sequence data
##----------------------------------------------------------------------------------------------------------------------
setClass("SequenceDNAStringSetTrack",
         representation=representation(sequence="DNAStringSet"),
         contains="SequenceTrack",
         prototype=prototype(sequence=DNAStringSet()))

setMethod("initialize", "SequenceDNAStringSetTrack", function(.Object, sequence, ...) {
    if(missing(sequence) || is.null(sequence))
        sequence <- DNAStringSet()
    .Object@sequence <- sequence
    .Object <- callNextMethod(.Object, ...)
     return(.Object)
})
##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## SequenceBSgenomeTrack:
## 
## A track to visualize nucleotide sequences that are stored in a BSgenome package
## Slots:
##    o sequence: a DNAStringSet object that contains all the sequence data
##    o pointerCache: an environemnt to hold pointers to the BSgenome sequences to prevent garbage collection. This
##       will only be filled once the individual sequences have been accessed for the first time
##----------------------------------------------------------------------------------------------------------------------
setClass("SequenceBSgenomeTrack",
         representation=representation(sequence="BSgenomeOrNULL", pointerCache="environment"),
         contains="SequenceTrack",
         prototype=prototype(sequence=NULL))

setMethod("initialize", "SequenceBSgenomeTrack", function(.Object, sequence=NULL, ...) {
    .Object@sequence <- sequence
    .Object@pointerCache <- new.env()
    .Object <- callNextMethod(.Object, ...)
     return(.Object)
})

##----------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------
## ReferenceSequenceTrack:
## 
## The file-based version of the ReferenceTrack class. This will mainly provide a means to dispatch to
## a special 'subseq' method which should stream the necessary data from disk. 
setClass("ReferenceSequenceTrack", contains=c("SequenceDNAStringSetTrack", "ReferenceTrack"))

## This just needs to set the appropriate slots that are being inherited from ReferenceTrack because the
## multiple inheritence has some strange features with regards to method selection
setMethod("initialize", "ReferenceSequenceTrack", function(.Object, stream, reference, ...) {
    .Object <- selectMethod("initialize", "ReferenceTrack")(.Object=.Object, reference=reference, stream=stream)
    .Object <- callNextMethod()
    return(.Object)
})
