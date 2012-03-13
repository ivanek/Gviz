## ============================================================================================================================
## All generic functions defined in the package
## ============================================================================================================================

## Display parameter accessors
setGeneric("setPar", function(x, value, ...) standardGeneric("setPar"))
setGeneric("displayPars<-", function(x, value) standardGeneric("displayPars<-"))
setGeneric("getPar", def = function(x, name, ...) standardGeneric("getPar"))
setGeneric("displayPars", function(x, name, ...) standardGeneric("displayPars"))

## Annotation accessors
setGeneric("gene", function(GdObject, ...) standardGeneric("gene"))
setGeneric("gene<-", function(GdObject, value) standardGeneric("gene<-"))
setGeneric("symbol", function(GdObject, ...) standardGeneric("symbol"))
setGeneric("symbol<-", function(GdObject, value) standardGeneric("symbol<-"))
setGeneric("transcript", function(GdObject, ...) standardGeneric("transcript"))
setGeneric("transcript<-", function(GdObject, value) standardGeneric("transcript<-"))
setGeneric("exon", function(GdObject, ...) standardGeneric("exon"))
setGeneric("exon<-", function(GdObject, value) standardGeneric("exon<-"))
setGeneric("feature", function(GdObject, ...) standardGeneric("feature"))
setGeneric("feature<-", function(GdObject, value) standardGeneric("feature<-"))
setGeneric("group", function(GdObject, ...) standardGeneric("group"))
setGeneric("group<-", function(GdObject, value) standardGeneric("group<-"))
setGeneric("identifier", function(GdObject, ...) standardGeneric("identifier"))
setGeneric("identifier<-", function(GdObject, value) standardGeneric("identifier<-"))


## General accessors
setGeneric("chromosome", function(GdObject, ...) standardGeneric("chromosome"))
setGeneric("chromosome<-", function(GdObject, value) standardGeneric("chromosome<-"))
setGeneric("[")
setGeneric("position", function(GdObject, ...) standardGeneric("position"))
setGeneric("imageMap", function(GdObject, ...) standardGeneric("imageMap"))
setGeneric("imageMap<-",  function(GdObject, value) standardGeneric("imageMap<-"))
##setGeneric("subset",  function(x, ...) standardGeneric("subset"))
setGeneric("coords",  function(ImageMap, ...) standardGeneric("coords"))
setGeneric("tags",  function(ImageMap, ...) standardGeneric("tags"))

## Prepare tracks for plotting
setGeneric("consolidateTrack", function(GdObject, ...) standardGeneric("consolidateTrack"))
setGeneric("collapseTrack", function(GdObject, ...) standardGeneric("collapseTrack"))
setGeneric("stacking", function(GdObject, ...) standardGeneric("stacking"))
setGeneric("stacking<-", function(GdObject, value) standardGeneric("stacking<-"))
setGeneric("stacks", function(GdObject, ...) standardGeneric("stacks"))
setGeneric("setStacks", function(GdObject, ...) standardGeneric("setStacks"))
setGeneric("setCoverage", function(GdObject, ...) standardGeneric("setCoverage"))

## Plotting methods
setGeneric("drawAxis", function(GdObject, ...) standardGeneric("drawAxis"))
setGeneric("drawGrid", function(GdObject, ...) standardGeneric("drawGrid"))
setGeneric("drawGD", function(GdObject, ...) standardGeneric("drawGD"))

## We may need those for dispatch in the name space
##if(!isGeneric("lapply"))
##    setGeneric("lapply")
##if(!isGeneric("sapply"))
##    setGeneric("sapply")
##if(!isGeneric("head"))
##    setGeneric("head")
##if(!isGeneric("split"))
##    setGeneric("split")

## Internal methods 
setGeneric(".buildRange",  function(range, start, end, width, ...) standardGeneric(".buildRange"))
