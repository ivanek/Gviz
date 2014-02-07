;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.


.gassign <- function(what) assignInNamespace(ns="Gviz", x=what, value=get(what, envir=globalenv()))

assignInNamespace(ns="Gviz", x=".computeGroupRange", value=.computeGroupRange)
.dpOrDefault <- Gviz:::.dpOrDefault
.getStringDims <- Gviz:::.getStringDims
.fontGp <- Gviz:::.fontGp

assignInNamespace(ns="Gviz", x=".barsAndLabels", value=.barsAndLabels)
.pxResolution <- Gviz:::.pxResolution
.getBiotypeColor = Gviz:::.getBiotypeColor

assignInNamespace(value=.dpOrDefault, x=".dpOrDefault", ns="Gviz")
assignInNamespace(value=.dpOrDefaultFont, x=".dpOrDefaultFont", ns="Gviz")
assignInNamespace(value=.fontGp, x=".fontGp", ns="Gviz")

library(Biobase)
assignInNamespace(ns="Gviz", x=".defaultRange", value=.defaultRange)
.dpOrDefault <- Gviz:::.dpOrDefault

assignInNamespace(ns="Gviz", x=".arrowBar", value=.arrowBar)
.pxResolution <- Gviz:::.pxResolution

assignInNamespace(ns="Gviz", x=".boxes", value=.boxes)
.dpOrDefault <- Gviz:::.dpOrDefault
.getBiotypeColor = Gviz:::.getBiotypeColor
.getAnn = Gviz:::.getAnn
.getImageMap <- Gviz:::.getImageMap

assignInNamespace(ns="Gviz", x=".filledBoxes", value=.filledBoxes)
.handleComposite <- Gviz:::.handleComposite

assignInNamespace(ns="Gviz", x=".filledArrow", value=.filledArrow)
.handleComposite <- Gviz:::.handleComposite

assignInNamespace(ns="Gviz", x=".updatePars", value=.updatePars)

assignInNamespace(ns="Gviz", x=".handleComposite", value=.handleComposite)

assignInNamespace(ns="Gviz", x=".setupTextSize", value=.setupTextSize)
.dpOrDefault <- Gviz:::.dpOrDefault
vpLocation <- Gviz:::vpLocation
.verticalSpace <- Gviz:::.verticalSpace
.needsAxis <- Gviz:::.needsAxis
.PLOT_TYPES <- Gviz:::.PLOT_TYPES

assignInNamespace(ns="Gviz", x="plotTracks", value=plotTracks)
.whichStrand <- Gviz:::.whichStrand
.recChromosome <- Gviz:::.recChromosome
.needsTitle <- Gviz:::.needsTitle
.fontGp <- Gviz:::.fontGp
drawAxis <- Gviz:::drawAxis
.getImageMap <- Gviz:::.getImageMap
drawGrid <- Gviz:::drawGrid
ImageMap <- Gviz:::ImageMap

## This is all needed for sourcing plotTracks
.defaultRange <- Gviz:::.defaultRange
setStacks <- Gviz:::setStacks
vpLocation <- Gviz:::vpLocation
.setupTextSize <- Gviz:::.setupTextSize
.needsAxis <- Gviz:::.needsAxis
drawAxis <- Gviz:::drawAxis
.getImageMap <- Gviz:::.getImageMap
drawGrid <- Gviz:::drawGrid
ImageMap <- Gviz:::ImageMap
.needsTitle <- Gviz:::.needsTitle
.whichStrand <- Gviz:::.whichStrand
.recChromosome <- Gviz:::.recChromosome


assignInNamespace(ns="Gviz", x=".needsTitle", value=.needsTitle)
assignInNamespace(ns="Gviz", x=".needsAxis", value=.needsAxis)

library(Gviz)
data(cyp2b10)
dt<-GeneRegionTrack(cyp2b10)
plotTracks(dt)

plotTracks(dt, showId=T, just.group="left", cex.group=0.8, reverseStrand=F, title.width=1)



plotTracks(c(dt), just.group="left", cex.group=1, reverseStrand=F, title.width=1, geneSymbols="symbol", fontface.group="bold.italic", fontcolor.group="red", showExonId=F, fontcolor.item=1, from=25897620-500, to=25897855, fontfamily="serif", fontfamily.item="sans", showId=T, just.group="left")


plotTracks(c(dt), showId=T, alpha=1, lwd=1, col=1, debug="draw", collapseTranscripts=F)


plotTracks(dt, col=1, shape=c("smallArrow", "arrow"))


plotTracks(list(GenomeAxisTrack(), dt), transcriptAnnotation="symbol", fill="#FFD58A", just.group="below", extend.left=-0.9)

## This is broken:
plotTracks(list(GenomeAxisTrack(), dt), fill="#FFD58A", showId=T, just.group="below", extend.left=-0.99, extend.right=-0.99)


st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
str <- c("-", "+", "-", "-")
gr <- c("Group1","Group2","Group1", "Group3")

annTrack <- AnnotationTrack(start=st, end=ed, strand=str, chromosome=7, genome="hg19", feature="test",
                            group=gr, id=paste("annTrack item", 1:4), name="generic annotation", stacking="squish")
plotTracks(annTrack)




