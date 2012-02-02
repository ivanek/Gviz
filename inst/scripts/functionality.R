library(Gviz)

## -------------------------------------------------
## GenomeAxisTrack class
## -------------------------------------------------
ax <- GenomeAxisTrack()
dev.new(height=1.5)
par(ask=TRUE)

## basic plotting
plotTracks(ax, from=1, to=100)

## automatic label setting
plotTracks(ax, from=1, to=10000)
plotTracks(ax, from=1, to=10000000)
plotTracks(ax, from=1, to=10000000000)

## manual label setting
plotTracks(ax, from=1, to=10000000000, exponent=8)
plotTracks(ax, from=-3000000, to=10000000, exponent=4)

## direction indicators
plotTracks(ax, from=1, to=10000, add53=TRUE)
plotTracks(ax, from=1, to=10000, add53=TRUE, add35=TRUE)

## text size
plotTracks(ax, from=1, to=10000, add53=TRUE, add35=TRUE, cex=1.5)
plotTracks(ax, from=1, to=10000, add53=TRUE, add35=TRUE, cex=0.6)

## second level labels
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE)

## axis spacing
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, distFromAxis=2)

## ranges
library(IRanges)
ax <- GenomeAxisTrack(range=IRanges(start=c(1,670), end=c(321, 912)), id=c("foo", "bar"))
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE)

## range annotion
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE)
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE, cex.id=1.5)

## backgrounds
plotTracks(ax, from=1, to=1000, background.panel="#DAE3E6")
plotTracks(ax, from=1, to=1000, background.title="#DAE3E6")

## frame
plotTracks(ax, from=1, to=1000, frame=TRUE)
plotTracks(ax, from=1, to=1000, background.title="#DAE3E6", col.frame="darkblue", background.panel="#FFDDB5", frame=TRUE)

## colors
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE, col="darkgreen")
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE, fontcolor="darkred")
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE, col.id="darkblue")
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE, fill.range="salmon2")

## font settings
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, showId=TRUE, fontfamily="serif", fontface=3)

## line width
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, lwd=3)

## tick mark orientation
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, labelPos="revAlternating")
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, labelPos="above")
plotTracks(ax, from=1, to=1000, add53=TRUE, add35=TRUE, littleTicks=TRUE, labelPos="below")


## -------------------------------------------------
## IdeogramTrack class
## -------------------------------------------------

## basic plotting
plotTracks(id, from=10000000, to=20000000)

## colors
plotTracks(id, from=10000000, to=20000000, col="darkblue", fill="lightblue")
plotTracks(id, from=10000000, to=20000000, background.panel="lightgreen")
plotTracks(id, from=10000000, to=20000000, fontcolor="red", background.title="lightgreen")

## line type
plotTracks(id, from=10000000, to=20000000, lty="dotted", lwd=2)

## font settings
plotTracks(id, from=10000000, to=20000000, fontsize=8)
plotTracks(id, from=10000000, to=20000000, cex=2)
plotTracks(id, from=10000000, to=20000000, fontface=3, fontcolor="red", fontfamily="serif")

## shape
plotTracks(id, from=10000000, to=20000000, bevel=0)


## ------------------
## Object creation:
## ------------------


## Annotation track:
## As individual arguments
annTrack <- AnnotationTrack(chromosome=7, feature="test", group=c(1,2,1), 
                            start=c(2000000, 2070000, 2100000), ID=paste("annTrack item", 1:3),  
                            end=c(2050000, 2130000, 2150000), genome="hg19", name="annTrack", 
                            stacking="squish", strand=c("-", "+", "-"), size=1.5)
## Or as data.frame
annTrack2 <- AnnotationTrack(data.frame(start=c(2000000, 2070000, 2100000), 
                                        ID=paste("annTrack item", 1:3), end=c(2050000, 2130000, 2170000), 
                                        feature="test", group=c("Group1","Group2","Group1"), strand=c("-", "+", "-")),
                             genome="hg19", name="annTrack2", 
                             stacking="squish",chromosome=7, size=1.5, showId=TRUE)
names(annTrack) <- "testAnnTrack"
plotTracks(list(annTrack, annTrack2), extend.left=20000)

## The track name


from <- 2000000
to <- 2100000

it <- IdeogramTrack(genome="hg19", name="chromosome 4", chromosome=4)
plotTracks(list(it, annTrack2))


## GeneRegion track
geneTrack <- BiomartGeneRegionTrack(chromosome=4, start=2000000, end=2100000, genome="hg19",
                                    stacking="squish", strand="+-", name="geneTrack")
## Todo: filter testing
plotTracks(geneTrack)

gt1 <- BiomartGeneRegionTrack(chromosome=4, start=2000000, end=2100000, genome="hg19",
                              stacking="squish", strand="+", name="+ Strand")
gt2 <- BiomartGeneRegionTrack(chromosome=4, start=2000000, end=2100000, genome="hg19",
                              stacking="squish", strand="-", name="- Strand")
gt3 <- BiomartGeneRegionTrack(chromosome=4, start=2000000, end=2100000, genome="hg19",
                              stacking="squish", strand="+-", name="RefSeq", filter=Gviz:::.refseqFilter)
plotTracks(list(geneTrack, gt1, gt2, gt3))

displayPars(gt1) <- list(showId=TRUE, fontsize=8)
displayPars(gt2) <- list(showId=TRUE, fontsize=8)
displayPars(gt3) <- list(showId=TRUE, fontsize=8)
plotTracks(list(gt1, gt2, gt3), extend.left=10000)




ga <- GenomeAxisTrack()
ga <- GenomeAxisTrack(range=GenomicRanges::GRanges(seqnames=letters[1:3], range=IRanges::IRanges(start=c(2045000, 2061000, 2151000), end=c(2050000, 2066000, 2221000))),
                      add53=TRUE, add35=TRUE)
plotTracks(list(gt1, ga, gt2, gt3), extend.left=10000)
displayPars(ga) <- list(littleTicks=TRUE)
plotTracks(list(annTrack2, gt1, ga, gt2, gt3), extend.left=10000)



bases <- seq(min(IRanges::start(range(gt2))), max(IRanges::end(range(gt2))), len=50)
start <- bases-(min(diff(bases))/5)*runif(50,1,3)
end <- bases+(min(diff(bases))/5)*runif(50,1,3)
data <- matrix(rnorm(length(bases)*6)+rep(runif(length(bases), -4,4), each=6), nrow=6)

dt <- DataTrack(start=start, end=end, data=data)

plotTracks(list(gt2, ga, dt), extend.left=10000)


dtp <- DataTrack(start=start, end=end, data=data, type="p", name="points")
dtl <- DataTrack(start=start, end=end, data=data, type="l", name="lines")
dtb <- DataTrack(start=start, end=end, data=data, type="b", name="both")
dta <- DataTrack(start=start, end=end, data=data, type="a", name="average")
dts <- DataTrack(start=start, end=end, data=data, type="s", name="steps")
dtg <- DataTrack(start=start, end=end, data=data, type="g", name="grid")
dtr <- DataTrack(start=start, end=end, data=data, type="r", name="regression")
dtsm <- DataTrack(start=start, end=end, data=data, type="smooth", name="loess")
dth <- DataTrack(start=start, end=end, data=data, type="h", name="hist")
dtm <- DataTrack(start=start, end=end, data=data, type="mountain", name="mountain")
dtbp <- DataTrack(start=start, end=end, data=data, type="boxplot", name="boxplot")
dthi <- DataTrack(start=start, end=end, data=data, type="histogram", name="histogram")

plotTracks(list(dtp, dtl, dtb, dta, dts, dtg, dtr, dtsm, dth))

plotTracks(list(dtm, dtbp, dthi))

trunc <- function(vals)
{
    trunc <- c(30,70)
    vals[vals<trunc[1]] <- trunc[1]
    vals[vals>trunc[2]] <- trunc[2]
    return(vals)
}

ut1 <- UcscTrack(genome="hg19", chromosome=4, track="GC", from=1985000, to=2150000,
                 trackType="DataTrack", start="start", end="end", data="score", type="gradient",
                 transformation=trunc, window="auto")

ut2 <- UcscTrack(genome="hg19", chromosome=4, track="Ensembl", from=1900000, to=2300000,
                 table="ensGene", trackType="Gene", rstarts="exonStarts", rends="exonEnds", gene="name2",
                 symbol="name2", transcript="name", strand="strand", ID="name",
                 feature="default", exonList=TRUE, showId=TRUE, fontsize=8)


ut4 <- UcscTrack(genome="canFam2", chromosome=10, track="Other RefSeq", from=13910391-10000, to=13956580+10000,
                 trackType="Gene", rstarts="exonStarts", rends="exonEnds", gene="name2",
                 symbol="name2", transcript="name", strand="strand", ID="name",
                 feature="default", exonList=TRUE, showId=TRUE, fontsize=8)


ut3 <-  UcscTrack(genome="hg19", chromosome=4, track="Repeat", from=1900000, to=2300000,
                  trackType="Annotation", start="genoStart", end="genoEnd", ID="repName",
                  feature="repClass", strand="strand", shape="box", stacking="dense", size=0.5, col="black")

ut4 <- UcscTrack(genome="hg19", chromosome=4, track="knownGene", from=1900000, to=2300000,
                 trackType="Gene", rstarts="exonStarts", rends="exonEnds", gene="name",
                 symbol="name", transcript="name", strand="strand", ID="name",
                 feature="default", exonList=TRUE, showId=TRUE, fontsize=8)

plotTracks(list(ut1, ga, ut2), from=1950000, to=2050000)

plotTracks(list(it, geneTrack, ga, ut2), fontsize=10, showId=TRUE, from=1985000, to=2150000)

plotTracks(list(ut1, ga, ut2, ut3), type="h", from=1985000, to=2150000)


                  


ut4 <-  UcscTrack(genome="hg19", chromosome=4, track="snp131", from=1900000, to=2200000,
                  trackType="Annotation", start="chromStart", end="chromEnd", ID="name",
                  feature="func", strand="strand", shape="box", stacking="squish", size=0.6, fill="black",
                  "coding-synon"="green", nonsense="red", missense="red", frameshift="red",
                  "untranslated-3"="blue", "untranslated-5"="blue", unknown="gray")

plotTracks(list(it, ut1, ga, ut2, ut4, ut3), type="h", from=1985000, to=2150000)
