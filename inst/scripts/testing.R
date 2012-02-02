##load("~muellar2/projects/Epigenetics/Methylation/HD2.1-MM_13-Week_Liver/exMA.Rdata")
##load("~muellar2/projects/Epigenetics/Methylation/HD2.1-MM_13-Week_Liver/probeAnno.Rdata")





panel.genomeGraphs <- function (x, y, fullData, subscripts, coord.factor, tracks, sizes, ...) {
	xy <- split(y,x)
	if(length(unique(listLen(xy))) != 1)
		stop("Supplied data is not a rectangular array.")
	data <- t(sapply(xy, function(x) x))
	dataTrack <- makeMatrixRangeTrack(start=sort(unique(fullData[subscripts, "ProbeStart"])), 
									  end=sort(unique(fullData[subscripts, "ProbeEnd"])), data=data)
	args <- list(...)
	setPar(dataTrack, args)
	cp <- current.panel.limits()
	setPar(dataTrack, "ylim", cp$ylim)
	##browser()
	tracks <- c(list(dataTrack), tracks)
	lt <- length(tracks)
	if(length(sizes) != lt)
		sizes <- rep(sizes, lt)[seq_len(lt)]
	plotTracks(tracks, panel.only=TRUE, from=cp$xlim[1]*coord.factor,  to=cp$xlim[2]*coord.factor, 
			   coord.factor=coord.factor, sizes=sizes)	
}





### objects exProbeAnno and exMA are required

### plots a trellis graphics for the 13 week study for a given gene
### using all time points and control and treated
###
### chr: the chromosome name (e.g. "7", "X")
### start: where to start (e.g. nucleotide 26682683)
### exProbeAnno: a probeAnno object from the Ringo package
### maList: the MAList object with the normalized data
### upstream: how much to extend upstream from the start (in nucleotides), must be a negative integer
### downstream: how much to extend downstream of the start (in nucleotides), must be a positive integer
plotAll = function(exProbeAnno, maList, gene=NULL, start=NULL, chr=NULL, upstream=10000, downstream=500, title=NULL,
		ylim, type=c("mountain", "g"), smooth.span=1/12, coord.factor=1000, tracks=list(), track.sizes=1, ...) {
	
	require(IRanges)
	require(lattice)
	require(org.Mm.eg.db) # for chromosomal location
	if ( is.null(start) ) {
		if ( is.null(gene) ) {
			stop("must either specify gene symbol (gene argument) or start (start argument)\n")
		}
		entrezId = suppressWarnings(as.numeric(gene))
		if ( is.na(entrezId) ) {
			entrezId = sym2entrezId(gene)
			if ( is.null(entrezId) ) {
				stop(paste("no entrezId found for gene", gene, "\n"))
			} else if ( length(entrezId) > 1 ) {
				stop("ambiguous mappings for gene ", gene, ": ", paste(entrezId, collapse=", "))
			}
			loc = mget(entrezId, org.Mm.egCHRLOC)
			if ( max(loc[[1]]) < 0 ) {
				loc = mget(entrezId, org.Mm.egCHRLOCEND)
			}
			chr = unique(names(loc[[entrezId]]))
			if ( length(chr) != 1 ) {
				stop(chr, " - cannot handle with chromosome!")
			}
			start = loc[[entrezId]]
			if ( length(start) > 1 ) { # several transcripts
				cat(paste("found several transcripts: ", paste(start, collapse=", "), "\n"))
				if ( min(start) < 0 && max(start) > 0 ) {
					stop(start, " - cannot handle two transcripts on opposite strands")
				}               
				if ( start[1] < 0 ) { # negative strand
					start = max(start)
				} else {
					start = min(start)
				}
			}           
			if ( is.na(start) ) {
				stop("no chromosomal location found for entrezId ", entrezId)
			}                     
			cat(paste("found", gene, "with entrezId", entrezId, "on chr", chr, "at position", start, ": up =", upstream,
							"to", downstream, "\n"))
		}
	} else {
		if ( is.null(chr) ) {
			stop("chromosome cannot be null")
		}
		gene = paste("position", as.character(start))
	}
	
	if ( start < 0 ) {
		start = abs(start)
		tmp = downstream
		downstream = upstream 
		upstream = tmp
	}
	
	cat("finding overlapping probes: ")
	gene.ir = IRanges(start=start-upstream, end=start+downstream)
	chr.ir = IRanges(start=exProbeAnno[paste(chr, "start", sep=".")], end=exProbeAnno[paste(chr, "end", sep=".")])
	gene.probes = findOverlaps(gene.ir, chr.ir)
	hits = subjectHits(gene.probes)
	hit.names = exProbeAnno[paste(chr, "index", sep=".")][hits] # names/indexes are not necessarily unique (why?)!
	probe.idx = (maList$genes$PROBE_ID %in% hit.names)
	m = maList[probe.idx,]
	probe.pos = start(chr.ir[hits]) + (end(chr.ir[hits]) - start(chr.ir[hits]))/2
	probe.start = start(chr.ir[hits])
	probe.end = end(chr.ir[hits])
	
	### create a data frame for lattice
	cat("creating data frame ...\n")
	d <- data.frame(Mvalue=as.numeric(m$M), ProbeStart=probe.start, ProbeEnd=probe.end,
					 Probe=m$genes$PROBE_ID, Slide=factor(rep(m$targets$SlideNumber, each=nrow(m))),
					 Treatment=factor(rep(m$targets[,"Treatment"], each=nrow(m))),
					 Day=ordered(rep(m$targets[,"Day"], each=nrow(m))), stringsAsFactors=FALSE)
	d$ProbeCenter <- rowMeans(cbind(d$ProbeStart, d$ProbeEnd))
	d$Group <- factor(paste(d$Treatment, d$Day, d$Probe, sep=":"))
	
	ir = range(chr.ir[hits])
	from = min(ir)
	to = max(ir)
	len = width(ir)
	if ( is.null(title) ) {
		title = paste("Averaged and smoothed methylation signals for ", gene, " within -", upstream,
				"/+", downstream, " bp", sep="")
	}
	
	cat(start, start+downstream/2, "\n")
	xyplot(Mvalue ~ (ProbeCenter/coord.factor)|Day+Treatment, data=d, type=type,
		   span=smooth.span, xlim=c((start-upstream)/coord.factor, max(c(extendrange(end(chr.ir[hits])), (start+downstream)))/coord.factor), 
		   ylim=ylim, main=title, xlab="probe position", ylab="M-value (log2(IP) - log2(total))", coord.factor=coord.factor,
		   scales=list(x=list(rot=90)), panel=panel.genomeGraphs, fullData=d, tracks=tracks, sizes=track.sizes, ...)
}

## library(biomaRt)
## bm <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
## gene <- makeGeneRegion(biomart=bm, chromosome="chr7", genome="mm9", start=26670000, end=26690000)
## an <- makeAnnotationTrack(start=c(26678000, 26690000), end=c(26680000, 26900000), ID=c("ann A", "ann B"))
## setPar(an, "plotId", TRUE)
## plotAll(exProbeAnno, exMA, "Cyp2b10", layout=c(7,2), cex=2,pch=".")
## plotAll(exProbeAnno, exMA, "Cyp2b10", type=c("mountain", "boxplot"), tracks=list(an, gene), track.sizes=c(5,1,1), downstream=20000)

### convert a gene symbol into it's entrez Id
sym2entrezId <- function(sym) {
	eg <- NULL
	require(org.Mm.eg.db)
	eg <- org.Mm.egSYMBOL2EG[[sym]]
	if(is.null(eg)) ### nothing found - lookup synonyms
		eg <- org.Mm.egALIAS2EG[[sym]]
	return(eg)
}


library(BSgenome.Mmusculus.UCSC.mm9)
xx <- as(subseq(Mmusculus[[7]], 26670000, 26690000), "DNAString")
start <- 26670000
stop <- 26690000
window <- 100
by <- IRanges(start=seq(start, stop-window-1, by=10), width=window)
cg <- DNAString("CG")
tmp =  aggregate(as(Mmusculus[[7]], "DNAString"), by, FUN=function(x) countPattern(cg, x))
xx <- (end(by)-start(by))/2
xxs <- loess.smooth(xx, tmp, span=1/20, evaluation=1000)
cpg <- makeBaseTrack(base=xxs$x, value=xxs$y)
setPar(cpg, "type", "l")
plotAll(exProbeAnno, exMA, "Cyp2b10", type=c("mountain"), tracks=list(cpg, gene), track.sizes=c(5,1,1), downstream=2000, ylim=c(-0.3, 0.3))



plotAll(exProbeAnno, exMA, "Cyp2b10", type=c("p"), tracks=list(cons, cpg, gene), track.sizes=c(5,1,1), downstream=2000)



