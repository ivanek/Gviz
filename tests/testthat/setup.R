## libs ---------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)

## general ------------------------------------------------------------------
set.seed(789)
ir <- IRanges(1L, 5L)
gr <- GRanges("chr1", ir, score = 1)
ir2 <- IRanges(2L, 6L)
gr2 <- GRanges("chr1", ir2, score = 2)
dna.sq <- DNAStringSet(c(chr1 = paste(sample(DNA_BASES, 100, replace = TRUE), collapse = "")))
rna.sq <- RNAStringSet(c(chr1 = paste(sample(RNA_BASES, 100, replace = TRUE), collapse = "")))
cyto.bands <- data.frame(
    chrom = rep(c("chrI", "chrII"), each = 4),
    chromStart = rep(c(1L, 148071L, 151524L, 154977L), 2),
    chromEnd = rep(c(148071L, 151524L, 154977L, 230218L), 2),
    name = rep(c(NA, "CEN1", "CEN1", NA), 2),
    gieStain = rep(c("gneg", "acen", "acen", "gneg"), 2),
    stringsAsFactors = FALSE
)
## internet access ----------------------------------------------------------
# hasUcscConnection <- !is(try(rtracklayer::browserSession(), silent = TRUE), "try-error")
hasUcscConnection <- FALSE
check_ucsc <- function() {
    if (!hasUcscConnection) {
        skip("UCSC not available")
    }
}

#
# oto <- options(timeout = 5)
# hasBiomartConnection <- (!is(try(download.file("http://www.biomart.org", tempfile(), quiet = TRUE)), "try-error") &&
#     !is(try(biomaRt::listMarts(), silent = TRUE), "try-error"))
#options(timeout = oto)
hasBiomartConnection <- FALSE
check_biomart <- function() {
    if (!hasBiomartConnection) {
        skip("Biomart not available")
    }
}

## Uncommenting this helps when the UCSC server has a hickup but still lets you connect:
## hasUcscConnection <- !is(try(rtracklayer::browserSession(), silent=TRUE), "try-error") &&
##   !is(try(IdeogramTrack(genome="hg19", chromosome=7), silent=TRUE), "try-error")

## tmp files ----------------------------------------------------------------
## BAM file
sam <- c(
    "@HD	VN:1.0	GO:none	SO:coordinate",
    "@SQ	SN:chr1	LN:267910886",
    "@PG	ID:SpliceMap	VN:3.3.5.2 (55)",
    "@CO	file merged using Picard tools",
    "NRCHBS-WDL30299:125:D1415ACXX:8:2215:20868:43279	65	chr1	189891483	255	15S51M6207N10M	=	189897797	47	CGGGACCGTGGTGAGTCAGCGAGTAGGAACTACTCAGGAACTACTCACTTCACTGACAGCCGTAAGTCACTCTGAC	;99;;;5((2@:-)()2:?<88;196?<??1))<>=9?9>>?<>??8>1>>=?;?>?;3=?<9===9=>==>>77=	XS:A:+	NM:i:0	XC:i:15",
    "NRCHBS-WDL30299:126:D14UTACXX:8:2113:17433:81932	177	chr1	189892202	255	76M	=	189892390	113	GCCAAATCAATTCTTTCTTCTCCTAAGCTGCTTTGGTTAGGGTCGTTTATCACAGCAACAGAAAAGTAAGTAAAAC	CC@FFFFFHHHHHJJJIJJIJJJJIJJJJIJJJJJIJJJIIJDGHHIFJJJJIJJJIJJJJJFIJJHHHHHHHHFC	NM:i:0	XC:i:0",
    "NRCHBS-WDL30299:126:D14UTACXX:8:2113:17433:81932	113	chr1	189892390	255	76M	=	189892202	113	GGGCTTTTGCAGGGTTGGGCTGTCCTCCCCTGCCTCTTCTCCAGCCCTTTAATCTTTCAGCACCCAGATCAAGCCG	CCCFFFFFHHHHHJCGHGIIJJIIJJIJJJIJJIJJGHHJIJIIJJJJJIIIGEHHFFFFDDDDDDDDDDDDCBD9	MD:Z:14C61	NM:i:1	XC:i:0",
    "NRCHBS-WDL30299:125:D1415ACXX:8:1313:16113:75579	65	chr1	189893347	255	76M	=	189893470	48	CTGCAGCGGCTCATACAAGACCTTTGCTTAATAGGGTTTTCTTTACCCTAGAAAAATTCCCCTGTGATTTAAGGGT	CCCFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJJJFHIJJJJJJJJJJJJJJJJIJJJJJJGHHHHFFFFFFDEA	NM:i:0	XC:i:0",
    "NRCHBS-WDL30299:126:D14UTACXX:8:2215:4526:88225	65	chr1	189893352	255	76M	=	189893457	30	GCGGCTCATACAAGACCTTTGCTTAATAGGGTTTTCTTTACCCTAGAAAAATTCCCCTGTGATTTAAGGGTACAGA	@@@DDADDHAF3CGBFH@HIIIEGIGEEHII?CBGDAFHIGFFGIIEGCH@EHGBFHEGCGGCCHICEGG6ACEHE	NM:i:0	XC:i:0"
)
samfile <- tempfile(fileext = ".sam")
cat(paste(sam, collapse = "\n"), file = samfile)
bamfile <- Rsamtools::asBam(samfile, indexDestination = FALSE)
Rsamtools::indexBam(bamfile)

## FASTA
fastafile <- system.file("extdata/test.fa", package = "Gviz")

bamgr <- GRanges("chr1", IRanges(
    c(189892390L, 189892202L, 189893347L, 189891483L, 189893352L),
    c(189892465L, 189892277L, 189893422L, 189891558L, 189893427L)
),
strand = rep(c("-", "+"), c(2, 3))
)
mcols(bamgr) <- DataFrame(
    id = c(
        "NRCHBS-WDL30299:126:D14UTACXX:8:2113:17433:81932",
        "NRCHBS-WDL30299:126:D14UTACXX:8:2113:17433:81932",
        "NRCHBS-WDL30299:125:D1415ACXX:8:1313:16113:75579",
        "NRCHBS-WDL30299:125:D1415ACXX:8:2215:20868:43279",
        "NRCHBS-WDL30299:126:D14UTACXX:8:2215:4526:88225"
    ),
    cigar = c("76M", "76M", "76M", "15S51M6207N10M", "76M"),
    mapq = rep(255L, 5),
    flag = c(113L, 177L, 65L, 65L, 65L),
    md = c("14C61", NA, NA, NA, NA),
    seq = DNAStringSet(c(
        `1` = paste(rep("+", 2600), collapse = ""),
        `2` = paste(rep("+", 2600), collapse = ""),
        `3` = paste(rep("+", 2600), collapse = ""),
        `4` = paste(rep("+", 2600), collapse = ""),
        `5` = paste(rep("+", 2600), collapse = "")
    )),
    isize = c(113L, 113L, 48L, 47L, 30L),
    groupid = c(1L, 1L, 2L, 3L, 4L),
    status = factor(rep(c("mated", "unmated"), c(2, 3)),
        levels = c("mated", "ambiguous", "unmated")
    )
)



selFun <- function(identifier, start, end, track, GdObject, ...) {
    gcount <- table(group(GdObject))
    ## This computes the width of 2 pixels in genomic coordinates
    pxRange <- Gviz:::.pxResolution(min.width = 20, coord = "x")
    return((end - start) < pxRange && gcount[identifier] == 1)
}
detFun <- function(identifier, GdObject.original, ...) {
    plotTracks(list(
        GenomeAxisTrack(scale = 0.3, size = 0.2, cex = 0.7),
        GdObject.original[group(GdObject.original) == identifier]
    ),
    add = TRUE, showTitle = FALSE
    )
}

data(geneDetails)
data(geneModels)

## classes ------------------------------------------------------------------

## IdeogramTrack
ideoTrack <- IdeogramTrack(chromosome = "chrI", genome = "sacCer3", bands = cyto.bands)

## GenomeAxisTrack
axisTrack <- GenomeAxisTrack(gr)

## DataTrack
dataTrack <- DataTrack(gr)

## AnnotationTrack
annoTrack <- AnnotationTrack(gr)

## GeneRegionTrack
geneTrack <- GeneRegionTrack(geneModels, genome = "hg19", chromosome = "chr7", name = "foo")
## BiomartGeneRegionTrack

## DetailsAnnotationTrack
detTrack <- DetailsAnnotationTrack(geneDetails,
    fun = detFun, selectFun = selFun,
    groupDetails = TRUE, details.size = 0.5,
    detailsConnector.cex = 0.5, detailsConnector.lty = "dotted",
    shape = c("smallArrow", "arrow"), groupAnnotation = "group"
)

## SequenceTrack
seqTrack.dna <- SequenceTrack(dna.sq)
seqTrack.rna <- RNASequenceTrack(rna.sq)
seqTrack.bs <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg19)

## AlignmentsTrack
alnTrack <- AlignmentsTrack(bamfile)

## HighlightTrack
highTrack <- HighlightTrack(dataTrack, ranges = gr, chromosome = "chr1")

## OverlayTrack
overTrack <- OverlayTrack(c(dataTrack, dataTrack))
