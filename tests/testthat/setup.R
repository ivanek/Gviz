## libs -----------------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)

## general --------------------------------------------------------------------
set.seed(789)
ir <- IRanges(1L,5L)
gr <- GRanges("chr1", ir, score=1)
ir2 <- IRanges(2L,6L)
gr2 <- GRanges("chr1", ir2, score=2)
dna.sq <- DNAStringSet(c(chr1=paste(sample(DNA_BASES, 100, replace=T), collapse="")))
rna.sq <- RNAStringSet(c(chr1=paste(sample(RNA_BASES, 100, replace=T), collapse="")))
cyto.bands <- data.frame(chrom = rep("chrI", 4),
                         chromStart = c(1L, 148071L, 151524L, 154977L),
                         chromEnd = c(148071L, 151524L, 154977L, 230218L),
                         name = c(NA, "CEN1", "CEN1", NA),
                         gieStain = c("gneg", "acen", "acen", "gneg"),
                         stringsAsFactors = FALSE)
## tmp files ------------------------------------------------------------------
## BAM file
sam <- c("@HD	VN:1.0	GO:none	SO:coordinate",
         "@SQ	SN:chr1	LN:267910886",
         "@PG	ID:SpliceMap	VN:3.3.5.2 (55)",
         "@CO	file merged using Picard tools",
         "NRCHBS-WDL30299:125:D1415ACXX:8:2215:20868:43279	65	chr1	189891483	255	15S51M6207N10M	=	189897797	47	CGGGACCGTGGTGAGTCAGCGAGTAGGAACTACTCAGGAACTACTCACTTCACTGACAGCCGTAAGTCACTCTGAC	;99;;;5((2@:-)()2:?<88;196?<??1))<>=9?9>>?<>??8>1>>=?;?>?;3=?<9===9=>==>>77=	XS:A:+	NM:i:0	XC:i:15",
         "NRCHBS-WDL30299:126:D14UTACXX:8:2113:17433:81932	177	chr1	189892202	255	76M	=	189892390	113	GCCAAATCAATTCTTTCTTCTCCTAAGCTGCTTTGGTTAGGGTCGTTTATCACAGCAACAGAAAAGTAAGTAAAAC	CC@FFFFFHHHHHJJJIJJIJJJJIJJJJIJJJJJIJJJIIJDGHHIFJJJJIJJJIJJJJJFIJJHHHHHHHHFC	NM:i:0	XC:i:0",
         "NRCHBS-WDL30299:126:D14UTACXX:8:2113:17433:81932	113	chr1	189892390	255	76M	=	189892202	113	GGGCTTTTGCAGGGTTGGGCTGTCCTCCCCTGCCTCTTCTCCAGCCCTTTAATCTTTCAGCACCCAGATCAAGCCG	CCCFFFFFHHHHHJCGHGIIJJIIJJIJJJIJJIJJGHHJIJIIJJJJJIIIGEHHFFFFDDDDDDDDDDDDCBD9	MD:Z:14C61	NM:i:1	XC:i:0",
         "NRCHBS-WDL30299:125:D1415ACXX:8:1313:16113:75579	65	chr1	189893347	255	76M	=	189893470	48	CTGCAGCGGCTCATACAAGACCTTTGCTTAATAGGGTTTTCTTTACCCTAGAAAAATTCCCCTGTGATTTAAGGGT	CCCFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJJJFHIJJJJJJJJJJJJJJJJIJJJJJJGHHHHFFFFFFDEA	NM:i:0	XC:i:0",
         "NRCHBS-WDL30299:126:D14UTACXX:8:2215:4526:88225	65	chr1	189893352	255	76M	=	189893457	30	GCGGCTCATACAAGACCTTTGCTTAATAGGGTTTTCTTTACCCTAGAAAAATTCCCCTGTGATTTAAGGGTACAGA	@@@DDADDHAF3CGBFH@HIIIEGIGEEHII?CBGDAFHIGFFGIIEGCH@EHGBFHEGCGGCCHICEGG6ACEHE	NM:i:0	XC:i:0")
samfile <- tempfile(fileext = ".sam")
cat(paste(sam, collapse="\n"), file=samfile)
bamfile <- Rsamtools::asBam(samfile, indexDestination = TRUE)
## classes --------------------------------------------------------------------

## IdeogramTrack
ideoTrack <- IdeogramTrack(chromosome="chrI", genome="sacCer3", bands=cyto.bands)

## GenomeAxisTrack
axisTrack <- GenomeAxisTrack(gr)

## DataTrack
dataTrack <- DataTrack(gr)

## AnnotationTrack
annoTrack <- AnnotationTrack(gr)

## GeneRegionTrack

## BiomartGeneRegionTrack

## DetailsAnnotationTrack

## SequenceTrack
seqTrack.dna <- SequenceTrack(dna.sq)
seqTrack.rna <- RNASequenceTrack(rna.sq)
seqTrack.bs <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg19)

## AlignmentsTrack

## HighlightTrack
highTrack <- HighlightTrack(dataTrack, ranges=gr, chromosome="chr1")

## OverlayTrack
overTrack <- OverlayTrack(c(dataTrack, dataTrack))