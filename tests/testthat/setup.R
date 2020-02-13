## libs -----------------------------------------------------------------------


library(BSgenome.Celegans.UCSC.ce2)

## general --------------------------------------------------------------------
set.seed(789)
ir <- IRanges(1L,5L)
gr <- GRanges("chr1", ir, score=1)
ir2 <- IRanges(2L,6L)
gr2 <- GRanges("chr1", ir2, score=2)
dna.sq <- DNAStringSet(c(chr1=paste(sample(DNA_BASES, 100, replace=T), collapse="")))
rna.sq <- RNAStringSet(c(chr1=paste(sample(RNA_BASES, 100, replace=T), collapse="")))
cyto.bands <- data.frame(chrom = rep("chrI", 4),
                         chromStart = c(1, 148071, 151524, 154977),
                         chromEnd = c(148071, 151524, 154977, 230218),
                         name = c(NA, "CEN1", "CEN1", NA),
                         gieStain = c("gneg", "acen", "acen", "gneg"))

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
seqTrack.rna <- SequenceTrack(rna.sq)
seqTrack.bs <- SequenceTrack(BSgenome.Celegans.UCSC.ce2)

## AlignmentsTrack

## HighlightTrack
highTrack <- HighlightTrack(dataTrack, ranges=gr, chromosome="chr1")

## OverlayTrack
overTrack <- OverlayTrack(c(dataTrack, dataTrack))