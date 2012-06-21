library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx = exons(TxDb.Hsapiens.UCSC.hg19.knownGene, 
   vals=list(gene_id=3845), columns=c("exon_id", "tx_id"))
kr1 = as.data.frame(ranges(tx)) 
kr1$gene = "3845"
kr1$exon = values(tx)$exon_id
ntx = sapply(values(tx)$tx_id, length)
rep(1:8, ntx)
kr2 = kr1[rep(1:8,ntx),] 
kr2$transcript = unlist(values(tx)$tx_id)
bad = which(names(kr2) == "width")
if (length(bad) > 0) kr2 = kr2[,-bad]
kr2$rank = 1:14
library(Gviz)
krg = GeneRegionTrack(kr2, chrom=12, genome="hg19")
krg@name = "KRAS"
values(krg@range)$symbol = kr2$transcript
displayPars(krg)$showFeatureId = TRUE
displayPars(krg)$background.title = "black"
plotTracks(list(IdeogramTrack(genome="hg19", chr="chr12"), 
    GenomeAxisTrack(), krg),
    showId=TRUE, extend.left=3000, showTitle=TRUE)
