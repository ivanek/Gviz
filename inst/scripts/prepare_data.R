library(tidyverse)
library(Gviz)
library(ensembldb)

# mouse genome
mm_genome <- "mm39"

# human genome
hs_genome <- "hg38"


# cpgIslands -----
my_url <- str_c(
    "https://hgdownload.gi.ucsc.edu/goldenPath/",
    hs_genome,
    "/database/cpgIslandExt.txt.gz"
)

cpgIslands <- read_tsv(gzcon(url(my_url)), col_names = FALSE) |>
    dplyr::filter(X2 == "chr7") |>
    dplyr::mutate(X3 = X3 + 1) |>
    dplyr::select(seqname = X2, start = X3, end = X4) |>
    top_n(n = 10) |>
    as(Class = "GRanges")
genome(cpgIslands) <- "chr38"

save(cpgIslands, file = "data/cpgIslands.rda")

# itrack -----
itrack <- IdeogramTrack(chromosome = "chr7", genome = hs_genome)
save(itrack, file = "data/itrack.rda")

# geneModels -----
# need to fix EnsDb
edb <- EnsDb("~/Downloads/hg38_ensdb112_toUCSC.sqlite")
mycolumns <- c("gene_biotype", "gene_id", "tx_id", "exon_id", "symbol")
myfilters <- GRangesFilter(range(cpgIslands))

# myfilters <- GeneIdFilter(c(
#     "ENSG00000005020",
#     "ENSG00000122548",
#     "ENSG00000213787",
#     "ENSG00000222004",
#     "ENSG00000226059",
#     "ENSG00000233760"
# ))

exons_all <- exonsBy(
    edb,
    by = "tx",
    columns = mycolumns,
    filter = myfilters
) |>
    unlist()

cds <- cdsBy(
    edb,
    by = "tx",
    columns = mycolumns,
    filter = myfilters
) |>
    unlist()

utr5 <- fiveUTRsByTranscript(
    edb,
    columns = mycolumns,
    filter = myfilters
) |>
    unlist()
utr5$gene_biotype <- "utr5"

utr3 <- threeUTRsByTranscript(
    edb,
    columns = mycolumns,
    filter = myfilters
) |>
    unlist()
utr3$gene_biotype <- "utr3"

exons_pc <- c(cds, utr5, utr3)
#ov <- findOverlaps(exons_all, exons_pc)
#exons_all <- exons_all[-unique(queryHits(ov))]
geneModels <- c(
    exons_pc,
    exons_all[-which(exons_all$exon_id %in% exons_pc$exon_id)]
) |>
    sort()
names(geneModels) <- NULL
geneModels <- geneModels |>
    as.data.frame() |>
    dplyr::select(
        chromosome = seqnames,
        start,
        end,
        width,
        strand,
        feature = gene_biotype,
        gene = gene_id,
        exon = exon_id,
        transcript = tx_id,
        symbol
    )
save(geneModels, file = "data/geneModels.rda")

# ideoTrack -----
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chrX")
save(ideoTrack, file = "data/ideoTrack.rda")


# denseAnnTrack -----

data("denseAnnTrack")
dp <- displayPars(denseAnnTrack)
denseAnn <- exons(edb, columns = mycolumns,
                  filter = GRangesFilter(GRanges("chr7", IRanges(1.0019e8, 1.00e8))))
                  #filter = SeqNameFilter("chr7"))
denseAnn <- denseAnn[denseAnn$gene_biotype == "protein_coding"]
denseAnn$group <- denseAnn$tx_id
denseAnn$id <- denseAnn$exon_id
denseAnnTrack <- AnnotationTrack(denseAnn)
displayPars(denseAnnTrack) <- dp
save(denseAnnTrack, file = "data/denseAnnTrack.rda")

