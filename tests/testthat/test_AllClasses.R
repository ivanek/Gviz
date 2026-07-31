## Classes in Gviz

test_that("ImageMap works", {
    mat <- matrix(1, ncol = 4, dimnames = list("a", NULL))
    tags <- list(a = c(a = "tag"))
    expect_s4_class(ImageMap(mat, tags = tags), "ImageMap")
    expect_identical(ImageMap(mat, tags = tags)@coords, mat)
    expect_identical(ImageMap(mat, tags = tags)@tags, tags)
    expect_error(ImageMap(matrix(1, ncol = 3), tags = tags), "must be a numeric matrix with 4 column")
    expect_error(ImageMap(matrix(1, ncol = 4), tags = tags), "Rownames must be set for the matrix in")
    expect_error(ImageMap(mat, tags = list(c(a = "tag1"))), "must be a named list with character vector items.")
    expect_error(ImageMap(mat, tags = list(a = c(a = "tag1", b = "tag2"))), "following values in the")
    expect_error(ImageMap(mat, tags = list(a = c("tag1"))), "items in the 'tags' list must be named character")
})

test_that("DisplayPars works", {
    dp <- DisplayPars(a = 1)

    expect_s4_class(dp, "DisplayPars")
    expect_identical(dp@pars$a, 1)
    expect_error(DisplayPars(1), "All supplied arguments must be named.")
    expect_error(DisplayPars(1, b = 2), "All supplied arguments must be named.")
    # .updateDp
    expect_message(.updateDp(dp), "Note that the behaviour of the")
    expect_identical(.updateDp(dp, interactive = FALSE)@pars, list(a = 1))
    dp@pars <- new.env()
    dp@pars[["a"]] <- 1
    expect_identical(.updateDp(dp, interactive = FALSE)@pars, list(a = 1))
    expect_message(.updateDp(dp, interactive = FALSE), "The DisplayPars object has been updated")
    ## getPar
    expect_identical(getPar(dp), list(a = 1))
    expect_identical(getPar(dp, "a"), 1)
    ## setPar
    expect_identical(setPar(dp, "b", 2, interactive = FALSE)@pars$b, 2)
    expect_error(setPar(dp, "b", c(1, 2), interactive = FALSE), "equal length")
    ## displayPars <-
    displayPars(dp) <- list(b = 2)
    expect_identical(dp@pars, list(a = 1, b = 2))
    ## getPar
    expect_identical(getPar(dp), list(a = 1, b = 2))
    expect_identical(getPar(dp, "a"), 1)
    ## displayPars
    expect_identical(displayPars(dp), list(a = 1, b = 2))
    expect_identical(displayPars(dp, "a"), 1)
    ## as.list
    expect_identical(as.list(dp), list(a = 1, b = 2))
})

test_that("IdeogramTrack works", {
    ideoTrack <- IdeogramTrack(chromosome = "chrI", genome = "sacCer3", bands = cyto.bands)
    expect_s4_class(ideoTrack, "GdObject")
    expect_s4_class(ideoTrack, "RangeTrack")
    expect_s4_class(ideoTrack, "IdeogramTrack")
})

test_that("GenomeAxisTrack works", {
    expect_s4_class(GenomeAxisTrack(), "GdObject")
    expect_s4_class(GenomeAxisTrack(), "GenomeAxisTrack")
})

test_that("DataTrack works", {
    expect_s4_class(DataTrack(), "GdObject")
    expect_s4_class(DataTrack(), "RangeTrack")
    expect_s4_class(DataTrack(), "NumericTrack")
    expect_s4_class(DataTrack(), "DataTrack")
})

test_that("AnnnotationTrack works", {
    expect_s4_class(AnnotationTrack(), "GdObject")
    expect_s4_class(AnnotationTrack(), "RangeTrack")
    expect_s4_class(AnnotationTrack(), "StackedTrack")
    expect_s4_class(AnnotationTrack(), "AnnotationTrack")

    expect_error(AnnotationTrack(stacking = "allover"), "following values for 'stacking'")
})

test_that("GeneRegionTrack works", {
    expect_s4_class(GeneRegionTrack(), "AnnotationTrack")
    expect_s4_class(GeneRegionTrack(), "GeneRegionTrack")
})

test_that("BiomartGeneRegionTrack works", {
    ## expect_s4_class(BiomartGeneRegionTrack(), "BiomartGeneRegionTrack")

    biomartMapping <- list(
        gene_id = "ensembl_gene_id", transcript_id = "ensembl_transcript_id", exon_id = "ensembl_exon_id",
        start = "exon_chrom_start", end = "exon_chrom_end", rank = "rank", strand = "strand",
        symbol = c("external_gene_name", "external_gene_id"), feature = "gene_biotype", chromosome = "chromosome_name",
        u5s = "5_utr_start", u5e = "5_utr_end", u3s = "3_utr_start", u3e = "3_utr_end", cdsl = c("cds_length", "cds_start"),
        phase = "phase"
    )
    expect_identical(.getBMFeatureMap(), biomartMapping)
})

test_that("SequenceTrack works", {
    expect_s4_class(SequenceTrack(), "GdObject")
    expect_s4_class(SequenceTrack(), "SequenceTrack")
    expect_s4_class(SequenceTrack(), "SequenceDNAStringSetTrack")
    expect_s4_class(RNASequenceTrack(), "SequenceRNAStringSetTrack")
})

test_that("AlignmentsTrack works", {
    expect_s4_class(AlignmentsTrack(), "GdObject")
    expect_s4_class(AlignmentsTrack(), "RangeTrack")
    expect_s4_class(AlignmentsTrack(), "StackedTrack")
    expect_s4_class(AlignmentsTrack(), "AlignmentsTrack")
    expect_s4_class(AlignmentsTrack(showIndels = TRUE), "AlignmentsTrack")

    expect_s4_class(AlignmentsTrack(bamfile), "ReferenceTrack")
    expect_s4_class(AlignmentsTrack(bamfile), "AlignmentsTrack")
})

test_that("CustomTrack works", {
    expect_s4_class(CustomTrack(), "GdObject")
    expect_s4_class(CustomTrack(), "CustomTrack")
})

test_that("HighlightTrack works", {
    expect_s4_class(HighlightTrack(), "GdObject")
    expect_s4_class(HighlightTrack(), "RangeTrack")
    expect_s4_class(HighlightTrack(), "HighlightTrack")

    expect_true(is(HighlightTrack()@range, "GRanges"))
    expect_error(
        HighlightTrack(c(annoTrack, geneModels)),
        "All elements in 'trackList' must inherit from 'GdObject'"
    )

    displayPars(highTrack) <- list(col = "black")
    expect_identical(displayPars(highTrack)$col, "black")

    displayPars(highTrack, recursive = TRUE) <- list(col = "red")
    expect_identical(vapply(highTrack@trackList, function(x) displayPars(x)$col,
        FUN.VALUE = character(1)
    ), "red")
})

test_that("OverlayTrack works", {
    expect_s4_class(OverlayTrack(), "GdObject")
    expect_s4_class(OverlayTrack(), "OverlayTrack")

    expect_true(is.list(OverlayTrack(annoTrack)@trackList))
    expect_error(
        OverlayTrack(c(annoTrack, geneModels)),
        "All elements in 'trackList' must inherit from 'GdObject'"
    )

    displayPars(overTrack) <- list(col = "black")
    expect_identical(displayPars(overTrack)$col, "black")

    displayPars(overTrack, recursive = TRUE) <- list(col = "red")
    expect_identical(vapply(overTrack@trackList, function(x) displayPars(x)$col,
        FUN.VALUE = character(1)
    ), c("red", "red"))
})

test_that("DetailsAnnotationTrack works", {
    expect_s4_class(DetailsAnnotationTrack(), "GdObject")
    expect_s4_class(DetailsAnnotationTrack(), "StackedTrack")
    expect_s4_class(DetailsAnnotationTrack(), "RangeTrack")
    expect_s4_class(DetailsAnnotationTrack(), "AnnotationTrack")

    expect_s4_class(detTrack, "AnnotationTrack")
    expect_identical(detTrack@selectFun, selFun)
    expect_identical(detTrack@fun, detFun)
})

test_that("md5 hash for BiomartTrack works", {
    data("biomTrack")
    expect_identical(.bmGuid(biomTrack), "015394ebfc9e60c183d9bf4f4013c92a")
})

test_that("interaction with biomart works", {
    check_biomart()
    bm <- useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl"
    )
    biomartTrack <- BiomartGeneRegionTrack(
        chromosome = "chr7",
        start = 20e6, end = 21e6,
        name = "ENSEMBL",
        biomart = bm

    )
    expect_s4_class(biomartTrack, "GdObject")
    expect_s4_class(biomartTrack, "StackedTrack")
    expect_s4_class(biomartTrack, "RangeTrack")
    expect_s4_class(biomartTrack, "AnnotationTrack")
    expect_s4_class(biomartTrack, "GeneRegionTrack")
    expect_s4_class(biomartTrack, "BiomartGeneRegionTrack")
})


## UcscTrack --------------------------------------------------------------

test_that("the default UCSC URL is https", {
    expect_identical(Gviz:::.gvizUcscUrl(), "https://genome.ucsc.edu/cgi-bin/")
})

test_that(".flattenTableNames passes through an already-atomic vector", {
    x <- c(knownGene = "knownGene", cytoBand = "cytoBandIdeo")
    expect_identical(Gviz:::.flattenTableNames(x), x)
})

test_that(".flattenTableNames drops NULL (protected) entries", {
    tn <- list(knownGene = "knownGene", someProtectedTrack = NULL, cytoBand = "cytoBandIdeo")
    flat <- Gviz:::.flattenTableNames(tn)

    expect_true(is.character(flat))
    expect_false("someProtectedTrack" %in% names(flat))
    expect_identical(unname(flat[c("knownGene", "cytoBand")]), c("knownGene", "cytoBandIdeo"))
})

test_that(".flattenTableNames replicates the track name across multi-table entries", {
    tn <- list(knownGene = "knownGene", geneHancer = c("geneHancerClusteredInteractions", "geneHancerRegulatoryElements"))
    flat <- Gviz:::.flattenTableNames(tn)

    expect_true(is.character(flat))
    expect_identical(unname(flat), c("knownGene", "geneHancerClusteredInteractions", "geneHancerRegulatoryElements"))
    expect_identical(names(flat), c("knownGene", "geneHancer", "geneHancer"))
})

test_that(".flattenTableNames returns an empty character vector when everything is protected", {
    tn <- list(someProtectedTrack = NULL, anotherOne = NULL)
    expect_identical(Gviz:::.flattenTableNames(tn), character(0))
})

test_that(".ucscTableQueryCompat reports an informative error on failure", {
    fakeQuery <- function(session, ...) stop("nope")
    testthat::local_mocked_bindings(ucscTableQuery = fakeQuery, .package = "Gviz")

    expect_error(
        Gviz:::.ucscTableQueryCompat("session", "knownGene"),
        "Unable to query UCSC track/table 'knownGene'"
    )
})

test_that(".bigGenePredToGenePredCompat synthesizes exonStarts/exonEnds from bigGenePred columns", {
    ## Real row shape returned by UCSC's REST-backed hg38 knownGene track
    ## (GENCODE-based bigGenePred), from the bug report this is fixing.
    tableDat <- data.frame(
        chrom = "chr12",
        chromStart = 6534011L,
        chromEnd = 6538371L,
        name = "ENST00000920777.1",
        blockCount = 9L,
        blockSizes = "297,52,100,107,91,116,82,413,271,",
        chromStarts = "0,798,2482,2672,2908,3089,3297,3572,4089,",
        thickStart = 6534832L,
        thickEnd = 6538170L,
        stringsAsFactors = FALSE
    )

    out <- Gviz:::.bigGenePredToGenePredCompat(tableDat)

    expect_true(all(c("exonStarts", "exonEnds", "exonCount", "txStart", "txEnd", "cdsStart", "cdsEnd") %in% colnames(out)))

    starts <- as.integer(strsplit(sub(",$", "", out$exonStarts), ",")[[1]])
    ends <- as.integer(strsplit(sub(",$", "", out$exonEnds), ",")[[1]])

    expect_length(starts, 9)
    expect_length(ends, 9)
    expect_identical(starts[1], tableDat$chromStart)
    ## the last exon should end exactly at chromEnd
    expect_identical(ends[length(ends)], tableDat$chromEnd)
    expect_identical(out$exonCount, tableDat$blockCount)
    expect_identical(out$txStart, tableDat$chromStart)
    expect_identical(out$txEnd, tableDat$chromEnd)
    expect_identical(out$cdsStart, tableDat$thickStart)
    expect_identical(out$cdsEnd, tableDat$thickEnd)
})

test_that(".bigGenePredToGenePredCompat is a no-op when exonStarts/exonEnds already exist", {
    tableDat <- data.frame(
        chrom = "chr1", chromStart = 1L, chromEnd = 100L,
        blockSizes = "10,", chromStarts = "0,",
        exonStarts = "1,", exonEnds = "11,",
        stringsAsFactors = FALSE
    )
    expect_identical(Gviz:::.bigGenePredToGenePredCompat(tableDat), tableDat)
})

test_that(".bigGenePredToGenePredCompat is a no-op when bigGenePred columns are absent", {
    tableDat <- data.frame(chrom = "chr1", start = 1L, end = 100L)
    expect_identical(Gviz:::.bigGenePredToGenePredCompat(tableDat), tableDat)
})

test_that(".bigGenePredToGenePredCompat passes non-data.frame input through unchanged", {
    expect_identical(Gviz:::.bigGenePredToGenePredCompat(NULL), NULL)
    expect_identical(Gviz:::.bigGenePredToGenePredCompat(data.frame()), data.frame())
})

test_that("UcscTrack fetches a hg38 GeneRegionTrack from UCSC", {
    check_ucsc()

    ## Small region around GAPDH (hg38 coordinates) to keep the download small
    from <- 6534405
    to <- 6538375

    knownGenes <- UcscTrack(
        genome = "hg38", chromosome = "chr12", track = "knownGene",
        from = from, to = to, trackType = "GeneRegionTrack",
        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
        symbol = "name", transcript = "name", strand = "strand",
        fill = "#8282d2", name = "UCSC Genes"
    )

    expect_s4_class(knownGenes, "GeneRegionTrack")
    expect_identical(unname(genome(knownGenes)), "hg38")
    expect_identical(unname(as.character(seqnames(knownGenes))[1]), "chr12")
    expect_true(length(knownGenes) > 0)
})

test_that("UcscTrack uses an https UCSC session for hg38", {
    check_ucsc()

    sessionInfo <- Gviz:::.cacheTracks(genome = "hg38", chromosome = "chr12", track = "knownGene")
    expect_true(startsWith(sessionInfo$session@url, "https://"))
})

test_that("IdeogramTrack / .cacheGenomes uses an https UCSC session for hg38", {
    check_ucsc()

    ideo <- IdeogramTrack(genome = "hg38", chromosome = "chr12")
    expect_s4_class(ideo, "IdeogramTrack")

    ## .cacheGenomes() caches the session under "session_<genome>" in the
    ## same .ucscCache environment used by .cacheTracks(), keyed off
    ## .gvizUcscUrl() at creation time.
    session <- get("session_hg38", envir = Gviz:::.ucscCache)
    expect_true(startsWith(session@url, "https://"))
})
