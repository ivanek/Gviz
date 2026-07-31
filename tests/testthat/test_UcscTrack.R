## UcscTrack

test_that("the default UCSC URL is https", {
    expect_identical(Gviz:::.gvizUcscUrl(), "https://genome.ucsc.edu/cgi-bin/")
})

test_that(".flattenTableNames passes through an already-atomic vector", {
    x <- c(knownGene = "knownGene", cytoBand = "cytoBandIdeo")
    expect_identical(Gviz:::.flattenTableNames(x), x)
})

test_that(".flattenTableNames drops NULL (protected) entries", {
    tn <- list(
        knownGene = "knownGene",
        someProtectedTrack = NULL,
        cytoBand = "cytoBandIdeo"
    )
    flat <- Gviz:::.flattenTableNames(tn)

    expect_true(is.character(flat))
    expect_false("someProtectedTrack" %in% names(flat))
    expect_identical(
        unname(flat[c("knownGene", "cytoBand")]),
        c("knownGene", "cytoBandIdeo")
    )
})

test_that(".flattenTableNames replicates the track name across multi-table entries", {
    tn <- list(
        knownGene = "knownGene",
        geneHancer = c(
            "geneHancerClusteredInteractions",
            "geneHancerRegulatoryElements"
        )
    )
    flat <- Gviz:::.flattenTableNames(tn)

    expect_true(is.character(flat))
    expect_identical(
        unname(flat),
        c(
            "knownGene",
            "geneHancerClusteredInteractions",
            "geneHancerRegulatoryElements"
        )
    )
    expect_identical(names(flat), c("knownGene", "geneHancer", "geneHancer"))
})

test_that(".flattenTableNames returns an empty character vector when everything is protected", {
    tn <- list(someProtectedTrack = NULL, anotherOne = NULL)
    expect_identical(Gviz:::.flattenTableNames(tn), character(0))
})

test_that(".ucscTableQueryCompat reports an informative error on failure", {
    fakeQuery <- function(session, ...) stop("nope")
    testthat::local_mocked_bindings(
        ucscTableQuery = fakeQuery,
        .package = "Gviz"
    )

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

    expect_true(all(
        c(
            "exonStarts",
            "exonEnds",
            "exonCount",
            "txStart",
            "txEnd",
            "cdsStart",
            "cdsEnd"
        ) %in%
            colnames(out)
    ))

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
        chrom = "chr1",
        chromStart = 1L,
        chromEnd = 100L,
        blockSizes = "10,",
        chromStarts = "0,",
        exonStarts = "1,",
        exonEnds = "11,",
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
    expect_identical(
        Gviz:::.bigGenePredToGenePredCompat(data.frame()),
        data.frame()
    )
})

test_that("UcscTrack fetches a hg38 GeneRegionTrack from UCSC", {
    check_ucsc()

    ## Small region around GAPDH (hg38 coordinates) to keep the download small
    from <- 6534405
    to <- 6538375

    knownGenes <- UcscTrack(
        genome = "hg38",
        chromosome = "chr12",
        track = "knownGene",
        from = from,
        to = to,
        trackType = "GeneRegionTrack",
        rstarts = "exonStarts",
        rends = "exonEnds",
        gene = "name",
        symbol = "name",
        transcript = "name",
        strand = "strand",
        fill = "#8282d2",
        name = "UCSC Genes"
    )

    expect_s4_class(knownGenes, "GeneRegionTrack")
    expect_identical(unname(genome(knownGenes)), "hg38")
    expect_identical(unname(as.character(seqnames(knownGenes))[1]), "chr12")
    expect_true(length(knownGenes) > 0)
})

test_that("UcscTrack uses an https UCSC session for hg38", {
    check_ucsc()

    sessionInfo <- Gviz:::.cacheTracks(
        genome = "hg38",
        chromosome = "chr12",
        track = "knownGene"
    )
    expect_true(startsWith(sessionInfo$session@url, "https://"))
})
