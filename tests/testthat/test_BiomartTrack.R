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
        start = 20e6,
        end = 21e6,
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
