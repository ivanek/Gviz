test_that("GeneRegionTrack works", {
    expect_s4_class(GeneRegionTrack(), "AnnotationTrack")
    expect_s4_class(GeneRegionTrack(), "GeneRegionTrack")
})
