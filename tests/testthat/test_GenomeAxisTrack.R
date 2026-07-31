test_that("GenomeAxisTrack works", {
    expect_s4_class(GenomeAxisTrack(), "GdObject")
    expect_s4_class(GenomeAxisTrack(), "GenomeAxisTrack")
})
