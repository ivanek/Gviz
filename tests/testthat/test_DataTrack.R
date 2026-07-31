test_that("DataTrack works", {
    expect_s4_class(DataTrack(), "GdObject")
    expect_s4_class(DataTrack(), "RangeTrack")
    expect_s4_class(DataTrack(), "NumericTrack")
    expect_s4_class(DataTrack(), "DataTrack")
})
