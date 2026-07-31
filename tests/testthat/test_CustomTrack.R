test_that("CustomTrack works", {
    expect_s4_class(CustomTrack(), "GdObject")
    expect_s4_class(CustomTrack(), "CustomTrack")
})
