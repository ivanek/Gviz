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
    expect_identical(
        vapply(
            highTrack@trackList,
            function(x) displayPars(x)$col,
            FUN.VALUE = character(1)
        ),
        "red"
    )
})
