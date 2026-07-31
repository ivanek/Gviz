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
    expect_identical(
        vapply(
            overTrack@trackList,
            function(x) displayPars(x)$col,
            FUN.VALUE = character(1)
        ),
        c("red", "red")
    )
})
