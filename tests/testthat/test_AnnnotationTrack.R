## AnnotationTrack

test_that("AnnnotationTrack works", {
    expect_s4_class(AnnotationTrack(), "GdObject")
    expect_s4_class(AnnotationTrack(), "RangeTrack")
    expect_s4_class(AnnotationTrack(), "StackedTrack")
    expect_s4_class(AnnotationTrack(), "AnnotationTrack")

    expect_error(
        AnnotationTrack(stacking = "allover"),
        "following values for 'stacking'"
    )
})

## DetailsAnnotationTrack

test_that("DetailsAnnotationTrack works", {
    expect_s4_class(DetailsAnnotationTrack(), "GdObject")
    expect_s4_class(DetailsAnnotationTrack(), "StackedTrack")
    expect_s4_class(DetailsAnnotationTrack(), "RangeTrack")
    expect_s4_class(DetailsAnnotationTrack(), "AnnotationTrack")

    expect_s4_class(detTrack, "AnnotationTrack")
    expect_identical(detTrack@selectFun, selFun)
    expect_identical(detTrack@fun, detFun)
})
