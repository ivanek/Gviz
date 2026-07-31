test_that("AlignmentsTrack works", {
    expect_s4_class(AlignmentsTrack(), "GdObject")
    expect_s4_class(AlignmentsTrack(), "RangeTrack")
    expect_s4_class(AlignmentsTrack(), "StackedTrack")
    expect_s4_class(AlignmentsTrack(), "AlignmentsTrack")
    expect_s4_class(AlignmentsTrack(showIndels = TRUE), "AlignmentsTrack")

    expect_s4_class(AlignmentsTrack(bamfile), "ReferenceTrack")
    expect_s4_class(AlignmentsTrack(bamfile), "AlignmentsTrack")
})
