test_that("SequenceTrack works", {
    expect_s4_class(SequenceTrack(), "GdObject")
    expect_s4_class(SequenceTrack(), "SequenceTrack")
    expect_s4_class(SequenceTrack(), "SequenceDNAStringSetTrack")
    expect_s4_class(RNASequenceTrack(), "SequenceRNAStringSetTrack")
})
