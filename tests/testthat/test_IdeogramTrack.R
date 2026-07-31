test_that("IdeogramTrack / .cacheGenomes uses an https UCSC session for hg38", {
    check_ucsc()

    ideo <- IdeogramTrack(genome = "hg38", chromosome = "chr12")
    expect_s4_class(ideo, "IdeogramTrack")
    expect_s4_class(ideoTrack, "GdObject")
    expect_s4_class(ideoTrack, "RangeTrack")

    ## .cacheGenomes() caches the session under "session_<genome>" in the
    ## same .ucscCache environment used by .cacheTracks(), keyed off
    ## .gvizUcscUrl() at creation time.
    session <- get("session_hg38", envir = Gviz:::.ucscCache)
    expect_true(startsWith(session@url, "https://"))
})
