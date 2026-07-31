test_that("DisplayPars works", {
    dp <- DisplayPars(a = 1)

    expect_s4_class(dp, "DisplayPars")
    expect_identical(dp@pars$a, 1)
    expect_error(DisplayPars(1), "All supplied arguments must be named.")
    expect_error(DisplayPars(1, b = 2), "All supplied arguments must be named.")
    # .updateDp
    expect_message(.updateDp(dp), "Note that the behaviour of the")
    expect_identical(.updateDp(dp, interactive = FALSE)@pars, list(a = 1))
    dp@pars <- new.env()
    dp@pars[["a"]] <- 1
    expect_identical(.updateDp(dp, interactive = FALSE)@pars, list(a = 1))
    expect_message(
        .updateDp(dp, interactive = FALSE),
        "The DisplayPars object has been updated"
    )
    ## getPar
    expect_identical(getPar(dp), list(a = 1))
    expect_identical(getPar(dp, "a"), 1)
    ## setPar
    expect_identical(setPar(dp, "b", 2, interactive = FALSE)@pars$b, 2)
    expect_error(setPar(dp, "b", c(1, 2), interactive = FALSE), "equal length")
    ## displayPars <-
    displayPars(dp) <- list(b = 2)
    expect_identical(dp@pars, list(a = 1, b = 2))
    ## getPar
    expect_identical(getPar(dp), list(a = 1, b = 2))
    expect_identical(getPar(dp, "a"), 1)
    ## displayPars
    expect_identical(displayPars(dp), list(a = 1, b = 2))
    expect_identical(displayPars(dp, "a"), 1)
    ## as.list
    expect_identical(as.list(dp), list(a = 1, b = 2))
})
