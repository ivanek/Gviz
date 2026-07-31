## ImageMap

test_that("ImageMap works", {
    mat <- matrix(1, ncol = 4, dimnames = list("a", NULL))
    tags <- list(a = c(a = "tag"))
    expect_s4_class(ImageMap(mat, tags = tags), "ImageMap")
    expect_identical(ImageMap(mat, tags = tags)@coords, mat)
    expect_identical(ImageMap(mat, tags = tags)@tags, tags)
    expect_error(
        ImageMap(matrix(1, ncol = 3), tags = tags),
        "must be a numeric matrix with 4 column"
    )
    expect_error(
        ImageMap(matrix(1, ncol = 4), tags = tags),
        "Rownames must be set for the matrix in"
    )
    expect_error(
        ImageMap(mat, tags = list(c(a = "tag1"))),
        "must be a named list with character vector items."
    )
    expect_error(
        ImageMap(mat, tags = list(a = c(a = "tag1", b = "tag2"))),
        "following values in the"
    )
    expect_error(
        ImageMap(mat, tags = list(a = c("tag1"))),
        "items in the 'tags' list must be named character"
    )
})
