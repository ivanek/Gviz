library(Gviz)

test_that("checking the class and structure works", {
    expect_null(Gviz:::.checkClass(data.frame(), "data.frame"))
    expect_error(Gviz:::.checkClass(data.frame(), "matrix"), "of class")
}
)
