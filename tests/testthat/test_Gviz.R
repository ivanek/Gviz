library(Gviz)

test_that("checking the class and structure works", {
    expect_null(Gviz:::.checkClass(data.frame(), "data.frame"))
    expect_error(Gviz:::.checkClass(data.frame(), "matrix"), "of class")
    expect_error(Gviz:::.checkClass(class="matrix"), "is missing with no default")
    expect_error(Gviz:::.checkClass(data.frame(), "data.frame", length=2), "of length")
}
)
