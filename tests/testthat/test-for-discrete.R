context("Discrete data")

dat <- data.frame(
    f1 = as.factor(rbinom(10, 1, 0.5)),
    x1 = rnorm(10),
    z1 = as.ordered(rpois(10, 1))
)

test_that("allows ordered input", {
    fit <- kdevine(dat[, 2:3])
    expect_is(dkdevine(dat[, 2:3], fit), "numeric")
})
