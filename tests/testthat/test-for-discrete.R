context("Discrete data")

dat <- data.frame(
    f1 = as.factor(rbinom(10, 1, 0.5)),
    x1 = rnorm(10),
    z1 = as.ordered(rpois(10, 1))
)

test_that("allows ordered input", {
    fit <- kde1d(dat$z1)
    expect_is(dkde1d(data.frame(dat$z1), fit), "numeric")
    expect_is(pkde1d(data.frame(dat$z1), fit), "numeric")

    expect_is(qkde1d(runif(10), fit), "ordered")
    fit <- kde1d(dat$x1)
    expect_is(qkde1d(runif(10), fit), "numeric")

    fit <- kdevine(dat[, 2:3])
    expect_is(dkdevine(dat[, 2:3], fit), "numeric")
})


test_that("kde1d does not allow factor input", {
    expect_error(fit <- kde1d(dat$f1))
    fit <- kde1d(dat$z1)
    expect_error(dkde1d(dat$f1, fit))
    expect_error(dkde1d(data.frame(dat$f1), fit))
})
