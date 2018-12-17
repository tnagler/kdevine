#' Univariate kernel density estimation for bounded and unbounded support
#'
#' Discrete variables are convoluted with the uniform distribution (see, Nagler,
#' 2017). If a variable should be treated as discrete, declare it as
#' [ordered()].
#'
#' @param x vector of length \eqn{n}.
#' @param mult numeric; the actual bandwidth used is \eqn{bw*mult}.
#' @param xmin lower bound for the support of the density.
#' @param xmax upper bound for the support of the density.
#' @param bw bandwidth parameter; has to be a positive number or \code{NULL};
#'   the latter calls [`KernSmooth::dpik()`].
#' @param bw_min minimum value for the bandwidth.
#' @param ... unused.
#'
#' @return An object of class \code{kde1d}.
#'
#' @details If \code{xmin} or \code{xmax} are finite, the density estimate will
#'   be 0 outside of \eqn{[xmin, xmax]}. Mirror-reflection is used to correct
#'   for boundary bias. Discrete variables are convoluted with the uniform
#'   distribution (see, Nagler, 2017).
#'
#' @seealso \code{\link{dkde1d}}, \code{\link{pkde1d}}, \code{\link{qkde1d}},
#'   \code{\link{rkde1d}} \code{\link{plot.kde1d}} , \code{\link{lines.kde1d}}
#'
#' @references Nagler, T. (2017). *A generic approach to nonparametric function
#'   estimation with mixed data.* [arXiv:1704.07457](https://arxiv.org/abs/1704.07457)
#'
#' @examples
#' data(wdbc, package = "kdecopula")  # load data
#' fit <- kde1d(wdbc[, 5])            # estimate density
#' dkde1d(1000, fit)                  # evaluate density estimate
#'
#' @importFrom KernSmooth dpik
#' @importFrom MASS bandwidth.nrd
#' @importFrom cctools cont_conv
#' @export
kde1d <- function(x, mult = 1, xmin = -Inf, xmax = Inf, bw = NULL, bw_min = 0, ...) {
    if (missing(x) & !is.null(list(...)$data))  # for backwards compatibility
        x <- list(...)$data

    ## check/complete function call
    stopifnot(NCOL(x) == 1)
    if (is.null(xmin))
        xmin <- NaN
    if (is.null(xmax))
        xmax <- NaN
    if (!is.finite(xmin))
        xmin <- NaN
    if (!is.finite(xmax))
        xmax <- NaN
    if (!is.nan(xmax) & !is.nan(xmin)) {
        if (xmin > xmax)
            stop("'xmin' is larger than 'xmax'")
        if (any(x < xmin) || any(x > xmax))
            stop("Not all data are contained in the interval [xmin, xmax].")
    } else if (!is.nan(xmin)) {
        if (any(x < xmin))
            stop("Not all data are larger than xmin.")
    } else if (!is.nan(xmax)) {
        if (any(x > xmax))
            stop("Not all data are samller than xmax")
    }

    ## make continuous if discrete
    if (!is.ordered(x) & is.factor(x))
        stop("Factors not allowed; use kdevine() or cctools::cckde().")
    lvls <- levels(x)
    x_cc <- cctools::cont_conv(x, ...)

    ## bandwidth selection
    if (is.null(bw))
        bw <- NA
    if (is.na(bw)) {
        # plug in method
        bw <- try(KernSmooth::dpik(x_cc))
        # if it fails: normal rule of thumb
        if (inherits(bw, "try-error"))
            bw <- MASS::bandwidth.nrd(x_cc)
        # if everything fails:
        if (all(bw == 0))
            bw <- 1
    }
    # bound from below
    if (length(attr(x_cc, "i_disc")) == 1) {
        # for discrete use 1 - theta
        bw <- max(bw, 0.5 - attr(x_cc, "theta"))
    } else {
        # for continuous use bw_min
        bw <- max(bw, bw_min)
    }

    ## return kde1d object
    res <- list(x_cc = x_cc,
                levels = lvls,
                xmin = xmin,
                xmax = xmax,
                bw   = bw * mult)
    class(res) <- "kde1d"
    res
}


#' Working with a kde1d object
#'
#' The density, cdf, or quantile function of a kernel density estimate are
#' evaluated at arbitrary points with \code{\link{dkde1d}}, \code{\link{pkde1d}},
#' and \code{\link{qkde1d}} respectively.
#'
#' @aliases pkde1d, qkde1d, rkde1d
#'
#' @param x vector of evaluation points.
#' @param obj a \code{kde1d} object.
#'
#'
#' @return The density or cdf estimate evaluated at \code{x}.
#'
#' @seealso
#' \code{\link{kde1d}}
#'
#' @examples
#' data(wdbc)  # load data
#' fit <- kde1d(wdbc[, 5])  # estimate density
#' dkde1d(1000, fit)        # evaluate density estimate
#' pkde1d(1000, fit)        # evaluate corresponding cdf
#' qkde1d(0.5, fit)         # quantile function
#' hist(rkde1d(100, fit))   # simulate
#'
#' @useDynLib kdevine
#' @importFrom Rcpp evalCpp
#' @importFrom cctools expand_as_numeric
#' @export
dkde1d <- function(x, obj) {
    if (is.data.frame(x))
        x <- x[[1]]
    if (!is.ordered(x))
        stopifnot(!is.factor(x))
    x <- cctools::expand_as_numeric(x)
    f <- eval_kde1d(sort(obj$x_cc), x, obj$xmin, obj$xmax, obj$bw)
    if (length(attr(obj$x_cc, "i_disc") == 1)) {
        # for discrete variables we can normalize
        x_all_num <- expand_as_numeric(as.ordered(obj$levels))
        f_all <- eval_kde1d(sort(obj$x_cc), x_all_num, obj$xmin, obj$xmax, obj$bw)
        f <- f / sum(f_all)
    }

    f
}

#' @rdname dkde1d
#' @export
pkde1d <- function(x, obj) {
    if (is.data.frame(x))
        x <- x[[1]]
    if (!is.ordered(x))
        stopifnot(!is.factor(x))
    x <- cctools::expand_as_numeric(x)
    if (length(attr(obj$x_cc, "i_disc") == 1)) {
        # for discrete variables we have to add the missing probability mass
        x_all_num <- expand_as_numeric(as.ordered(obj$levels))
        f_all <- dkde1d(x_all_num, obj)
        p <- sapply(x, function(y) sum(f_all[x_all_num <= y]))
    } else {
        p <- eval_pkde1d(sort(obj$x_cc), x, obj$xmin, obj$xmax, obj$bw)
    }

    p
}

#' @rdname dkde1d
#' @export
qkde1d <- function(x, obj) {
    if (is.data.frame(x))
        x <- x[[1]]
    stopifnot(all((x >= 0) & (x <= 1)))
    q <- eval_qkde1d(sort(obj$x_cc), x, obj$xmin, obj$xmax, obj$bw)

    ## for discrete variables compute quantile from the density
    if (length(attr(obj$x_cc, "i_disc") == 1)) {
        x_all_num <- expand_as_numeric(as.ordered(obj$levels))

        # pdf at all possible values of x
        dd <- eval_kde1d(sort(obj$x_cc), x_all_num, obj$xmin, obj$xmax, obj$bw)
        pp <- c(cumsum(dd)) / sum(dd)

        # generalized inverse
        q <- x_all_num[vapply(x, function(y) which(y <= pp)[1], integer(1))]
        q <- ordered(obj$levels[q + 1], levels = obj$levels)
    }

    q
}

#' @param n integer; number of observations.
#' @param quasi logical; the default (\code{FALSE}) returns pseudo-random
#' numbers, use \code{TRUE} for quasi-random numbers (generalized Halton, see
#' \code{\link[qrng:ghalton]{ghalton}}).
#'
#' @rdname dkde1d
#' @importFrom qrng ghalton
#' @importFrom stats runif
#' @export
rkde1d <- function(n, obj, quasi = FALSE) {
    # simulate (psuedo/quasi) uniform random variables
    if (!quasi) {
        U <- runif(n)
    } else {
        U <- ghalton(n, d = 1)
    }
    # simulated data from KDE is the quantile transform of U
    qkde1d(U, obj)
}

#' Plotting kde1d objects
#'
#' @aliases lines.kde1d
#' @method plot kde1d
#'
#' @param x \code{kde1d} object.
#' @param ... further arguments passed to \code{\link{plot.default}}.
#'
#' @seealso
#' \code{\link{kde1d}}
#' \code{\link{lines.kde1d}}
#'
#' @examples
#' data(wdbc)  # load data
#' fit <- kde1d(wdbc[, 7])  # estimate density
#' plot(fit)  # plot density estimate
#'
#' fit2 <- kde1d(as.ordered(wdbc[, 1])) # discrete variable
#' plot(fit2, col = 2)
#'
#' @importFrom graphics plot
#' @importFrom utils modifyList
#' @export
plot.kde1d <- function(x, ...) {
    p.l <- if (is.nan(x$xmin)) min(x$x_cc) - x$bw else x$xmin
    p.u <- if (is.nan(x$xmax)) max(x$x_cc) + x$bw else x$xmax
    ev <- seq(p.l, p.u, l = 100)
    plot_type <- "l"  # for continuous variables, use a line plot
    if (length(attr(x$x_cc, "i_disc")) == 1) {
        ev <- as.ordered(x$levels)
        plot_type <- "h"  # for discrete variables, use a histrogram
    }
    fhat <- dkde1d(expand_as_numeric(ev), x)

    pars <- list(x = ev,
                 y = fhat,
                 type = plot_type,
                 xlab = "x",
                 ylab = "density",
                 ylim = c(0, 1.1 * max(fhat)))

    do.call(plot, modifyList(pars, list(...)))
}

#' @method lines kde1d
#'
#' @rdname plot.kde1d
#' @importFrom graphics lines
#' @importFrom utils modifyList
#' @export
lines.kde1d <- function(x, ...) {
    if (length(attr(x$x_cc, "i_disc") == 1))
        stop("lines does not work for discrete estimates.")
    p.l <- if (is.nan(x$xmin)) min(x$x_cc) - x$bw else x$xmin
    p.u <- if (is.nan(x$xmax)) max(x$x_cc) + x$bw else x$xmax
    ev <- seq(p.l, p.u, l = 100)

    fhat <- dkde1d(ev, x)

    pars <- list(x = ev, y = fhat)
    do.call(lines, modifyList(pars, list(...)))
}



#' #' Bandwidth selection for kde1d
#' #'
#' #' @param data data vector.
#' #' @param xmin lower bound for the support of the density.
#' #' @param xmax upper bound for the support of the density.
#' #' @param n.subs number of subsamples used to estimate the cross-validation
#' #' criterion.
#' #'
#' #' @importFrom stats integrate mad sd optimize
#' bw_kde1d <- function(data, xmin, xmax, n.subs) {
#'     n <- length(data)
#'     cK <- 1
#'     if (missing(xmin))
#'         xmin <- NaN
#'     if (missing(xmax))
#'         xmax <- NaN
#'     ## bandwidth selection
#'     if (is.nan(xmin) & is.nan(xmax)) {
#'         ## for unbounded data calll hpi from ks package
#'         bw <- hpi(data)
#'     } else {
#'         ## LSCV for bounded data
#'         subs <- sample(n.subs)
#'         M1 <- data[subs] %*% t(rep(1, n.subs))
#'         M2 <- rep(1, n.subs) %*% t(data[subs])
#'         difs <- M1 - M2
#'         sums <- M1 + M2
#'         diag(difs) <- rep(NA, n.subs)
#'         lsc <- function(bw) {
#'             K <- kern_gauss(difs / bw) / bw
#'             if (!is.nan(xmin))
#'                 K <- K + kern_gauss((2 * xmin - sums) / bw) / bw
#'             if (!is.nan(xmax))
#'                 K <- K + kern_gauss((2 * xmax - sums) / bw) / bw
#'             lowr <-    if (is.list(data))
#'                 data <- data[[1]]
#'             if (is.nan(xmin)) min(data) - bw else xmin
#'             upr  <- if (is.nan(xmax)) max(data) + bw else xmax
#'             integrand <- function(x) eval_kde1d(sort(data), x, xmin, xmax, bw)^2
#'
#'             tmp1 <- try(integrate(integrand, lowr, upr, rel.tol = 5e-2)$value,
#'                         silent = TRUE)
#'             tmp2 <- 2 * sum(K, na.rm = TRUE) / (n.subs*(n.subs-1))
#'
#'             if ("try-error" %in% class(tmp1)) Inf else tmp1 - tmp2
#'         }
#'         # starting parameter with normal reference rule
#'         bw.start <- cK * min(sd(data), mad(data)) * n^(-1/5)
#'         if (bw.start == 0)
#'             bw.start = 0.5
#'         opt <- suppressWarnings(optimize(lsc,
#'                                          lower = bw.start / 5,
#'                                          upper = bw.start * 5))
#'         bw <- opt$minimum * (n / n.subs)^(-1/5)
#'
#'         # correction for discreteness
#'         bw <- min(bw * (1 + log(n/length(levels(as.factor(data))))),
#'                   0.5 * diff(range(data)))
#'     }
#'
#'     ## return results
#'     bw
#' }
