#' Univariate kernel density estimation for bounded and unbounded support
#'
#' @param data vector of length \eqn{n}.
#' @param xmin lower bound for the support of the density.
#' @param xmax upper bound for the support of the density.
#' @param bw bandwidth parameter; has to be a positive number or \code{NULL};
#' the latter calls an automatic selection routine.
#' @param mult numeric; the actual bandwidth used is \eqn{bw*mult}.
#'
#' @return An object of class \code{kde1d}.
#'
#' @details
#' If \code{xmin} or \code{xmax} are finite, the density estimate will be 0
#' outside of \eqn{[xmin, xmax]}. Mirror-reflection is used to correct for
#' boundary bias.
#'
#' @seealso
#' \code{\link{dkde1d}},
#' \code{\link{pkde1d}},
#' \code{\link{qkde1d}},
#' \code{\link{rkde1d}}
#' \code{\link{plot.kde1d}} ,
#' \code{\link{lines.kde1d}}
#'
#' @examples
#' data(wdbc, package = "kdecopula")  # load data
#' fit <- kde1d(wdbc[, 5])            # estimate density
#' dkde1d(1000, fit)                  # evaluate density estimate
#'
#' @importFrom ks hpi
#' @export
kde1d <- function(data, xmin = -Inf, xmax = Inf, bw = NULL, mult = 1) {
    ## check/complete function call
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
        if (any(data < xmin) || any(data > xmax))
            stop("Not all data are contained in the interval [xmin, xmax].")
    } else if (!is.nan(xmin)) {
        if(any(data < xmin))
            stop("Not all data are larger than xmin.")
    } else if (!is.nan(xmax)) {
        if(any(data > xmax))
            stop("Not all data are samller than xmax")
    }

    ## bandwidth selection
    if (is.null(bw))
        bw <- NA
    if (is.na(bw))
        bw <- hpi(data)

    ## return kde1d object
    res <- list(data = data,
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
#' dkde1d(1000, fit)  # evaluate density estimate
#' pkde1d(1000, fit)  # evaluate corresponding cdf
#'
#' @useDynLib kdevine
#' @import Rcpp
#' @export
dkde1d <- function(x, obj) {
    eval_kde1d(sort(obj$data), x, obj$xmin, obj$xmax, obj$bw)
}

#' @rdname dkde1d
#' @export
pkde1d <- function(x, obj) {
    eval_pkde1d(obj$data, x, obj$xmin, obj$xmax, obj$bw)
}

#' @rdname dkde1d
#' @useDynLib kdevine
#' @export
qkde1d <- function(x, obj) {
    stopifnot(all((x >= 0) & (x <= 1)))
    eval_qkde1d(obj$data, x, obj$xmin, obj$xmax, obj$bw)
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
    qkde1d(U)
}

#' Plotting kde1d objects
#'
#' @aliases lines.kde1d
#' @method plot kde1d
#'
#' @param x \code{kde1d} object.
#' @param ev gridpoints for the plot (optional).
#' @param ... further arguments passed to \code{\link{plot.default}}.
#'
#' @seealso
#' \code{\link{kde1d}}
#' \code{\link{lines.kde1d}}
#'
#' @examples
#' data(wdbc)  # load data
#' fit <- kde1d(wdbc[, 6])  # estimate density
#' plot(fit)  # plot density estimate
#'
#' fit2 <- kde1d(wdbc[, 7])  # estimate density for another variable
#' lines(fit2, col = 2)  # add second estimate to the plot
#'
#' @importFrom graphics plot
#' @importFrom utils modifyList
#' @export
plot.kde1d <- function(x, ev = NULL, ...) {
    if (is.null(ev)) {
        p.l <- if (is.nan(x$xmin)) min(x$data) - x$bw else x$xmin
        p.u <- if (is.nan(x$xmax)) max(x$data) + x$bw else x$xmax
        ev <- seq(p.l, p.u, l = 100)
    }
    fhat <- dkde1d(ev, x)

    pars <- list(x = ev,
                 y = fhat,
                 type = "l",
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
lines.kde1d <- function(x, ev = NULL, ...) {
    if (is.null(ev)) {
        p.l <- if (is.nan(x$xmin)) min(x$data) - x$bw else x$xmin
        p.u <- if (is.nan(x$xmax)) max(x$data) + x$bw else x$xmax
        ev <- seq(p.l, p.u, l = 100)
    }
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
#'             lowr <- if (is.nan(xmin)) min(data) - bw else xmin
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
