#' Kernel density estimatior based on simplified vine copulas
#'
#' Implements the vine-copula based estimator of Nagler and Czado (2016). The
#' marginal densities are estimated by \code{\link{kde1d}}, the vine copula
#' density by \code{\link{kdevinecop}}. Discrete variables are convoluted with
#' the uniform distribution (see, Nagler, 2017). If a variable should be treated
#' as discrete, declare it as [ordered()]. Factors are expanded into binary
#' dummy codes.
#'
#' @param x (\eqn{n x d}) data matrix.
#' @param mult_1d numeric; all bandwidhts for marginal kernel density estimation
#'   are multiplied with \code{mult_1d}. Defaults to `log(1 + d)` where `d` is
#'   the number of variables after applying [cctools::expand_as_numeric()].
#' @param xmin numeric vector of length d; see \code{\link{kde1d}}.
#' @param xmax numeric vector of length d; see \code{\link{kde1d}}.
#' @param copula.type either \code{"kde"} (default) or \code{"parametric"} for
#'   kernel or parametric estimation of the vine copula.
#' @param ... further arguments passed to \code{\link{kde1d}} or
#'   \code{\link{kdevinecop}}.
#'
#' @return An object of class \code{kdevine}.
#'
#' @seealso \code{\link{dkdevine}} \code{\link{kde1d}} \code{\link{kdevinecop}}
#'
#' @references Nagler, T., Czado, C. (2016) *Evading the curse of
#'   dimensionality in nonparametric density estimation with simplified vine
#'   copulas.* Journal of Multivariate Analysis 151, 69-89
#'   (doi:10.1016/j.jmva.2016.07.003) \cr \cr
#'   Nagler, T. (2017). *A generic approach to nonparametric function
#'   estimation with mixed data.* [arXiv:1704.07457](https://arxiv.org/abs/1704.07457)

#'
#' @examples
#' # load data
#' data(wdbc, package = "kdecopula")
#' \dontshow{wdbc <- wdbc[1:30, ]}
#' # estimate density (use xmin to indicate positive support)
#' fit <- kdevine(wdbc[, 5:7], xmin = rep(0, 3))
#'
#' # evaluate density estimate
#' dkdevine(c(1000, 0.1, 0.1), fit)
#'
#' # plot simulated data
#' pairs(rkdevine(nrow(wdbc), fit))
#'
#' @importFrom VineCopula RVineStructureSelect RVineCopSelect
#' @importFrom cctools expand_vec
#' @export
kdevine <- function(x, mult_1d = NULL, xmin = NULL,
                    xmax = NULL, copula.type = "kde", ...) {
    if (missing(x) & !is.null(list(...)$data))  # for backwards compatibilitiy
        x <- list(...)$data
    if (!is.null(list(...)$mult.1d))  # for backwards compatibilitiy
        mult_1d <- list(...)$mult.1d
    x_cc <- cont_conv(x)
    if (NCOL(x_cc) == 1)
        stop("x must be multivariate or a factor.")
    d <- ncol(x_cc)

    ## sanity checks
    if (!is.null(xmin)) {
        xmin <- cctools::expand_vec(xmin, x)
    }
    if (!is.null(xmax)) {
        xmax <- cctools::expand_vec(xmax, x)
    }
    bw <- list(...)$bw
    if (!is.null(bw)) {
        bw <- cctools::expand_vec(bw, x)
    }

    ## estimation of the marginals
    i_disc <- attr(x_cc, "i_disc")
    if (is.null(mult_1d))
        mult_1d <- log(1 + ncol(x_cc))
    marg.dens <- as.list(numeric(d))
    for (k in 1:d) {
        marg.dens[[k]] <- kde1d(
            x_cc[, k],
            xmin = xmin[k],
            xmax = xmax[k],
            bw   = bw[k],
            mult = mult_1d,
            bw_min = ifelse(k %in% i_disc, 0.5 - attr(x_cc, "theta"), 0)
        )
    }
    res <- list(x_cc = x_cc, marg.dens = marg.dens)

    ## estimation of the R-vine copula (only if d > 1)
    if (d > 1) {
        # transform to copula data
        u <- sapply(1:d, function(k) pkde1d(x_cc[, k], marg.dens[[k]]))

        if (copula.type == "kde") {
            res$vine  <- suppressWarnings(
                kdevinecop(
                    u,
                    matrix      = list(...)$matrix,
                    method      = list(...)$method,
                    mult        = list(...)$mult,
                    info        = list(...)$info,
                    test.level  = list(...)$test.level,
                    trunc.level = list(...)$trunc.level,
                    treecrit    = list(...)$treecrit,
                    cores       = list(...)$cores
                )
            )
        } else if (copula.type == "parametric") {
            # get family and matrix if available
            fam <- list(...)$familyset
            if (is.null(fam))
                fam <- NA
            mat <- list(...)$Matrix

            # fit parametric vine
            res$vine <- if (is.null(mat) & d > 2) {
                # structure selection if no matrix is provided and d > 2
                RVineStructureSelect(u, familyset = fam)
            } else if (d == 2) {
                # select copula for default structure if d = 2
                RVineCopSelect(u,
                               familyset = fam,
                               Matrix = matrix(c(2, 1, 0, 1), 2, 2))
            } else {
                # or select copulas to provided structure
                RVineCopSelect(u,
                               familyset = fam,
                               Matrix = mat)
            }
        } else {
            stop("copula.type not implemented.")
        }
    }

    ## return results
    res$copula.type <- copula.type
    class(res) <- "kdevine"
    res
}

#' Evaluate the density of a kdevine object
#'
#' @param x (\eqn{m x d}) matrix of evaluation points (or vector of length \eqn{d}).
#' @param obj a \code{kdevine} object.
#'
#' @return The density estimate evaluated at \code{x}.
#'
#' @seealso
#' \code{\link{kdevine}}
#'
#' @examples
#' # load data
#' data(wdbc)
#' \dontshow{wdbc <- wdbc[1:30, ]}
#' # estimate density (use xmin to indicate positive support)
#' fit <- kdevine(wdbc[, 5:7], xmin = rep(0, 3))
#'
#' # evaluate density estimate
#' dkdevine(c(1000, 0.1, 0.1), fit)
#'
#' @importFrom VineCopula RVinePDF
#'
#' @export
dkdevine <- function(x, obj) {
    stopifnot(class(obj) == "kdevine")
    if (NCOL(x) == 1)
        x <- t(x)
    nms <- colnames(x)
    # must be numeric, factors are expanded
    x <- expand_as_numeric(x)
    # variables must be in same order
    if (!is.null(nms))
        x <- x[, colnames(obj$x_cc), drop = FALSE]

    ## evaluate marginals
    d <- ncol(x)
    margvals <- u <- x
    for (k in 1:d) {
        x_k <- x[, k]
        if (k %in% attr(obj$x_cc, "i_disc")) {
            # use normalization if discrete
            attr(x_k, "i_disc") <- 1
            obj$marg.dens[[k]]$levels <- attr(obj$x_cc, "levels")[[k]]
        }
        margvals[, k] <- dkde1d(x_k, obj$marg.dens[[k]])
    }

    if (!is.null(obj$vine)) {
        # PIT to copula level
        for (k in 1:d) {
            if (k %in% attr(obj$x_cc, "i_disc")) {
                # use continuous variant for PIT
                attr(x_k, "i_disc") <- integer(0)
                obj$marg.dens[[k]]$levels <- NULL
            }
            u[, k] <- pkde1d(x[, k], obj$marg.dens[[k]])
        }
        if (inherits(obj$vine, "kdevinecop")) {
            vinevals <- dkdevinecop(u, obj = obj$vine, stable = TRUE)
        } else if (inherits(obj$vine, "RVineMatrix")) {
            vinevals <- RVinePDF(u, obj$vine)
        } else {
            stop("vine has incompatible type")
        }
    } else {
        vinevals <- rep(1, nrow(x))
    }

    ## final density estimate is product of marginals and copula density
    apply(cbind(margvals, vinevals), 1, prod)
}


#' Simulate from a kdevine object
#'
#' @param n number of observations.
#' @param obj a \code{kdevine} object.
#' @param quasi logical; the default (\code{FALSE}) returns pseudo-random
#' numbers, use \code{TRUE} for quasi-random numbers (generalized Halton, only
#' works for fully nonparametric fits).
#'
#' @return An \eqn{n x d} matrix of simulated data from the \code{kdevine}
#' object.
#'
#' @seealso
#' \code{\link{kdevine}},
#' \code{\link{rkdevinecop}},
#' \code{\link{rkde1d}}
#'
#' @examples
#' # load and plot data
#' data(wdbc)
#' \dontshow{wdbc <- wdbc[1:30, ]}
#' # estimate density
#' fit <- kdevine(wdbc[, 5:7], xmin = rep(0, 3))
#'
#' # plot simulated data
#' pairs(rkdevine(nrow(wdbc), fit))
#'
#' @importFrom VineCopula pobs RVineSim
#' @export
rkdevine <- function(n, obj, quasi = FALSE) {
    # simulate from copula
    usim <- switch(obj$copula.type,
                   "kde" = rkdevinecop(n, obj$vine, quasi = quasi),
                   "parametric" = RVineSim(n, obj$vine))
    # use quantile transformation for marginals
    sapply(seq_len(ncol(usim)), function(i)
        qkde1d(usim[, i], obj$marg.dens[[i]])
    )
}




