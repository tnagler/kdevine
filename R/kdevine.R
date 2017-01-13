#' Kernel density estimatior based on simplified vine copulas
#'
#' Implements the vine-copula based estimator of Nagler and Czado (2016). The
#' marginal densities are estimated by \code{\link{kde1d}}, the vine copula
#' density by \code{\link{kdevinecop}}.
#'
#' @param data (\eqn{n x d}) data matrix.
#' @param mult.1d numeric; all bandwidhts for marginal kernel density estimation
#' are multiplied with \code{mult.1d}.
#' @param xmin numeric vector of length d; see \code{\link{kde1d}}.
#' @param xmax numeric vector of length d; see \code{\link{kde1d}}.
#' @param copula.type either \code{"kde"} (default) or \code{"parametric"} for
#' kernel or parametric estimation of the vine copula.
#' @param ... further arguments passed to \code{\link{kde1d}} or
#'  \code{\link{kdevinecop}}.
#'
#' @return An object of class \code{kdevine}.
#'
#' @seealso
#' \code{\link{dkdevine}}
#' \code{\link{kde1d}}
#' \code{\link{kdevinecop}}
#'
#' @references
#' Nagler, T., Czado, C. (2016) \cr
#' Evading the curse of dimensionality in nonparametric density estimation with
#' simplified vine copulas. \cr
#' \emph{Journal of Multivariate Analysis 151, 69-89 (doi:10.1016/j.jmva.2016.07.003)}
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
#' @export
kdevine <- function(data, mult.1d = 1, xmin = NULL, xmax = NULL, copula.type = "kde", ...) {
    data <- as.matrix(data)
    d <- ncol(data)

    ## sanity checks
    if (!is.null(xmin)) {
        if (length(xmin) != d)
            stop("'xmin' has to be of length d")
    }
    if (!is.null(xmax)) {
        if (length(xmax) != d)
            stop("'xmin' has to be of length d")
    }
    if (length(list(...)$bw) != d && !is.null(list(...)$bw))
        stop("'bw' hast to be of length d")
    if (is.null((list(...)$copula.type))) {
        copula.type <- "kde"
    } else {
        copula.type <- list(...)$copula.type
    }
    if (ncol(data) != d)
        data <- t(data)
    stopifnot(ncol(data) == d)

    ## estimation of the marginals
    marg.dens <- as.list(numeric(d))
    for (k in 1:d) {
        marg.dens[[k]] <- kde1d(data[, k],
                                xmin = list(...)$xmin[k],
                                xmax = list(...)$xmax[k],
                                bw   = list(...)$bw[k],
                                mult = mult.1d)
    }
    res <- list(marg.dens = marg.dens)

    ## estimation of the R-vine copula (only if d > 1)
    if (d > 1) {
        # transform to copula data
        u <- sapply(1:d, function(k) pkde1d(data[, k], marg.dens[[k]]))

        if (copula.type == "kde") {
            res$vine  <- suppressWarnings(
                kdevinecop(u,
                           matrix      = list(...)$matrix,
                           method      = list(...)$method,
                           mult        = list(...)$mult,
                           info        = list(...)$info,
                           test.level  = list(...)$test.level,
                           trunc.level = list(...)$trunc.level,
                           treecrit    = list(...)$treecrit,
                           cores       = list(...)$cores)
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
    x <- as.matrix(x)
    n <- length(obj$marg.dens[[1]]$data)
    if (ncol(x) == 1)
        x <- t(x)
    d <- ncol(x)

    stopifnot(class(obj) == "kdevine")
    if (length(obj$marg.dens) != d)
        stop("'x' has incorrect dimension")

    ## evaluate marginals
    margvals <- u <- x
    for(i in 1:d){
        margvals[, i] <- dkde1d(x[, i], obj$marg.dens[[i]])
    }

    ## evaluate copula density (if necessary)
    if (!is.null(obj$vine)) {
        # PIT to copula level
        for (i in 1:d)
            u[, i] <- pkde1d(x[, i], obj$marg.dens[[i]])
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
rkdevine <- function(n, obj) {
    # simulate from copula
    usim <- switch(obj$copula.type,
                   "kde" = rkdevinecop(n, obj$vine),
                   "parametric" = RVineSim(n, obj$vine))
    # use quantile transformation for marginals
    sapply(seq_len(ncol(usim)),
           function(i) qkde1d(usim[, i], obj$marg.dens[[i]]))
}




