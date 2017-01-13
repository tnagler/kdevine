#' Working with a \code{kdevinecop} object
#'
#' A vine copula density estimate (stored in a \code{kdevinecop} object)
#' can be evaluated on arbitrary points with \code{dkevinecop}. Furthermore,
#' you can simulate from the estimated density with \code{rkdevinecop}.
#'
#' @aliases dkevinecop rkdevinecop
#'
#' @param u \eqn{m x 2} matrix of evaluation points.
#' @param obj \code{kdevinecop} object.
#' @param stable logical; option for stabilizing the estimator: the estimated
#' pair copula density is cut off at \eqn{50}.
#'
#' @return A numeric vector of the density/cdf or a \eqn{n x 2} matrix of
#' simulated data.
#'
#' @author Thomas Nagler
#'
#' @seealso
#' \code{\link{kdevinecop}},
#' \code{\link[kdecopula:dkdecop]{dkdecop}},
#' \code{\link[kdecopula:rkdecop]{rkdecop}},
#' \code{\link[qrng:ghalton]{ghalton}}
#'
#' @references
#' Nagler, T., Czado, C. (2016) \cr
#' Evading the curse of dimensionality in nonparametric density estimation. \cr
#' Journal of Multivariate Analysis 151, 69-89 (doi:10.1016/j.jmva.2016.07.003)
#'
#' Dissmann, J., Brechmann, E. C., Czado, C., and Kurowicka, D. (2013). \cr
#' Selecting and estimating regular vine copulae and application to financial returns. \cr
#' Computational Statistics & Data Analysis, 59(0):52--69.
#'
#' @examples
#' data(wdbc, package = "kdecopula")                    # load data
#' u <- VineCopula::pobs(wdbc[, 5:7], ties = "average") # rank-transform
#'
#' fit <- kdevinecop(u)                # estimate density
#' dkdevinecop(c(0.1, 0.1, 0.1), fit)  # evaluate density estimate
#'
#' @importFrom kdecopula dkdecop hkdecop
#' @export
dkdevinecop <- function(u, obj, stable = FALSE) {
    stopifnot(is.numeric(u))
    stopifnot(ncol(u) == ncol(obj$matrix))
    if (any(u > 1) || any(u < 0))
        stop("Values have to be in the interval (0,1).")

    u <- as.matrix(u)
    if (ncol(u) == 1)
        u <- matrix(u, 1, nrow(u))
    n <- nrow(u)
    d <- ncol(u)
    if (ncol(u) != ncol(obj$matrix))
        stop("Dimensions of 'u' and 'matrix' do not match.")

    ##
    M <- as.matrix(obj$matrix)
    Mold <- M
    o <- diag(M)
    M <- reorderRVineMatrix(M)
    uev <- as.matrix(u[, o[length(o):1]])
    if (ncol(uev) == 1)
        uev <- matrix(uev, 1, nrow(uev))

    MaxMat <- createMaxMat(M)
    CondDistr <- neededCondDistr(M)

    ## initialize objects
    res <- as.list(numeric(d - 1))
    for (i in 1:(d - 1))
        res[[i]] <- as.list(numeric(d - i))
    V <- list()
    V$direct <- array(NA, dim = c(d, d, n))
    V$indirect <- array(NA, dim = c(d, d, n))
    V$direct[d, , ] <- t(uev[, d:1])

    val <- array(1, dim = c(d, d, n))

    for (k in d:2) {
        for (i in 1:(k - 1)) {
            # get names and match object
            name <- split_num(naming(Mold[c(i, k:d), i]))
            nums <- lapply(extract_nums(obj[[d - k + 1]]), split_num)
            match <- sapply(nums, function(x) all(name %in% x))
            if (any(match)) {
                ## get object and pseudo-observations
                res.ki <- obj[[d - k + 1]][[which(match)]]
                m <- MaxMat[k, i]
                zr1 <- V$direct[k, i, ]
                zr2 <- if (m == M[k, i]) {
                    V$direct[k, (d - m + 1), ]
                } else {
                    V$indirect[k, (d - m + 1), ]
                }
                ev <- cbind(zr2, zr1)
                cfit <- res.ki$c

                ## flip grid if arguments in object are in different order
                if (!all(name[1:2] == nums[[which(match)]][1:2]))
                    cfit$flip <- TRUE

                ## evaluate density and h-functions
                if (any(class(cfit) == "indep.copula")) {
                    val[k-1, i, ] <- rep(1, nrow(ev))
                } else {
                    val[k-1, i, ] <- dkdecop(ev,
                                             obj = cfit,
                                             stable = stable)
                }

                if (k > 2) {
                    if(CondDistr$direct[k - 1, i]) {
                        if (any(class(cfit) == "indep.copula")) {
                            V$direct[k - 1, i, ] <- ev[,2]
                        } else {
                            V$direct[k - 1, i, ] <- hkdecop(ev,
                                                            obj = cfit,
                                                            cond.var = 1L)
                        }
                    }

                    if(CondDistr$indirect[k - 1, i]) {
                        if (any(class(cfit) == "indep.copula")) {
                            V$indirect[k - 1, i, ] <- ev[,1]
                        } else {
                            V$indirect[k - 1, i, ] <- hkdecop(ev,
                                                              obj = cfit,
                                                              cond.var = 2L)
                        }
                    }

                }
            }
        }
    }
    out <- exp(apply(log(val), 3, sum))
    if (stable)
        out <- pmin(out, 10^(1 + d/2))
    out
}

#' @param n integer; number of observations.
#' @param U (optional) \eqn{n x d} matrix of independent uniform random
#'  variables.
#' @param quasi logical; the default (\code{FALSE}) returns pseudo-random
#' numbers, use \code{TRUE} for quasi-random numbers (generalized Halton, see
#' \code{\link[qrng:ghalton]{ghalton}}).
#'
#' @rdname dkdevinecop
#' @importFrom kdecopula hkdecop
#' @importFrom stats runif
#' @export
rkdevinecop <- function(n, obj, U = NULL, quasi = FALSE) {
    n <- round(n)
    stopifnot(n > 0)
    stopifnot(is.logical(quasi))

    ## get structurte matrix and helper matrices,
    ## convert all to uper triangular form for simpler indexing
    M <- obj$matrix
    o <- diag(M)
    d <- length(o)
    Mold <- M[d:1, d:1]
    M <- reorderRVineMatrix(M)
    MaxMat <- createMaxMat(M)[d:1, d:1]
    CondDistr <- lapply(neededCondDistr(M), function(m) m[d:1, d:1])
    M <- M[d:1, d:1]

    ## initialize V matrices with independent uniform data
    ## (see Scherer, Mai (2014). "Simulating Copulas", Chapter 4)
    stopifnot(is.logical(quasi))
    if (!quasi) {
        # simulate independent uniform random variables
        if (is.null(U)) {
            U <- matrix(runif(n * d), ncol = d)
        } else {
            U <- matrix(U, ncol = d)
            U <- U[, rev(o), drop = FALSE]
        }
    } else {
        # generate quasi random numbers
        U <- ghalton(n, d = d)
    }

    Vdirect <- Vindirect <- array(dim = c(d, d, n))
    for (i in 1:d) {
        Vdirect[i, i, ] <- U[, i]
    }
    Vindirect[1, 1, ] <- Vdirect[1, 1, ]

    ## simulation algorithm
    for (i in 2:d) {  # loop through variables
        for (k in (i - 1):1) {  # loop over pairs involving variable i
            # find the pair copula object
            name <- split_num(naming(Mold[c(i, k:1), i]))
            nums <- lapply(extract_nums(obj[[k]]), split_num)
            match <- sapply(nums, function(x) all(name %in% x))
            res.ki <- obj[[k]][[which(match)]]

            # find the corresponding pseudo-data
            mm <- MaxMat[k, i]
            if (mm == M[k, i]) {
                zz <- Vdirect[k, mm, ]
            } else {
                zz <- Vindirect[k, mm, ]
            }

            # flip grid if arguments in object are in different order
            if (!all(name[1:2] == nums[[which(match)]][1:2]))
                res.ki$c$flip <- TRUE

            # invert h-function to create pseudo-data for next step
            Vdirect[k, i, ] <- hkdecop(cbind(zz, Vdirect[k + 1, i, ]),
                                       res.ki$c,
                                       cond.var = 1,
                                       inverse = TRUE)
            if ((i < d) & CondDistr$indirect[k + 1, i]) {
                Vindirect[k + 1, i, ] <- hkdecop(cbind(zz, Vdirect[k, i, ]),
                                                 res.ki$c,
                                                 cond.var = 2)
            }
        }
    }

    ## return result in initial ordering of the variables
    out <- t(Vdirect[1, , ])
    out[, sort(o[length(o):1], index.return = TRUE)$ix, drop = FALSE]
}
