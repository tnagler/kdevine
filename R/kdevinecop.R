#' Kernel estimation of vine copula densities
#'
#' The function estimates a vine copula density using kernel estimators for the
#' pair copulas (based on the \link{kdecopula} package).
#'
#' @param data (\eqn{n x d}) matrix of copula data (have to lie in \eqn{[0,1^d]}).
#' @param matrix R-Vine matrix (\eqn{n x d}) specifying the structure of the vine;
#'  if \code{NA} (default) the structure selection heuristic of Dissman et al.
#'  (2013) is applied.
#' @param method see \code{\link[kdecopula:kdecop]{kdecop}}.
#' @param renorm.iter see \code{\link[kdecopula:kdecop]{kdecop}}.
#' @param mult see \code{\link[kdecopula:kdecop]{kdecop}}.
#' @param test.level significance level for independence test. If you provide a
#' number in \eqn{[0, 1]}, an independence test
#' (\code{\link[VineCopula:BiCopIndTest]{BiCopIndTest}}) will be performed for
#' each pair; if the null hypothesis of independence cannot be rejected, the
#' independence copula will be set for this pair. If \code{test.level = NA}
#' (default), no independence test will be performed.
#' @param trunc.level integer; the truncation level. All pair copulas in trees
#' above the truncation level will be set to independence.
#' @param treecrit criterion for structure selection; defaults to \code{"tau"}.
#' @param cores integer; if \code{cores > 1}, estimation will be parallized
#' within each tree (using \code{\link[foreach]{foreach}}).
#' @param info logical; if \code{TRUE}, additional information about the
#' estimate will be gathered (see \code{\link[kdecopula:kdecop]{kdecop}}).
#'
#' @return An object of class \code{kdevinecop}. That is, a list containing
#' \item{T1, T2, ...}{lists of the estimted pair copulas in each tree,}
#' \item{matrix}{the structure matrix of the vine,}
#' \item{info}{additional information about the fit (if \code{info = TRUE}).}
#'
#' @references
#' Nagler, T., Czado, C. (2016) \cr Evading the curse of
#' dimensionality in nonparametric density estimation with simplified vine
#' copulas. \cr \emph{Journal of Multivariate Analysis 151, 69-89
#' (doi:10.1016/j.jmva.2016.07.003)}
#'
#' Nagler, T., Schellhase, C. and Czado, C. (2017) \cr Nonparametric
#' estimation of simplified vine copula models: comparison of methods
#' arXiv:1701.00845
#'
#' Dissmann, J., Brechmann, E. C., Czado, C., and Kurowicka, D. (2013). \cr
#' Selecting and estimating regular vine copulae and application to financial
#' returns. \cr
#' Computational Statistics & Data Analysis, 59(0):52--69.
#'
#' @seealso
#' \code{\link{dkdevinecop}},
#' \code{\link[kdecopula:kdecop]{kdecop}},
#' \code{\link[VineCopula:BiCopIndTest]{BiCopIndTest}},
#' \code{\link[foreach]{foreach}}
#'
#' @examples
#' data(wdbc, package = "kdecopula")
#' # rank-transform to copula data (margins are uniform)
#' u <- VineCopula::pobs(wdbc[, 5:7], ties = "average")
#' \dontshow{u <- u[1:30, ]}
#' fit <- kdevinecop(u)                   # estimate density
#' dkdevinecop(c(0.1, 0.1, 0.1), fit)     # evaluate density estimate
#' contour(fit)                           # contour matrix (Gaussian scale)
#' pairs(rkdevinecop(500, fit))           # plot simulated data
#'
#' @importFrom kdecopula kdecop hkdecop
#' @importFrom VineCopula BiCopIndTest RVineMatrix TauMatrix
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
kdevinecop <- function(data, matrix = NA, method = "TLL2", renorm.iter = 3L,
                       mult = 1, test.level = NA, trunc.level = NA,
                       treecrit = "tau", cores = 1, info = FALSE) {
    ## adjust input
    if (is.null(info))
        info <- FALSE
    if (is.null(matrix))
        matrix <- NA
    if (is.null(method))
        method <- "TLL2"
    if (is.null(mult))
        mult <- 1
    if (is.null(test.level))
        test.level <- 1
    if (is.na(test.level))
        test.level <- 1
    if (is.null(trunc.level))
        trunc.level <- ncol(data)
    if (is.na(trunc.level))
        trunc.level <- ncol(data)
    if (is.null(treecrit))
        treecrit <- "tau"
    if (is.na(treecrit))
        treecrit <- "tau"
    if (is.null(cores))
        cores <- 1

    data <- as.matrix(data)
    data <- pobs(data, ties.method = "first")
    matrix <- as.matrix(matrix)

    ## sanity checks
    d <- ncol(data)
    n <- nrow(data)
    if (any(data >= 1) || any(data < 0))
        stop("Data has be in the interval [0,1].")
    if (n < 2)
        stop("Number of observations has to be at least 2.")
    if (d < 2)
        stop("Dimension has to be at least 2.")
    if (!(treecrit %in% c("tau", "AIC", "cAIC")))
        stop("'treecrit' not available; please choose either 'tau', 'AIC' or 'cAIC'")

    ## call structure selection routine if no matrix given
    if (any(is.na(matrix)) & d > 2) {
        return(structure_select2(data,
                                 type = 0,
                                 method = method,
                                 mult = mult,
                                 info = info,
                                 struct.crit = treecrit,
                                 test.level = test.level,
                                 trunc.level = trunc.level,
                                 renorm.iter = renorm.iter,
                                 cores = cores))
    } else if (any(is.na(matrix)) & d == 2) {
        # set standard matrix for d = 2 (there is only on possible structure)
        matrix <- matrix(c(2, 1, 0, 1), 2, 2)
    } else if (nrow(matrix) == 1 && matrix == 1) {  # select C-Vine
        return(structure_select2(data,
                                 type = 1,
                                 method = method,
                                 mult = mult,
                                 info = info,
                                 struct.crit = treecrit,
                                 test.level = test.level,
                                 trunc.level = trunc.level,
                                 renorm.iter = renorm.iter,
                                 cores = cores))
    }

    ## check matrix if provided
    if (nrow(matrix) != ncol(matrix))
        stop("Structure matrix has to be quadratic.")
    if (ncol(data) != ncol(matrix))
        stop("Dimensions of data and matrix don't match.")
    if (max(matrix) > nrow(matrix))
        stop("Error in the structure matrix.")

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
            on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
        }
    }

    ##
    M <- ToLowerTri(matrix)
    Mold <- M
    o <- diag(M)
    M <- reorderRVineMatrix(M)
    data <- data[, o[length(o):1]]
    MaxMat <- createMaxMat(M)
    CondDistr <- neededCondDistr(M)

    ## initialize objects
    res <- as.list(numeric(d - 1))
    for (i in 1:(d - 1))
        res[[i]] <- as.list(numeric(d - i))
    llikv <- array(0, dim = c(d, d, n))
    llik <- matrix(0, d, d)
    effp <- matrix(0, d, d)
    AIC <- matrix(0, d, d)
    cAIC <- matrix(0, d, d)
    BIC <- matrix(0, d, d)
    nms <- matrix("", d - 1, d - 1)
    V <- list()
    V$direct <- array(NA, dim = c(d, d, n))
    V$indirect <- array(NA, dim = c(d, d, n))
    V$direct[d, , ] <- t(data[, d:1])

    ## For independence pair-copulas
    indepinfo <- list(effp = 0,
                      likvalues = rep(1, n),
                      loglik = 0,
                      effp = 0,
                      AIC = 0,
                      cAIC = 0,
                      BIC = 0)

    for (k in d:2) {
        doEst <- function(i) {
            if (k > i) {
                m <- MaxMat[k, i]
                zr1 <- V$direct[k, i, ]

                zr2 <- if (m == M[k, i]) {
                    V$direct[k, (d - m + 1), ]
                } else {
                    V$indirect[k, (d - m + 1), ]
                }
                samples <- cbind(zr2, zr1)

                indep <- ifelse(test.level < 1,
                                BiCopIndTest(zr2, zr1)$p.value >= test.level,
                                FALSE)
                if (trunc.level <= (d-k))
                    indep <- TRUE

                if (indep) {
                    cfit <- list()
                    if (info)
                        cfit$info <- indepinfo
                    class(cfit) <- c("kdecopula", "indep.copula")
                } else {
                    cfit <- kdecop(samples,
                                   mult = mult,
                                   method = method,
                                   renorm.iter = renorm.iter,
                                   info = info)
                }

                hfit <- list()
                direct <- indirect <- NULL
                if (CondDistr$direct[k - 1, i]) {
                    direct <- hkdecop(samples,
                                      obj = cfit,
                                      cond.var = 1L)
                }
                if (CondDistr$indirect[k - 1, i]) {
                    indirect <- hkdecop(samples,
                                        obj = cfit,
                                        cond.var = 2L)
                }
                names <- naming(Mold[c(i, k:d), i])
                res.ki <- list(c = cfit, name = names)
                return(list(direct = direct,
                            indirect = indirect,
                            res.ki = res.ki))
            } else {
                return(NULL)
            }
        }
        res.k <- if (cores > 1) {
            foreach(i = 1:(k-1),
                    .export = c("naming"),
                    .packages = c("kdevine", "kdecopula")) %dopar% doEst(i)
        } else lapply(1:(k-1), doEst)

        ## clean up and finalize
        for (i in 1:(d-1)) {
            nums <- Mold[c(i, k:d), i]
            #nums[1:2] <- nums[2:1]
            name <- naming(nums)
            if (any(extract_nums(res.k) == name)) {
                ind <- which(extract_nums(res.k) == name)
                res.ki <- res.k[[ind]]
                res[[d + 1 - k]][[i]] <- res.ki$res.ki
                if (!is.null(res.ki$direct)) {
                    V$direct[k - 1, i, ] <- res.ki$direct
                }
                if (!is.null(res.ki$indirect)) {
                    V$indirect[k - 1, i, ] <- res.ki$indirect
                }

                if (info) {
                    cfit <- res.ki$res.ki$c
                    llikv[k, i, ] <- log(cfit$info$likvalues)
                    llik[k, i] <- cfit$info$loglik
                    effp[k, i] <- cfit$info$effp
                    AIC[k, i] <- -2 * cfit$info$loglik + 2 * effp[k, i]
                    cAIC[k, i] <-
                        AIC[k, i] + (2 * effp[k, i] * (effp[k, i] + 1)) /
                        (n - effp[k, i] - 1)
                    BIC[k, i] <- -2 * cfit$info$loglik + log(n) * effp[k, i]
                }
            }
        } # end i = 1:(d-1)
    } # end k = d:2

    ## finalize results
    res[[d]] <- data
    res[[d + 1]] <- Mold
    res[[d + 2]] <- list(llikv       = apply(llikv, 3, sum),
                         loglik      = sum(llik),
                         pair.loglik = llik,
                         effp        = sum(effp),
                         pair.effp   = effp,
                         AIC         = sum(AIC),
                         pair.AIC    = AIC,
                         cAIC        = sum(AIC) +
                             (2 * sum(effp) * (sum(effp) + 1)) /
                             (n - sum(effp) - 1),
                         pair.cAIC   = cAIC,
                         BIC         = sum(BIC),
                         pair.BIC    = BIC)
    names(res) <- vapply(1:(d - 1), function(x) paste("T", x, sep = ""), "")
    names(res)[d:(d + 2)] <- c("data", "matrix", "info")

    if (!info)
        res[[d + 2]] <- NULL
    ## return results
    class(res) <- "kdevinecop"
    res
}




## required functions from VineCopula package ----------------------------------

normalizeRVineMatrix <- function(RVM) {

    oldOrder <- diag(RVM$Matrix)
    Matrix <- reorderRVineMatrix(RVM$Matrix)

    return(RVineMatrix(Matrix,
                       RVM$family,
                       RVM$par,
                       RVM$par2,
                       names = rev(RVM$names[oldOrder])))
}

reorderRVineMatrix <- function(Matrix) {
    oldOrder <- diag(Matrix)

    O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)

    for (i in 1:nrow(Matrix)) {
        Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
    }

    return(Matrix)
}


createMaxMat <- function(Matrix) {
    if (dim(Matrix)[1] != dim(Matrix)[2])
        stop("Structure matrix has to be quadratic.")

    MaxMat <- reorderRVineMatrix(Matrix)
    n <- nrow(MaxMat)

    for (j in 1:(n - 1)) {
        for (i in (n - 1):j) {
            MaxMat[i, j] <- max(MaxMat[i:(i + 1), j])
        }
    }

    tMaxMat <- MaxMat
    tMaxMat[is.na(tMaxMat)] <- 0
    oldSort <- diag(Matrix)
    oldSort <- oldSort[n:1]

    for (i in 1:n) {
        MaxMat[tMaxMat == i] <- oldSort[i]
    }

    return(MaxMat)
}

neededCondDistr <- function(Vine) {
    if (dim(Vine)[1] != dim(Vine)[2])
        stop("Structure matrix has to be quadratic.")

    Vine <- reorderRVineMatrix(Vine)
    MaxMat <- createMaxMat(Vine)
    d <- nrow(Vine)

    M <- list()
    M$direct <- matrix(FALSE, d, d)
    M$indirect <- matrix(FALSE, d, d)
    M$direct[2:d, 1] <- TRUE

    for (i in 2:(d - 1)) {
        v <- d - i + 1
        bw <- as.matrix(MaxMat[i:d, 1:(i - 1)]) == v
        direct <- Vine[i:d, 1:(i - 1)] == v

        M$indirect[i:d, i] <- apply(as.matrix(bw & (!direct)), 1, any)
        M$direct[i:d, i] <- TRUE
        M$direct[i, i] <- any(as.matrix(bw)[1, ] & as.matrix(direct)[1, ])
    }

    return(M)
}

ToLowerTri <- function(x) {
    ## only change matrix if not already lower triagonal
    if(all(x[lower.tri(x)] == 0)) {
        x[nrow(x):1, ncol(x):1]
    } else {
        x
    }
}
