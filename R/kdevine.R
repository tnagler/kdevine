#' High-dimensional kernel density estimation based on vine copulas
#' 
#' @param data (\eqn{n x d}) data matrix.
#' @param mult.1d numeric; all bandwidhts for univariate kernel density estimation
#'  are multiplied with \code{mult.1d}.
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
#' @examples
#' data(wdbc)  # load data
#' fit <- kdevine(wdbc[, 5:7])  # estimate density
#' dkdevine(c(1000, 0.1, 0.1), fit)  # evaluate density estimate
#' 
#' @importFrom VineCopula RVineStructureSelect RVineCopSelect
#' @export
kdevine <- function(data, mult.1d = 1, ...) {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
    
    ## sanity checks
    if (!is.null(list(...)$xmin)) {
        if(length(list(...)$xmin) != d)
            stop("'xmin' has to be of length d")
    }
    if (!is.null(list(...)$xmax)) {
        if(length(list(...)$xmax) != d)
            stop("'xmin' has to be of length d")
    }
    if (length(list(...)$bw) != d && !is.null(list(...)$bw))
        stop("'bw' hast o be of length d")
    if (is.null((list(...)$copula.type)))
        copula.type <- "kernel"
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
                                mult = mult.1d,
                                fast = list(...)$fast)
    }
    res <- list(marg.dens = marg.dens)
    
    ## estimation of the R-vine copula (only if d > 1)
    if (d > 1) {
        # transform to copula data
        u <- sapply(1:d, function(k) pkde1d(data[, k], marg.dens[[k]]))
        
        # estimate vine for type "kernel"
        if (copula.type == "kde") {
            res$vine  <- suppressWarnings(
                kdevinecop(u,
                           matrix      = list(...)$matrix,
                           method      = list(...)$method,
                           info        = list(...)$info,
                           test.level  = list(...)$test.level,
                           trunc.level = list(...)$trunc.level,
                           cores       = list(...)$cores)
            )
        }
        # estimate vine for type "parametric"
        if (copula.type == "parametric") {
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
#' data(wdbc)  # load data
#' fit <- kdevine(wdbc[, 5:7])  # estimate density
#' dkdevine(c(1000, 0.1, 0.1), fit)  # evaluate density estimate
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
        # PIT with boundary correction
        for (i in 1:d)
            u[, i] <- pkde1d(x[, i], obj$marg.dens[[i]]) * n / (n + 1)
        vinevals <- dkdevinecop(u, obj = obj$vine, stable = TRUE) 
    } else {
        vinevals <- rep(1, nrow(x))
    }
    
    ## final density estimate is product of marginals and copula density
    apply(cbind(margvals, vinevals), 1, prod)
}

