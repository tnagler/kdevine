#' Univariate kernel density estimation for bounded and unbounded support
#'
#' Function is deprecated, use [`kde1d::kde1d()`] instead.
#'
#' @param x vdeprecated.
#' @param mult deprecated.
#' @param xmin deprecated.
#' @param xmax deprecated..
#' @param bw deprecated.
#' @param bw_min deprecated.
#' @param ... deprecated.
#' @export
kde1d <- function(x, mult = 1, xmin = -Inf, xmax = Inf, bw = NULL, bw_min = 0, ...) {
    stop("kdevine::kde1d() is deprecated; use the kde1d package instead.")
}


#' @rdname kde1d
#' @param obj deprecated.
#' @export
pkde1d <- function(x, obj) {
    stop("kdevine::kde1d() is deprecated; use the kde1d package instead.")
}

#' @rdname kde1d
#' @export
qkde1d <- function(x, obj) {
    stop("kdevine::kde1d() is deprecated; use the kde1d package instead.")
}

#' @param n deprecated.
#' @param quasi deprecated.
#'
#' @rdname kde1d
#' @export
rkde1d <- function(n, obj, quasi = FALSE) {
    stop("kdevine::kde1d() is deprecated; use the kde1d package instead.")
}
