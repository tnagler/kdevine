#' Plotting kde1d objects
#'
#' @aliases lines.kde1d
#' @method plot kde1d
#'
#' @param x \code{kde1d} object.
#' @param ev gridpoints for the plot (optional)
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


