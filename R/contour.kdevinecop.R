#' Contour plots of pair copula kernel estimates
#'
#' @param x a \code{\link{kdevinecop}} object.
#' @param tree \code{"ALL"} or integer vector; specifies which trees are
#' plotted.
#' @param xylim numeric vector of length 2; sets \code{xlim} and \code{ylim}
#' for the contours.
#' @param cex.nums numeric; expansion factor for font of the numbers.
#' @param ... arguments passed to \code{\link{contour.kdecopula}}.
#'
#' @examples
#' data(wdbc, package = "kdecopula")                     # load data
#' u <- VineCopula::pobs(wdbc[, 5:7], ties = "average")  # rank-transform
#' \dontshow{wdbc <- wdbc[1:30, ]}
#' # estimate density
#' fit <- kdevinecop(u)
#'
#' # contour matrix
#' contour(fit)
#'
#' @importFrom graphics abline contour par plot.new plot.window polygon
#' @importFrom graphics strheight strwidth text
#' @importFrom grDevices col2rgb rgb
#'
#' @export
contour.kdevinecop <- function(x, tree = "ALL", xylim = NULL, cex.nums = 1, ...) {

    ## check input
    M <- x$matrix
    d <- nrow(M)
    if (all(tree == "ALL"))
        tree <- seq.int(d-1)
    n.tree <- length(tree)
    if (!is.null(list(...)$type))
        stop("Only contour plots allowed. Don't use the type argument!")

    ## set up for plotting windows (restore settings on exit)
    usr <- par(mfrow = c(n.tree, d - min(tree)), mar = rep(0, 4))  # dimensions of contour matrix
    on.exit(par(usr))

    ## default style --------------------------------------------------
    # headings: create blue color scale
    TUMblue <- rgb(0, 103/255, 198/255)
    tint.seq <- seq(0.5, 0.5, length.out = d - 1)
    clrs <- rev(sapply(tint.seq, function(x) tint(TUMblue, x, 0.7)))

    # contours: set limits for plots
    if (!is.null(list(...)$margins)) {
        margins <- list(...)$margins
        if (!(margins %in% c("norm", "unif", "exp", "flexp")))
            stop("margins have to be one of 'norm', 'unif', 'exp', 'flexp'")
    } else {
        margins <- "norm"
    }
    if (is.null(xylim))
        xylim <- switch(margins,
                        "norm"  = c(-3, 3),
                        "unif"  = c(0, 1 - 1e-2),
                        "exp"   = c(0, 10),
                        "flexp" = c(-10, 0))
    xlim <- ylim <- xylim

    # contours: adjust limits for headings
    offs <- 0.25
    mult <- 1.35
    ylim[2] <- ylim[2] + offs*diff(ylim)


    ## run through trees -----------------------------------------------
    # initialize check variables
    cnt <- 0
    k <- d
    maxnums <- get_num(1, tree = max(tree), x$matrix)
    for (i in rev(tree)) {
        for (j in 1:(d - min(tree))) {
            if (j <= d - i) {
                name <- split_num(naming(M[c(j, (d - i + 1):d), j]))
                nums <- lapply(extract_nums(x[[i]]), split_num)
                match <- sapply(nums, function(x) all(name %in% x))
                cfit <- x[[i]][[which(match)]]$c
                if (!all(name[1:2] == nums[[which(match)]][1:2]))
                    cfit$flip <- TRUE


                # set up list of contour arguments
                args <- list(x = cfit,
                             drawlabels = FALSE,
                             xlab = "",
                             ylab = "",
                             xlim = xlim,
                             ylim = ylim,
                             xaxt = "n",
                             yaxt = "n",
                             add  = TRUE)

                # create empty plot
                plot.new()
                plot.window(xlim = xlim, ylim = ylim,
                            xaxs = "i",  yaxs = "i")

                # call contour.BiCop with ... arguments
                do.call(contour, args)#modifyList(args, list(...)))

                # draw area for headings
                abline(h = ylim[2] - diff(ylim)/mult*offs)
                ci <- min(length(clrs) + 1 - i, 10)
                polygon(x = c(xlim[1] - diff(xlim),
                              xlim[1] - diff(xlim),
                              xlim[2] + diff(xlim),
                              xlim[2] + diff(xlim)),
                        y = c(ylim[2] + diff(ylim)/mult*offs,
                              ylim[2] - diff(ylim)/mult*offs,
                              ylim[2] - diff(ylim)/mult*offs,
                              ylim[2] + diff(ylim)/mult*offs),
                        col = clrs[ci])

                # add separating lines
                abline(v = xlim)
                abline(h = ylim)

                # add pair-copula ID
                cx1 <- 0.75 * diff(xlim) / strwidth(maxnums)
                ty <- ylim[2] - diff(ylim)/mult*offs
                cx2 <- 0.75 * (ylim[2] - ty) / strheight(maxnums)
                cx <- min(cx1, cx2)
                text(x = sum(xlim)/2,
                     y = ty + 0.225 / cex.nums * (ylim[2] - ty),
                     cex    = cex.nums * cx,
                     labels = get_num(j, tree = i, M),
                     pos    = 3,
                     offset = 0)
            } else {
                plot.new()
            }
        }
    }
    invisible(x)
}

get_num <-  function(j, tree, M) {
    d <- nrow(M)
    # get numbers from structure matrix
    nums <- as.character(M[c(j, (d - tree + 1):d), j])
    # conditioned set
    bef <- paste(nums[2],
                 nums[1],
                 sep = ",",
                 collapse = "")
    # conditioning set
    aft <- if (length(nums) > 2) {
        gsub(" ",
             ",",
             do.call(paste, as.list(as.character(nums[3:length(nums)]))))
    }  else ""
    # paste together
    sep <- if (length(nums) > 2) " ; " else ""
    paste(bef, aft, sep = sep, collapse = "")
}

tint <- function(x, fac, alpha = 1) {
    x <- c(col2rgb(x))
    x <- (x + (255 - x) * fac) / 255
    rgb(x[1], x[2], x[3], alpha)
}

