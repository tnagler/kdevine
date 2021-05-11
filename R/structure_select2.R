structure_select2 <- function(data, type, method, mult, struct.crit, test.level,
                              trunc.level, renorm.iter, cores, info, progress = FALSE) {
    if (type == 0) {
        type <- "RVine"
    } else if (type == 1) {
        type <- "CVine"
    }
    if (type != "RVine" & type != "CVine")
        stop("Vine model not implemented.")
    n <- nrow(data)
    d <- ncol(data)
    if (n < 2)
        stop("Number of observations has to be at least 2.")
    if (d < 2)
        stop("Dimension has to be at least 2.")
    if (any(data > 1) || any(data < 0))
        stop("Data has be in the interval [0,1].")
    if (is.na(test.level))
        test.level <- 1
    if (!(struct.crit %in% c("tau", "AIC", "cAIC", "hoeffd")))
        stop("'struct.crit' has to be one of 'tau', 'AIC', 'cAIC', or 'hoeffd'")

    # set names if not available
    if (is.null(colnames(data)))
        colnames(data) = paste("V", 1:d, sep = "")

    ## initialize objects
    RVine <- list(Tree = NULL, Graph = NULL)
    res <- as.list(numeric(d - 1))
    for (i in 1:(d - 1))
        res[[i]] <- as.list(numeric(d - i))
    llikv <- array(0, dim = c(d, d, n))
    llik <- matrix(0, d, d)
    effp <- matrix(0, d, d)
    AIC <- matrix(0, d, d)
    cAIC <- matrix(0, d, d)
    BIC <- matrix(0, d, d)

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
        }
    }

    ## build first tree
    g <- initializeFirstGraph2(data, struct.crit = struct.crit, weights = NA)
    mst <- findMaxTree2(g, mode = type)

    ## estimate copulas in first tree and store results
    est <- est.FirstTreeCopulas2(mst,
                                 data,
                                 method = method,
                                 mult = mult,
                                 info = info,
                                 test.level = test.level,
                                 renorm.iter = renorm.iter,
                                 parallel = cores > 1)
    res[[1]] <- est$est
    VineTree <- est$mst
    RVine$Tree[[1]] <- VineTree
    RVine$Graph[[1]] <- g

    ## higher trees
    for (k in 2:(d - 1)) {
        g <- buildNextGraph2(VineTree,
                             weights = NA,
                             struct.crit = struct.crit,
                             parallel = cores > 1)
        mst <- findMaxTree2(g, mode = type)
        est <- est.TreeCopulas2(mst,
                                k = k,
                                d2 = d,
                                data = data,
                                oldVineGraph = VineTree,
                                method = method,
                                mult = mult,
                                info = info,
                                test.level = test.level,
                                renorm.iter = renorm.iter,
                                parallel = cores > 1,
                                truncate = trunc.level <= k)

        res[[k]] <- est$est
        VineTree <- est$tree
        RVine$Tree[[k]] <- VineTree
        RVine$Graph[[k]] <- g
    }

    ## finalize results
    res[[d]] <- data
    res[[d + 1]] <- M <- as.RVMKernel2(RVine)$Matrix
    for (i in 1:(d - 1)) {
        for (k in (i + 1):d) {
            hfit <- res[[i]][[k - i]]$h
            if (info) {
                pcfit <- res[[i]][[k - i]]$c
                llikv[k, i, ] <- log(pcfit$info$likvalues)
                llik[k, i] <- pcfit$info$loglik
                effp[k, i] <- pcfit$info$effp
                AIC[k, i] <- -2 * pcfit$info$loglik + 2 * effp[k, i]
                cAIC[k, i] <-
                    AIC[k, i] + (2 * effp[k, i] * (effp[k, i] + 1)) /
                    (n - effp[k, i] - 1)
                BIC[k, i] <- -2 * pcfit$info$loglik + log(n) * effp[k, i]
            }
        }
    }
    res[[d + 2]] <- if (info) {
        list(llikv       = apply(llikv, 3, sum),
             loglik      = sum(llik),
             pair.loglik = llik,
             effp        = sum(effp),
             pair.effp   = effp,
             AIC         = sum(AIC),
             pair.AIC    = AIC,
             cAIC        =
                 sum(AIC) +
                 (2 * sum(effp) * (sum(effp) + 1)) /
                 (n - sum(effp) - 1),
             pair.cAIC   = cAIC,
             BIC         = sum(BIC),
             pair.BIC    = BIC)
    } else list(NULL)

    names(res) <- vapply(1:(d - 1), function(x) paste("T", x, sep = ""), "")
    names(res)[d:(d + 2)] <- c("data", "matrix", "info")

    ## return results
    class(res) <- "kdevinecop"
    res
}

est.FirstTreeCopulas2 <- function(mst, data.univ, method, mult, test.level,
                                  renorm.iter, info, parallel) {
    d <- nrow(mst$E$nums)

    ## estimation
    pkgs <- c("kdevine", "kdecopula")

    ## For independence pair-copulas
    indepinfo <- list(effp = 0,
                      likvalues = rep(1, nrow(data.univ)),
                      loglik = 0,
                      effp = 0,
                      AIC = 0,
                      cAIC = 0,
                      BIC = 0)

    doEst <- function(i) {
        a <- rev(mst$E$nums[i, ])

        ## store samples
        Copula.Data.1 =  list(data.univ[, a[1]])
        Copula.Data.2 =  list(data.univ[, a[2]])

        ## set names for this edge
        if (is.null(mst$V$names[a[1]])) {
            Copula.CondName.1 <- a[1]
        } else {
            Copula.CondName.1 <- mst$V$names[a[1]]
        }
        if (is.null(mst$V$names[a[2]])) {
            Copula.CondName.2 <- a[2]
        } else {
            Copula.CondName.2 <- mst$V$names[a[2]]
        }
        if (is.null(mst$V$names[a[1]]) || is.null(mst$V$names[a[2]])) {
            Copula.Name <- paste(a[1],
                                 a[2],
                                 sep = " , ")
        } else {
            Copula.Name <- paste(mst$V$names[a[1]],
                                 mst$V$names[a[2]],
                                 sep = " , ")
        }

        ## identify with numbers
        nums <- paste(a[2], a[1], sep = ",")

        ## estimation
        s <- cbind(data.univ[,a[2]], data.univ[,a[1]])
        indep <- ifelse(test.level < 1,
                        BiCopIndTest(s[, 1], s[, 2])$p.value >= test.level,
                        FALSE)
        if (is.na(indep))
            browser()

        if (indep) {
            pcfit <- list()
            if (info)
                pcfit$info <- indepinfo
            class(pcfit) <- c("kdecopula", "indep.copula")
        } else {
            pcfit <- kdecop(s,
                            method = method,
                            mult = mult,
                            renorm.iter = renorm.iter,
                            info = info)
            # pcfit$udata <- NULL
        }

        ## get and store new pseudo-samples
        if (indep == TRUE) {
            Copula.CondData.1 <- s[,1]
            Copula.CondData.2 <- s[,2]
        } else {
            Copula.CondData.1 <- list(hkdecop(s,
                                              obj = pcfit,
                                              cond.var = 2L))
            Copula.CondData.2 <- list(hkdecop(s,
                                              obj = pcfit,
                                              cond.var = 1L))
        }

        ## store estimates
        resi <- list(c = pcfit, name = nums)

        ## return results
        list(Copula.Data.1 = Copula.Data.1,
             Copula.Data.2 = Copula.Data.2,
             Copula.CondName.1 = Copula.CondName.1,
             Copula.CondName.2 = Copula.CondName.2,
             Copula.Name = Copula.Name,
             Copula.CondData.1 = Copula.CondData.1,
             Copula.CondData.2 = Copula.CondData.2,
             resi = resi)
    }

    res <- if (parallel) {
        foreach(i = 1:d,
                .export = c("d"),
                .packages = pkgs) %dopar% doEst(i)
    } else {
        lapply(1:d, doEst)
    }

    ## clean and return results
    for (i in 1:d) {
        mst$E$Copula.Data.1[[i]] =  res[[i]]$Copula.Data.1
        mst$E$Copula.Data.2[[i]] =  res[[i]]$Copula.Data.2
        mst$E$Copula.CondName.1[[i]] = res[[i]]$Copula.CondName.1
        mst$E$Copula.CondName.2[[i]] = res[[i]]$Copula.CondName.2
        mst$E$Copula.Name[[i]] = res[[i]]$Copula.Name
        mst$E$Copula.CondData.1[[i]] = res[[i]]$Copula.CondData.1
        mst$E$Copula.CondData.2[[i]] = res[[i]]$Copula.CondData.2
        mst$E$Copula.Nums.1 = res[[i]]$resi$name
        res[[i]] <- res[[i]]$resi
    }

    list(mst = mst, est = res)
}



est.TreeCopulas2 <- function(mst, k, d2, data, oldVineGraph, method, mult, info,
                             test.level, renorm.iter, weights = NA, parallel,
                             truncate) {

    d <- nrow(mst$E$nums) # number of nodes in the tree

    ## For independence pair-copulas
    indepinfo <- list(effp = 0,
                      likvalues = rep(1, nrow(data)),
                      loglik = 0,
                      effp = 0,
                      AIC = 0,
                      cAIC = 0,
                      BIC = 0)

    ## estimation
    exp <- c("split_name", "split_num", "naming")
    pkgs <- c("kdecopula")
    doEst <- function(i, mst) {
        ## get edge and corresponding data
        con <- rev(mst$E$nums[i, ])
        temp <- oldVineGraph$E$nums[con, ]

        ## find conditional variable
        if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
            same <- temp[2, 1]
        } else {
            if((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
                same <- temp[2, 2]
            }
        }

        ## find conditioned variables
        other1 <- temp[1, temp[1, ] != same]
        other2 <- temp[2, temp[2, ] != same]

        ## extract observations and names
        if (temp[1, 1] == same) {
            zr1 <- oldVineGraph$E$Copula.CondData.2[[con[1]]]
            n1  <- oldVineGraph$E$Copula.CondName.2[[con[1]]]
        } else {
            zr1 <- oldVineGraph$E$Copula.CondData.1[[con[1]]]
            n1  <- oldVineGraph$E$Copula.CondName.1[[con[1]]]
        }
        if (temp[2, 1] == same) {
            zr2 <- oldVineGraph$E$Copula.CondData.2[[con[2]]]
            n2  <- oldVineGraph$E$Copula.CondName.2[[con[2]]]
        } else {
            zr2 <- oldVineGraph$E$Copula.CondData.1[[con[2]]]
            n2  <- oldVineGraph$E$Copula.CondName.1[[con[2]]]
        }

        zr1a <- if (is.list(zr1)) as.vector(zr1[[1]]) else zr1
        zr2a <- if (is.list(zr2)) as.vector(zr2[[1]]) else zr2
        n1a <- if (is.list(n1)) as.vector(n1[[1]]) else n1
        n2a <- if (is.list(n2)) as.vector(n2[[1]]) else n2

        ## store pseudo-samples
        Copula.Data.1 <- zr1a
        Copula.Data.2 <- zr2a

        ## store names
        Copula.CondName.2 <- n1a
        Copula.CondName.1 <- n2a

        ## estimations
        samples <- cbind(zr1a, zr2a)
        indep <- ifelse(test.level < 1,
                        BiCopIndTest(samples[, 1], samples[, 2])$p.value >= test.level,
                        FALSE)

        if (truncate)
            indep <- TRUE
        if (indep) {
            pcfit <- list()
            if (info)
                pcfit$info <- indepinfo
            class(pcfit) <- c("kdecopula", "indep.copula")
        } else {
            pcfit <- kdecop(samples,
                            mult = mult,
                            method = method,
                            renorm.iter = renorm.iter,
                            info = info)
            # pcfit$udata <- NULL
        }

        ## naming
        tmpname <- mst$E$names[i]
        namesplt <- split_name(tmpname)
        numsplt  <- sapply(namesplt, function(x) which(colnames(data) == x))
        nums <- naming(sprintf("%d", numsplt))

        ## get and store new pseudo-samples
        if (indep == TRUE) {
            Copula.CondData.1 <- samples[,2]
            Copula.CondData.2 <- samples[,1]
        } else {
            Copula.CondData.1 <- hkdecop(samples,
                                         obj = pcfit,
                                         cond.var = 1L)
            Copula.CondData.2 <- hkdecop(samples,
                                         obj = pcfit,
                                         cond.var = 2L)
        }


        ## store estimates
        resi <- list(c = pcfit, name = nums)

        ## return results
        list(Copula.Data.1 = Copula.Data.1,
             Copula.Data.2 = Copula.Data.2,
             Copula.CondName.1 = Copula.CondName.1,
             Copula.CondName.2 = Copula.CondName.2,
             Copula.CondData.1 = Copula.CondData.1,
             Copula.CondData.2 = Copula.CondData.2,
             resi = resi)
    }

    res.k <- if (parallel) {
        foreach(i = 1:d, .export = exp, .packages = pkgs) %dopar% doEst(i, mst)
    } else {
        lapply(1:d, doEst, mst = mst)
    }

    ## clean and return results
    for (i in 1:d) {
        mst$E$Copula.Data.1[[i]] =  res.k[[i]]$Copula.Data.1
        mst$E$Copula.Data.2[[i]] =  res.k[[i]]$Copula.Data.2
        mst$E$Copula.CondName.1[[i]] = res.k[[i]]$Copula.CondName.1
        mst$E$Copula.CondName.2[[i]] = res.k[[i]]$Copula.CondName.2
        mst$E$Copula.CondData.1[[i]] = res.k[[i]]$Copula.CondData.1
        mst$E$Copula.CondData.2[[i]] = res.k[[i]]$Copula.CondData.2
        res.k[[i]] <- res.k[[i]]$resi
    }
    list(tree = mst, est = res.k)
}

hoeffd <- function(x) {
    n <- nrow(x)
    R <- apply(x,2,rank)-1
    Q <- sapply(1:n, function(i) sum(x[,1] < x[i,1] & x[,2] < x[i,2]))
    A <- sum(apply(R*(R-1),1,prod))
    B <- sum(apply((R-1),1,prod)*Q)
    C <- sum(Q*(Q-1))
    return((A - 2*(n-2)*B + (n-2)*(n-3)*C)/(n*(n-1)*(n-2)*(n-3)*(n-4)))
}

initializeFirstGraph2 <- function(data.univ, weights, struct.crit = "tau") {

    # C = cor(data.univ,method='kendall')
    q <- dim(data.univ)[2]
    C <- matrix(rep(1, q * q), ncol = q)

    for (i in 1:(q - 1)) {
        for (j in (i + 1):q) {
            if (struct.crit == "tau") {
                crit <- fasttau(data.univ[, i],
                                data.univ[, j],
                                weights)
            } else if (struct.crit == "AIC") {
                crit <- kdecop(data.univ[, c(i, j)], info = TRUE)$info$AIC
            } else if (struct.crit == "cAIC") {
                crit <- kdecop(data.univ[, c(i, j)], info = TRUE)$info$cAIC
            } else if (struct.crit == "hoeffd") {
                crit <- abs(hoeffd(data.univ[, c(i, j)]))
            }
            C[i, j] <- crit
            C[j, i] <- crit
        }
    }
    rownames(C) <- colnames(C) <- colnames(data.univ)
    rownames(C) <- colnames(C) <- colnames(data.univ)

    graphFromTauMatrix(C)
}


## functions for handling the tree structure -------------------------
graphFromTauMatrix2 <- function(tau) {
    d <- ncol(tau)
    # get variable names
    nms <- colnames(tau)
    # construct edge set
    E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
               do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
    # add edge names
    E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
    # set weights
    w <- tau[upper.tri(tau)]

    ## return results
    list(V = list(names = nms,
                  conditionedSet = NULL,
                  conditioningSet = NULL),
         E = list(nums = E,
                  names = E.names,
                  weights = w,
                  conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                  conditioningSet = NULL))
}


## initialize graph for next vine tree (possible edges)
buildNextGraph2 <- function(oldVineGraph, weights = NA, struct.crit = "tau", parallel) {
    d <- nrow(oldVineGraph$E$nums)

    ## initialize with full graph
    g <- makeFullGraph2(d)
    g$V$names <- oldVineGraph$E$names
    g$V$conditionedSet <- oldVineGraph$E$conditionedSet
    g$V$conditioningSet <- oldVineGraph$E$conditioningSet
    ## get info for all edges
    if (parallel) {
        i <- NULL  # dummy for CRAN check
        out <- foreach(i = 1:nrow(g$E$nums)) %dopar% getEdgeInfo2(i,
                                                                  g = g,
                                                                  oldVineGraph = oldVineGraph,
                                                                  weights = weights,
                                                                  struct.crit = struct.crit)
    } else {
        out <- lapply(1:nrow(g$E$nums),
                      getEdgeInfo2,
                      g = g,
                      oldVineGraph = oldVineGraph,
                      weights = weights,
                      struct.crit = struct.crit)
    }

    ## annotate graph (same order as in old version of this function)
    g$E$weights         <- sapply(out, function(x) x$w)
    g$E$names           <- sapply(out, function(x) x$name)
    g$E$conditionedSet  <- lapply(out, function(x) x$nedSet)
    g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
    g$E$todel           <- sapply(out, function(x) x$todel)

    ## delete edges that are prohibited by the proximity condition
    deleteEdges(g)
}

## function for obtaining edge information
getEdgeInfo2 <- function(i, g, oldVineGraph, weights, struct.crit = "tau") {

    ## get edge
    con <- g$E$nums[i, ]
    temp <- oldVineGraph$E$nums[con, ]

    ## check for proximity condition
    ok <- FALSE
    if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
        ok <- TRUE
        same <- temp[2, 1]
    } else {
        if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
            ok <- TRUE
            same <- temp[2, 2]
        }
    }

    ## dummy output
    w <- nedSet <- ningSet <- name <- NA
    todel <- TRUE

    # info if proximity condition is fulfilled ...
    if (ok) {
        ## get data
        if (temp[1, 1] == same) {
            zr1 <- oldVineGraph$E$Copula.CondData.2[[con[1]]]
        } else {
            zr1 <- oldVineGraph$E$Copula.CondData.1[[con[1]]]
        }
        if (temp[2, 1] == same) {
            zr2 <- oldVineGraph$E$Copula.CondData.2[[con[2]]]
        } else {
            zr2 <- oldVineGraph$E$Copula.CondData.1[[con[2]]]
        }
        zr1a <- if (is.list(zr1)) as.vector(zr1[[1]]) else zr1
        zr2a <- if (is.list(zr2)) as.vector(zr2[[1]]) else zr2

        ## calculate Kendall's tau
        keine_nas <- !(is.na(zr1a) | is.na(zr2a))
        if (struct.crit == "tau") {
            w <- fasttau(zr1a[keine_nas],
                                      zr2a[keine_nas],
                                      weights)
        } else if (struct.crit == "AIC") {
            w <- kdecop(cbind(zr1a[keine_nas], zr2a[keine_nas]),
                                     info = TRUE)$info$AIC
        } else if (struct.crit == "cAIC") {
            w <- kdecop(cbind(zr1a[keine_nas], zr2a[keine_nas]),
                                     info = TRUE)$info$cAIC
        }  else if (struct.crit == "hoeffd") {
            w <- abs(hoeffd(cbind(zr1a[keine_nas], zr2a[keine_nas])))
        }

        ## get names
        name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
        name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]

        ## infer conditioned set and conditioning set
        l1 <- c(g$V$conditionedSet[[con[1]]],
                g$V$conditioningSet[[con[1]]])
        l2 <- c(g$V$conditionedSet[[con[2]]],
                g$V$conditioningSet[[con[2]]])
        nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
        ningSet <- intersect(l1, l2)

        ## set edge name
        nmdiff <- c(setdiff(name.node1, name.node2),
                    setdiff(name.node2, name.node1))
        nmsect <- intersect(name.node1, name.node2)
        name <- paste(paste(nmdiff, collapse = ","),
                      paste(nmsect, collapse = ","),
                      sep = " ; ")

        ## mark as ok
        todel <- FALSE
    }

    ## return edge information
    list(w = w,
         nedSet = nedSet,
         ningSet = ningSet,
         name = name,
         todel = todel)
}

## functions for handling the tree structure -------------------------
graphFromTauMatrix <- function(tau) {
    d <- ncol(tau)
    # get variable names
    nms <- colnames(tau)
    # construct edge set
    E <- cbind(do.call(c, sapply(1:(d-1), function(i) seq.int(i))),
               do.call(c, sapply(1:(d-1), function(i) rep(i+1, i))))
    # add edge names
    E.names <- apply(E, 1, function(x) paste(nms[x[1]],  nms[x[2]], sep = ","))
    # set weights
    w <- tau[upper.tri(tau)]

    ## return results
    list(V = list(names = nms,
                  conditionedSet = NULL,
                  conditioningSet = NULL),
         E = list(nums = E,
                  names = E.names,
                  weights = w,
                  conditionedSet = lapply(1:nrow(E), function(i) E[i, ]),
                  conditioningSet = NULL))
}



makeFullGraph2 <- function(d) {
    ## create matrix of all combinations
    E <- cbind(do.call(c, lapply(1:(d-1), function(i) rep(i, d-i))),
               do.call(c, lapply(1:(d-1), function(i) (i+1):d)))
    E <- matrix(E, ncol = 2)

    ## output dummy list with edges set
    list(V = list(names = NULL,
                  conditionedSet = NULL,
                  conditioningSet = NULL),
         E = list(nums = E,
                  names = NULL,
                  weights = NULL,
                  conditionedSet = E,
                  conditioningSet = NULL))
}

findMaxTree2 <- function(g, mode = "RVine") {
    ## construct adjency matrix
    A <- adjacencyMatrix(g)
    d <- ncol(A)

    if (mode == "RVine") {
        ## initialize
        tree <- NULL
        edges <- matrix(NA, d - 1, 2)
        w <- numeric(d - 1)
        i <- 1

        ## construct minimum spanning tree
        for (k in 1:(d - 1)) {
            # add selected edge to tree
            tree <- c(tree, i)

            # find edge with minimal weight
            m <- apply(as.matrix(A[, tree]), 2, min)
            a <- apply(as.matrix(A[, tree]), 2, function(x) order(rank(x)))[1, ]
            b <- order(rank(m))[1]
            j <- tree[b]
            i <- a[b]

            # store edge and weight
            edges[k, ] <- c(j, i)
            w[k] <- A[i, j]

            ## adjust adjecency matrix to prevent loops
            for (t in tree)
                A[i, t] <- A[t, i] <- Inf
        }

        ## reorder edges for backwads compatibility with igraph output
        edges <- t(apply(edges, 1, function(x) sort(x)))
        edges <- edges[order(edges[, 2], edges[, 1]), ]

        ## delete unused edges from graph
        E <- g$E$nums
        in.tree <- apply(matrix(edges, ncol = 2), 1,
                         function(x) which((x[1] == E[, 1]) & (x[2] == E[, 2])))
        if (is.list(in.tree))
            in.tree <- unlist(in.tree)
        MST <- g
        g$E$todel <- rep(TRUE, nrow(E))
        if (any(g$E$todel)) {
            g$E$todel[in.tree] <- FALSE
            MST <- deleteEdges(g)
        }
    } else if (mode  == "CVine") {
        ## set root as vertex with minimal sum of weights
        A <- adjacencyMatrix(g)
        diag(A) <- 0
        sumtaus <- rowSums(A)
        root <- which.min(sumtaus)

        ## delete unused edges
        g$E$todel <- !((g$E$nums[, 2] == root) | (g$E$nums[, 1] == root))
        MST <- g
        if (any(g$E$todel ))
            MST <- deleteEdges(g)
    } else {
        stop("vine not implemented")
    }

    ## return result
    MST
}

adjacencyMatrix <- function(g) {
    ## create matrix of all combinations
    d <- length(g$V$names)
    v.all <- cbind(do.call(c, lapply(1:(d-1), function(i) seq.int(i))),
                   do.call(c, lapply(1:(d-1), function(i) rep(i+1, i))))

    ## fnd weight
    vals <- apply(v.all, 1, set_weight, E = g$E)

    ## create symmetric matrix of weights
    M <- matrix(0, d, d)
    M[upper.tri(M)] <- vals
    M <- M + t(M)
    diag(M) <- Inf

    ## return final matrix
    M
}

set_weight <- function(x, E) {
    is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
    w <- if (!any(is.edge)) Inf else (1 - abs(E$weights[which(is.edge)]))
    w
}


deleteEdges <- function(g) {
    ## reduce edge list
    keep <- which(!g$E$todel)
    E <- list(nums            = matrix(g$E$nums[keep, ], ncol = 2),
              names           = g$E$names[keep],
              weights         = g$E$weights[keep],
              conditionedSet  = g$E$conditionedSet[keep],
              conditioningSet = g$E$conditioningSet[keep])

    ## return reduced graph
    list(V = g$V, E = E)
}

# findMaximumAICTree <- function(g, mode = "RVine") {
#     if (mode == "RVine") {
#         return(minimum.spanning.tree(g, weights = -E(g)$weight))
#     } else if (mode == "CVine") {
#         stop("'struct.crit = AIC' is not yet available for C-Vines")
#     }
# }


# findMaxTree <- function(g, mode = "RVine", struct.crit = "tau") {
#     switch(struct.crit,
#            "tau"  = findMaximumTauTree(g, mode = mode),
#            "AIC"  = findMaximumAICTree(g, mode = mode),
#            "cAIC" = findMaximumAICTree(g, mode = mode))
# }

fasttau <- function(x, y, weights = NA) {
    TauMatrix(cbind(x, y))[1, 2]
}


as.RVMKernel2 <- function(RVine) {

    n <- length(RVine$Tree) + 1
    nam <- RVine$Tree[[1]]$V$names
    nedSets <- list()

    ## get selected pairs
    for (k in 1:(n - 2)) {
        nedSets[[k]]    <- RVine$Tree[[k]]$E$conditionedSet
    }
    if (is.list(RVine$Tree[[n - 1]]$E$conditionedSet)) {
        nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet[[1]])
    } else {
        nedSets[[n - 1]] <- list(RVine$Tree[[n - 1]]$E$conditionedSet)
    }

    M <- matrix(NA,n,n)

    for (k in 1:(n-1)) {
        w <- nedSets[[n-k]][[1]][1]

        M[k, k] <- w
        M[k+1, k] <- nedSets[[n-k]][[1]][2]

        if (k == (n-1)) {
            M[(k+1),(k+1)] <- nedSets[[n-k]][[1]][2]
        } else {
            for (i in (k+2):n) {
                for (j in 1:length(nedSets[[n-i+1]])) {
                    cs <- nedSets[[n-i+1]][[j]]
                    if (cs[1] == w) {
                        M[i, k] <- cs[2]
                        break
                    } else if (cs[2] == w) {
                        M[i, k] <- cs[1]
                        break
                    }
                }
                nedSets[[n-i+1]][[j]] <- NULL
            }
        }
    }

    M[is.na(M)] <- 0
    zrs <- matrix(0, n, n)

    RVineMatrix(M, family = zrs, par = zrs, par2 = zrs, names = nam)
}

