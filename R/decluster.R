#' A declustering scheme for dependent exceedances
#'
#' @description
#' Calculate clusters of exceedances for a sequence of dependent random variables
#'
#' @param x         Numeric, sequence of depedent random variables.
#' @param index     Numeric, vector specifying the index for the elements of x.
#' @param uq        Numeric in (0, 1), calculates the threshold as quantile(x, uq).
#' @param r         Numeric, nonnegative integer, the decluster parameter.
#' @param u         Numeric, if specified, the threshold is u and overrides what is
#'                  chosen for uq.
#' @param doplot    Logical, should a plot be generated?
#' @param ...       Additional arguments passed to the plot that is generated if
#'                  doplot = TRUE.
#'
#' @details
#' The argument index may be used to specify when each observation in x has occurred.
#' For example, suppose x contains daily observations in January on two separate years,
#' making length(x) = 62. Leaving index blank would treat January 31 of year one as
#' the day just before January 1 of year two, allowing the possibility of clusters to
#' exist that span both year. In this example, index should be c(1:31, 366:396).
#'
#' The decluster parameter r controls the cluster sizes. Two exceedances are part of
#' different clusters if there are at least r non-exceedances between the two obsverations.
#' This has the consequence that for r = 0, every observation will be its own cluster.
#'
#' When multiple observations are part of the same cluster, the maximum value in that
#' cluster is the value that is return.
#'
#' If doplot = TRUE (the default), then a plot of (index, x) is shown with a horizontal line
#' at the threshold, and clusters are marked by the shape and color. (Note: the color/shape
#' pattern of clusters is repeated every 77 clusters).
#'
#' @return
#' The maximum _excess_ for each cluster, the number of clusters, the number of exceedances,
#' the number of observations, the threshold used, the decluster parameter
#'
#' @export
#' @examples
#' set.seed(1)
#' x = rnorm(100)
#' decluster(x, uq = 0.5, r=0) # 50 clusters
#' decluster(x, uq = 0.5, r=1) # 26 clusters
#' decluster(x, uq = 0.5, r=2) # 12 clusters
#' decluster(x, uq = 0.5, r=3) # 4 clusters
#' decluster(x, uq = 0.5, r=4) # 4 clusters
#' decluster(x, uq = 0.5, r=5) # 3 clusters
#' decluster(x, uq = 0.5, r=6) # 2 clusters
#' decluster(x, uq = 0.5, r=7) # 1 cluster


decluster = function(x, index, uq = 0.90, r = 0, u, doplot = TRUE, ...){
    r = round(r)
    if (r < 0)
        stop("Number of observations below threshold between clusters, r, must be an integer at least 0.")
    if (uq <= 0 || uq >= 1)
        stop("uq must be in (0, 1).")
    if (missing(u)){
        threshold = quantile(x, uq)
    } else {
        threshold = u
        }
    if (missing(index))
        index = 1:length(x)
    ind1 = (which(c(x, c(rep(-Inf, r+1), Inf)) >= threshold))
    ind2 = which(diff(c(as.numeric(index), rep(Inf, r+2))[ind1]) > r)
    ind3 = apply(cbind(c(1, ind2[-length(ind2)]+1), ind2), 1,
        function(x) seq(x[1], x[2], by = 1))

    if (r != 0 && !is.list(ind3))
        ind3 = list(c(ind3))

    n_c = length(ind3)    # number of clusters
    n_u = length(ind1)-1  # number of observations above threshold
    n = length(x)

    cols = rep(rgb(0.7,0.7,0.7,0.3), length(x))
    pchs = rep(16, length(x))

    y = double(n_c)
    for (j in 1:n_c){
        tt = which(x >= threshold)[ind3[[j]]]
        cols[tt] = ((j-1) %% 7)+1
        pchs[tt] = ((j-1) %% 11) + (((j-1) %% 11) < 4)*15 - (((j-1) %% 11) > 4)*4
        y[j] = max(x[tt] - threshold)
        }

    if (doplot){
        plot(index, x, pch = pchs, col = cols,  ...)
        lines(index, x, col = rgb(0.7, 0.7, 0.7, 0.3))
        abline(h = threshold)
        }

    return (list("max_cluster_excess"=y, "n_clust"=n_c, "n_u"=n_u, "n"=n, "threshold"=threshold, "r"=r))
    }


