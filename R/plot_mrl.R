#' Generate a mean residual life plot
#'
#' @description
#' Plots the mean residual life for a vector of data, with 95 confidence bands
#'
#' @param x     Numeric, vector.
#' @param qs    Numeric, vector quantiles which will be marked on the x-axis. This
#'              does not affect the computation of the mean residual life.
#' @param noopt Logical, should the plot function be called with no additional options?
#'              Defaults to FALSE. This would need to be TRUE for full control over the plot.
#' @param ...   Additional arguments passed to plot().
#'
#' @details
#' If the generalized Pareto is a valid model for the tail of a distribution (i.e. given
#' a threshold), then the mean residual life plot should be approximately linear beyond
#' the threshold.
#'
#' The motivation for using this plot is in threshold selection. By using the plot, one
#' can observe the lowest threshold for which the generalized Pareto may be valid.
#'
#' If noopt = TRUE, anything regarding the plot can be specified by the user except ylim
#' and the lines.
#'
#' @export
#' @examples
#' plot_mrl(rnorm(1000))
#' plot_mrl(rgamma(1000, 1, 1))

plot_mrl = function(x, qs = c(0, 0.5, 0.75, seq(0.8, 1, by = 0.01)), noopt = FALSE, ...){
    u = head(sort(x), length(x)-2)
    mrl = double(length(u))
    mrlqq = double(length(u))
    for (i in 1:length(u)){
        excess = (x[which(x > u[i])] - u[i])
        mrl[i] = mean(excess)
        mrlqq[i] = qnorm(0.975) * sd(excess) / sqrt(length(excess))
        }
    mrl = ifelse(is.na(mrl), 0, mrl)
    mrlqq = ifelse(is.na(mrlqq), 0, mrlqq)

    if (noopt){
        plot(u, mrl, type='l', ylim = c(min(mrl-mrlqq), max(mrl+mrlqq)), ...)
        lines(u, mrl-mrlqq, type='l', col='gray50')
        lines(u, mrl+mrlqq, type='l', col='gray50')
    } else {
        plot(u, mrl, type='l', ylim = c(min(mrl-mrlqq), max(mrl+mrlqq)), axes = FALSE,
            main = "Mean residual life plot", xlab = "Threshold (u) Quantile",
            ylab = "Mean excess, m(u)", ...)
        lines(u, mrl-mrlqq, type='l', col='gray50')
        lines(u, mrl+mrlqq, type='l', col='gray50')
        axis(2)
        axis(1, at = quantile(x, qs), label = as.character(qs))
        }
    }
