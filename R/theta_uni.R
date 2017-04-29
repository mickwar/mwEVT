#' Extremal index estimation based on Ferro and Segers (2003)
#'
#' @description
#' Uses the likelihood provided in Ferro and Segers (2003) in a Bayesian setting
#' to obtain posterior samples of the extremal index. Separates exceedances into
#' clusters.
#'
#' @param y         Numeric, sequence of depedent random variables. May be a matrix.
#'                  See details.
#' @param u         Numeric, the threshold. Defaults to quantile(y, 0.90).
#' @param ord       Numeric, vector of length R = NCOL(y). Defaults to 1:R. See details.
#' @param prior     List.
#' @param ...       Additional arguments passed to mwBASE::mcmc_sampler.
#'
#' @details
#' The data y is expected to be observed values along an evenly-spaced grid of width 1.
#' If y is an n x R matrix, then it will be collapsed to y = c(y[,ord]). This may be
#' the case when working with multiple realizations from the same process or with
#' specific seasons on a time-series.
#'
#' @export
#' @examples
#' ### Generate dependent random variables with standard Frechet marginals
#' set.seed(1)
#' n = 365*50     # 50 years of data, no leap years
#'
#' theta = 0.25
#' W = -1/log(runif(n))
#' y = double(n)
#' y[1] = W[1] / theta
#' for (i in 2:n)
#'     y[i] = max((1-theta)*y[i-1], W[i])
#'
#' out = theta_uni(y, quantile(y, 0.95), nburn = 40000, nmcmc = 20000)
#'
#' mean(out$mcmc)     # 0.277
#' quantile(out$mcmc, c(0.025, 0.975))    # c(0.248, 0.306)
#'
#'
#' ### Generate dependent random variables with standard Frechet marginals,
#' ### but only for winter, say.
#' set.seed(1)
#' R = 50   # number of seasons (years)
#' n = 90   # number of observations (days) per season
#'
#' theta = 0.25
#' W = matrix(-1/log(runif(n*R)), n, R)
#' y = matrix(0, n, R)
#' y[1,] = W[1,] / theta
#' for (i in 2:n)
#'     y[i,] = apply(rbind((1-theta)*y[i-1,], W[i,]), 2, max)
#'
#' out = theta_uni(y, quantile(y, 0.95), nburn = 40000, nmcmc = 20000)
#'
#' mean(out$mcmc)     # 0.267
#' quantile(out$mcmc, c(0.025, 0.975))    # c(0.214, 0.325)

theta_uni = function(y, u, ord, prior, ...){
    require(mwBASE)

    n = NROW(y)
    R = NCOL(y)

    # ``ord'' is for re-ordering the columns when vec'ing
    # (Could result in different interexceedance times)
    if (missing(ord))
        ord = 1:R

    # Vec the matrix (if applicable)
    if (R > 1)
        y = c(y[,ord])

    # Default threshold to the 0.90 quantile
    if (missing(u))
        u = quantile(y, 0.90)

    
    exceed = which(y > u)
    Tu = diff(exceed)
    N = length(exceed)
    m1 = sum(Tu == 1)
    emp.p = mean(y <= u)

    # If multiple realizations, adjust m1 for the possibility of truncating
    # clusters of extremes
#   if (R > 1)
#       m1 = m1 + (sum(which(x > u) %% n == 0) + sum(which(x > u) %% n == 1)) / 2

    if (missing(prior))
        prior = list("theta_a" = 1/2, "theta_b" = 1,
            "p_a" = emp.p*100, "p_b" = (1-emp.p)*100)

    calc.post = function(dat, param){
        theta = param[1]
        p = param[2]
        if (theta <= 0 || theta > 1)
            return (-Inf)
        if (p <= 0 || p >= 1)
            return (-Inf)

        # Likelihood (Ferro, Eq. 3)
        out = m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
            theta*log(p)*sum(dat - 1)

        # Priors
        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
        out = out + dbeta(p, prior$p_a, prior$p_b, log = TRUE)
        return (out)
        }
    mcmc_out = mcmc_sampler(data = Tu, target = calc.post, nparam = 2, ...)

    theta.hat = mean(mcmc_out$param[,1])


    ### Decluster
    C = floor(theta.hat*N)+1
    C = min(C, N-1)
    tmp = sort(Tu, decreasing = TRUE)
    T_C = tmp[C]
    while (!(tmp[C-1] > T_C) && (C > 1)){
        C = C - 1
        T_C = tmp[C]
        }

    # The set of independent intercluster times
    inter.Clust = Tu[Tu > T_C]

    i_j = which(Tu > T_C)
    i_j = c(0, i_j, N)
    ind.seq = rep(list(NULL), C)
    intra.Clust = rep(list(NULL), C)    # The interexceedance times within each cluster
    nice.S = rep(list(NULL), C)     # List of independent clusters, marking when exceedances occur
    nice.C = rep(list(NULL), C)     # The observed value at the exceedance times

    for (k in 2:(C+1)){
        ind.seq[[k-1]] = seq(i_j[k-1]+1, i_j[k]-1)
        if (i_j[k-1]+1 == i_j[k]){
    #       nice.T[[j-1]] = NULL
        } else {
            intra.Clust[[k-1]] = Tu[seq(i_j[k-1]+1, i_j[k]-1)]
            }
        nice.S[[k-1]] = exceed[seq(i_j[k-1]+1, i_j[k])]
        nice.C[[k-1]] = y[nice.S[[k-1]]]
        }

    # Get the greatest value within each (independent) cluster
    ind.obs = sapply(nice.C, max)

    return (list("y" = ind.obs, "N" = N, "T_C" = T_C, "mcmc" = mcmc_out$param[,1]))
    }
