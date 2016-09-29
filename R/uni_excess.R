#' Univariate threshold excess model for climate observations and simulations
#'
#' @description
#' Fit a univariate threshold exceedance model for a processed data set. The
#' data should be in a matrix format where the first column is the time and
#' all subsequent columns contain the values of the time series.
#'
#' @param x                 character, the name of the data, e.g. "decadal1961",
#'                          "historical", or "obs19501999", corresponding to a
#'                          processed data set.
#' @param r                 numeric, nonnegative integer, specifies the degree
#'                          of separation between clusters in a declustering
#'                          scheme. A value of r=0 means each excess is given
#'                          its own cluster. If r is not specified, it will be
#'                          selected automatically.
#' @param uq                numeric, within (0, 1), the quantile of the complete data
#'                          (includes replicates if available) that determines the
#'                          threshold
#' @param threshold         numeric, values greater than threshold are considered
#'                          excess. Defaults to the threshold given by uq, but if
#'                          specified, uq is ignored.
#' @param m.ksi,s.ksi       numeric, prior mean and sd for ksi, ksi is normal
#' @param nburn             numeric, nonnegative integer, the first nburn samples
#'                          are discared. Defaults to 80000.
#' @param nmcmc             numeric, positive integer, the number of samples that are
#'                          kept _post_ burn-in. That is, a total of nburn + nmcmc
#'                          iterations of the chain are run. Defaults to 40000.
#' @param window            numeric, during the burn-in, the covariance matrix of the
#'                          proposal distribution is changed every window iterations
#'                          to improve the acceptance rate. Also, the iteration number
#'                          is displayed every window iterations. Defaults to 500.
#' @param chain_init        numeric vector of length 2, specifying the starting
#'                          location of the Markov chain for (sigma, ksi). Defaults to
#'                          c(1, 1e-6).
#'
#' @details
#'
#' @export
#'

uni_excess = function(x, uq = 0.90, threshold, r = 0, m.ksi = 0, s.ksi = 10,
    nburn = 10000, nmcmc = 10000, window = 200, chain_init){

   

    if (missing(threshold))
        threshold = quantile(x, uq)

    if (is.null(r))
        r = r.select(output)

    ## This group runs together
    message("Retrieving excesses ...")
    output$y = NULL
    output$n_c = 0
    output$n_u = 0
    output$n = 0
    for (j in 1:length(output$run.index)){
        tmp = decluster(x = output$varmat[,output$run.index[j]],
            times = output$time.dates, u = output$threshold,
            r = output$r, doplot = FALSE)
        output$y = c(output$y, tmp$max_cluster_excess)
        output$n_c = output$n_c + tmp$n_clust
        output$n_u = output$n_u + tmp$n_u
        output$n = output$n + tmp$n
        }
    rm(tmp)
    ## End group

    ### MCMC
    nparam = 2
    params = matrix(0, nburn + nmcmc, nparam)
    accept = double(nburn + nmcmc)
    cand.sig = diag(0.1, nparam)

    params[1,] = chain_start

    lower = c(0, -Inf)
    upper = c(Inf, Inf)

    calc.post = function(x, params){
        # x is a list where x$y is the exceedances, x$n_c is the number of clusters
        sigma = params[1]
        ksi = params[2]
        # (lower) boundary check
        if (any(1 + ksi*x$y/sigma < 0))
            return (-Inf)
        # Likelihood
        if (ksi != 0){
            out = -x$n_c*log(sigma) - (1 + 1/ksi)*sum(log(1+ksi*x$y/sigma))
        } else {
            out = -x$n_c*log(sigma) - sum(x$y)/sigma
            }
        # Priors
        out = out - log(sigma)
        out = out + dnorm(ksi, m.ksi, s.ksi, log = TRUE)
        return (out)
        }

    post = calc.post(output, params[1,])

    message("Obtaining posterior samples ...")
    for (i in 2:(nburn + nmcmc)){
        if (floor(i/window) == i/window)
            cat("\r   ", i, "/", nburn+nmcmc)
        params[i,] = params[i-1,]
        cand = mvrnorm(1, params[i-1,], cand.sig)
        if (all(cand > lower) && all(cand < upper)){
            cand.post = calc.post(output, cand)
            if (log(runif(1)) <= cand.post - post){
                post = cand.post
                params[i,] = cand
                accept[i] = 1
                }
            }
        if ((floor(i/window) == i/window) && (i <= nburn))
            cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window/50) *
                (cand.sig + window * var(params[(i-window+1):i,]) / i)
        if (i == (nburn + nmcmc))
            cat("\n")
        }

    params = tail(params, nmcmc)
    accept = tail(accept, nmcmc)

    ### Beta(1/2, 1/2) priors assumed for zeta and theta
    # output$n_u/output$n   # zeta mle
    zeta = rbeta(nmcmc, output$n_u + 1/2, output$n - output$n_u + 1/2)

    # output$n_c/output$n_u # theta mle
    theta = rbeta(nmcmc, output$n_c + 1/2, output$n_u - output$n_c + 1/2)

    params = cbind(params, zeta, theta)
    colnames(params)[1:2] = c("sigma", "ksi")

    output$params = params
    output$accept = accept
    output$cand.sig = cand.sig

    output$hpds = lapply(apply(params, 2, function(x) list(range(hpd_mult(x)))), "[[", 1)

    ### Calculations for diagnostic plots
    message("Computing return levels and other diagnostics ...")
    tmp = diag_computations(output, llu = 100, additional.return.periods = c(20, 30, 50))
    output$return.period = tmp$return.period
    output$lm1           = tmp$lm1
    output$lq1           = tmp$lq1
    output$lm2           = tmp$lm2
    output$lq2           = tmp$lq2
    output$Zm            = tmp$Zm
    output$Zq            = tmp$Zq              
    rm(tmp)

    ### Compute DIC
    d.theta = apply(output$params, 1, function(x) -2*sum(log(dgpd(output$y, 0, x[1], x[2])/length(output$y))))
    output$DIC = mean(d.theta) + 0.5*var(d.theta)

    obj.name = paste0(output$region, "_", output$variable, "_", output$model,
        "_", ifelse(output$anomaly, "anom_", ""), output$season, "_",
        gsub("-", "", date_begin), "_", gsub("-", "", date_end))
    assign(obj.name, output)

    if (!dir.exists("./RData/"))
        dir.create("./RData/")
    save(list = obj.name, file = paste0("./RData/", obj.name, ".RData"))
    message("Note: R object saved in:")
    message(paste0("    ", getwd(), "/RData/", obj.name, ".RData\n"))

    if (return_output)
        return (output)
    }
