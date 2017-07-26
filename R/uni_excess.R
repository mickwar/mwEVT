#' Univariate threshold excess model
#'
#' @description
#' Fit a standard univariate threshold exceedance model.
#'
#' @param x                 Numeric vector of data.
#' @param uq                numeric, within (0, 1), the quantile of the complete data
#'                          (includes replicates if available) that determines the
#'                          threshold
#' @param threshold         numeric, values greater than threshold are considered
#'                          excess. Defaults to the threshold given by uq, but if
#'                          specified, uq is ignored.
#' @param r                 Numeric, nonnegative integer, specifies the degree
#'                          of separation between clusters in a declustering
#'                          scheme. Defaults to r = 0 (no clustering).
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

uni_excess = function(x, uq = 0.95, threshold, r = 0, m.ksi = 0, s.ksi = 10,
    nburn = 10000, nmcmc = 10000, window = 200, chain_init){

    require(mwBASE)

    if (missing(threshold))
        threshold = quantile(x, uq)

    if (missing(chain_init))
        chain_init = c(1, 1e-6)

    dat = NULL
    tmp = decluster(x = x, r = r, doplot = FALSE)
    dat$y = tmp$max_cluster_excess
    dat$n_c = tmp$n_clust
    dat$n_u = tmp$n_u
    dat$n = tmp$n
    rm(tmp)

    ### MCMC
    calc.post = function(dat, params){
        # x is a list where x$y is the exceedances, x$n_c is the number of clusters
        sigma = params[1]
        ksi = params[2]
        # (lower) boundary check
        if (any(1 + ksi*dat$y/sigma < 0))
            return (-Inf)
        # Likelihood
        if (ksi != 0){
            out = -dat$n_c*log(sigma) - (1 + 1/ksi)*sum(log(1+ksi*dat$y/sigma))
        } else {
            out = -dat$n_c*log(sigma) - sum(dat$y)/sigma
            }
        # Priors
        out = out - log(sigma)
        out = out + dnorm(ksi, m.ksi, s.ksi, log = TRUE)
        return (out)
        }

    out = mcmc_sampler(dat, calc.post, 2, nmcmc = nmcmc, nburn = nburn,
        nthin = 1, window = window, bounds = list("lower" = c(0, -Inf),
        "upper" = c(Inf, Inf)), chain_init = chain_init)


    mean(out$accept)


    out$hpds = lapply(apply(out$params, 2, function(x) list(range(hpd_mult(x)))), "[[", 1)

    ### Calculations for diagnostic plots
    message("Computing return levels and other diagnostics ...")
    tmp = diag_computations(out, llu = 100, additional.return.periods = c(20, 30, 50))
    out$return.period = tmp$return.period
    out$lm1           = tmp$lm1
    out$lq1           = tmp$lq1
    out$lm2           = tmp$lm2
    out$lq2           = tmp$lq2
    out$Zm            = tmp$Zm
    out$Zq            = tmp$Zq              
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
