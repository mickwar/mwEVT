theta_hier = function(y, u, n, prior, ...){
    require(mwBASE)
    exceedance = apply(y, 2, function(x) which(x > u))

    if (missing(n))
        n = NROW(y)

    R = NCOL(y)

    Tu = lapply(exceedance, diff)
    N = sapply(exceedance, length)
    m1 = sapply(Tu, function(x) sum(x == 1))
    emp.p = mean(y <= u)

    dat = list("y" = y, "exceed" = exceedance, 
        "Tu" = Tu, "N" = N, "m1" = m1, "emp.p" = emp.p)

    if (missing(prior))
        prior = list("theta_a" = 1/2, "theta_b" = 1,
            "p_a" = emp.p*100, "p_b" = (1-emp.p)*100,
            "nu_a" = 1, "nu_b" = 1/10,
            "tau_a" = 1, "tau_b" = 1/10)

    calc.post = function(x, param){
        R = length(x$exceed)
        
        # theta_i, p_i for each ``break''
        theta_i = param[1:R]
        p_i = param[R + (1:R)]

        theta = param[2*R + 1]
        p = param[2*R + 2]
        nu = param[2*R + 3]
        tau = param[2*R + 4]

        if (any(theta_i <= 0 | theta_i > 1))
            return (-Inf)
        if (any(p_i <= 0 | p_i >= 1))
            return (-Inf)
        if (theta <= 0 || theta > 1)
            return (-Inf)
        if (p <= 0 || p > 1)
            return (-Inf)
        if (nu <= 0)
            return (-Inf)
        if (tau <= 0)
            return (-Inf)

        m1 = x$m1

        # Likelihood (from Eq. 3)
        out = sum(m1 * log(1 - theta_i*p_i^theta_i) + 
            ifelse(N - 1 - m1 >= 0, N - 1 - m1, 0) * (log(theta_i) + log(1-p_i^theta_i)) + 
            theta_i*log(p_i) * sapply(x$Tu, function(x) sum(x - 1)))

        # Priors
        out = out + sum(dbeta(theta_i, theta*nu, (1-theta)*nu, log = TRUE))
        out = out + sum(dbeta(p_i, p*tau, (1-p)*tau, log = TRUE))
        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
        out = out + dbeta(p, prior$p_a, prior$p_b, log = TRUE)
        out = out + dgamma(nu, prior$nu_a, prior$nu_b, log = TRUE)
        out = out + dgamma(tau, prior$tau_a, prior$tau_b, log = TRUE)
        return (out)
        }
    mcmc_out = mcmc_sampler(data = dat, target = calc.post, nparam = 2*R + 4,
        groups = list(1:R, R + (1:R), 2*R + (1:4)), ...)

    theta.hat = colMeans(mcmc_out$param[,1:R])


    ### Decluster
    C = floor(theta.hat*N)+1

    for (j in 1:R)
        C[j] = min(C[j], N[j]-1)

    tmp = lapply(Tu, sort, decreasing = TRUE)

    T_C = double(R)
    for (j in 1:R){
        T_C[j] = tmp[[j]][C[j]]
        while (!(tmp[[j]][C[j]-1] > T_C[j]) && (C[j] > 1)){
            C[j] = C[j] - 1
            T_C[j] = tmp[[j]][C[j]]
            }
        }

    ind.obs = rep(list(NULL), R)

    for (j in 1:R){
        # The set of independent intercluster times
        inter.Clust = Tu[[j]][Tu[[j]] > T_C[j]]

        i_j = which(Tu[[j]] > T_C[j])
        i_j = c(0, i_j, N[j])
        ind.seq = rep(list(NULL), C[j])
        intra.Clust = rep(list(NULL), C[j])    # The interexceedance times within each cluster
        nice.S = rep(list(NULL), C[j])     # List of independent clusters, marking when exceedances occur
        nice.C = rep(list(NULL), C[j])     # The observed value at the exceedance times

        for (k in 2:(C[j]+1)){
            ind.seq[[k-1]] = seq(i_j[k-1]+1, i_j[k]-1)
            if (i_j[k-1]+1 == i_j[k]){
        #       nice.T[[j-1]] = NULL
            } else {
                intra.Clust[[k-1]] = Tu[[j]][seq(i_j[k-1]+1, i_j[k]-1)]
                }
            nice.S[[k-1]] = exceedance[[j]][seq(i_j[k-1]+1, i_j[k])]
            nice.C[[k-1]] = y[nice.S[[k-1]], j]
            }

        ind.obs[[j]] = sapply(nice.C, max) # Get the greatest value within each (independent) cluster
        }

    par = mcmc_out$params[,c(1:R, 2*R+1)]

    return (list("y" = ind.obs, "N" = N, "T_C" = T_C, "mcmc" = par))
    }