% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta_uni.R
\name{theta_uni}
\alias{theta_uni}
\title{Extremal index estimation}
\usage{
theta_uni(y, u, ord, prior, likelihood, method, bsamp = 10000, ...)
}
\arguments{
\item{y}{Numeric, sequence of depedent random variables. May be a matrix.
See details.}

\item{u}{Numeric, the threshold. Defaults to quantile(y, 0.90).}

\item{ord}{Numeric, vector of length R = NCOL(y). Defaults to 1:R. See details.}

\item{prior}{List.}

\item{likelihood}{Character, either "ferro" or "suveges" (default).}

\item{method}{Character, either "classical" or "bayesian" (default).}

\item{bsamp}{Numeric, number of bootstrap samples for "classical" method.
Defaults to 10000.}

\item{...}{Additional arguments passed to mwBASE::mcmc_sampler.}
}
\description{
Uses the likelihood provided in Ferro and Segers (2003) in a Bayesian setting
to obtain posterior samples of the extremal index. Separates exceedances into
clusters.
}
\details{
The data y is expected to be observed values along an evenly-spaced grid of width 1.
If y is an n x R matrix, then it will be collapsed to y = c(y[,ord]). This may be
the case when working with multiple realizations from the same process or with
specific seasons on a time-series.

The 'classical' option for method produces the estimates given by either Ferro and
Segers (2003) or Suveges (2007). When method == 'bayesian' the likelihoods given in
the aforementioned papers are used to produce posterior samples for theta.
}
\examples{
### Generate dependent random variables with standard Frechet marginals
set.seed(1)
n = 365*50     # 50 years of data, no leap years

theta = 0.25
W = -1/log(runif(n))
y = double(n)
y[1] = W[1] / theta
for (i in 2:n)
    y[i] = max((1-theta)*y[i-1], W[i])

ferro_cl = theta_uni(y, quantile(y, 0.95), likelihood = "ferro",
    method = "classical")
suveges_cl = theta_uni(y, quantile(y, 0.95), likelihood = "suveges",
    method = "classical")
ferro_ba = theta_uni(y, quantile(y, 0.95), likelihood = "ferro",
    method = "bayesian", nburn = 40000, nmcmc = 20000)
suveges_ba = theta_uni(y, quantile(y, 0.95), likelihood = "suveges",
    method = "bayesian", nburn = 40000, nmcmc = 20000)

ferro_cl$theta      # 0.275
suveges_cl$theta    # 0.279

mean(ferro_ba$mcmc)     # 0.277
quantile(ferro_ba$mcmc, c(0.025, 0.975))    # c(0.248, 0.306)

mean(suveges_ba$mcmc)   # 0.279
quantile(suveges_ba$mcmc, c(0.025, 0.975))  # c(0.257, 0.301)


### Generate dependent random variables with standard Frechet marginals,
### but only for winter, say.
set.seed(1)
R = 50   # number of seasons (years)
n = 90   # number of observations (days) per season

theta = 0.25
W = matrix(-1/log(runif(n*R)), n, R)
y = matrix(0, n, R)
y[1,] = W[1,] / theta
for (i in 2:n)
    y[i,] = apply(rbind((1-theta)*y[i-1,], W[i,]), 2, max)

out = theta_uni(y, quantile(y, 0.95), nburn = 40000, nmcmc = 20000)

mean(out$mcmc)     # 0.267
quantile(out$mcmc, c(0.025, 0.975))    # c(0.214, 0.325)
}
