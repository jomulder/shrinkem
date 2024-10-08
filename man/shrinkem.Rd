% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrinkem_main_functions.R
\name{shrinkem}
\alias{shrinkem}
\title{Fast Bayesian regularization using Gaussian approximations}
\usage{
shrinkem(
  x,
  Sigma,
  type,
  group,
  iterations,
  burnin,
  store,
  cred.level,
  df1,
  df2,
  scale2,
  lambda2.fixed,
  lambda2,
  ...
)
}
\arguments{
\item{x}{A vector of estimates.}

\item{Sigma}{A covariance matrix capturing the uncertainty of the estimates (e.g., error covariance matrix).}

\item{type}{A character string which specifies the type of regularization method is used. Currently, the
types "ridge", "lasso", and "horseshoe", are supported.}

\item{group}{A vector of integers denoting the group membership of the estimates, where each group receives
a different global shrinkage parameter which is adapted to the observed data.}

\item{iterations}{Number of posterior draws after burnin. Default = 5e4.}

\item{burnin}{Number of posterior draws in burnin. Default = 1e3.}

\item{store}{Store every store-th draw from posterior. Default = 1 (implying that every draw is stored).}

\item{cred.level}{The significance level that is used to check whether a parameter is nonzero depending on whether
0 is contained in the credible interval. The default is \code{cred.level = 0.95}.}

\item{df1}{First hyperparameter (degrees of freedom) of the prior for a shrinkage parameter lambda^2, which follows a F(df1,df2,scale2)
distribution. The default is \code{df1 = 1}. For \code{df1 = 1}, this corresponds to half-t distribution for lambda with degrees of freedom \code{df2}
and scale parameter \code{sqrt(scale2/df2)}.}

\item{df2}{Second hyperparameter (degrees of freedom) of the prior for a shrinkage parameter lambda^2, which follows a F(df1,df2,scale2)
distribution. The default is \code{df2 = 1}.}

\item{scale2}{Second hyperparameter (scale parameter) of the prior for a shrinkage parameter lambda^2, which follows a F(df1,df2,scale2)
distribution. The default is \code{df2 = 1e3}.}

\item{lambda2.fixed}{Logical indicating whether the penalty parameters(s) is/are fixed. Default is FALSE.}

\item{lambda2}{Positive scalars of length equal to the number of groups in 'group'. The argument is only
used if the argument 'lambda2.fixed' is 'TRUE'.}

\item{...}{Parameters passed to and from other functions.}
}
\value{
The output is an object of class \code{shrinkem}. The object has elements:
\itemize{
\item \code{estimates}: A data frame with the input estimates, the shrunken posterior mean, median, and mode,
the lower and upperbound of the credbility interval based on the shrunken posterior, and a logical which indicates if
zero is contained in the credibility interval.
\item \code{draws}: List containing the posterior draws of the effects (\code{beta}), the prior parameters (\code{tau2}, \code{gamma2}),
and the penalty parameters (\code{psi2} and \code{lambda2}).
\item \code{dim.est}: The dimension of the input estimates of \code{beta}.
\item \code{input.est}: The input vector of the unshrunken estimates of \code{beta}.
\item \code{call}: Input call.
}
}
\description{
The \code{shrinkem} function can be used for regularizing a vector
of estimates using Bayesian shrinkage methods where the uncertainty of the estimates
are assumed to follow a Gaussian distribution.
}
\examples{
\donttest{
# EXAMPLE
estimates <- -5:5
covmatrix <- diag(11)
# Bayesian horseshoe where all beta's have the same global shrinkage
# (using default 'group' argument)
shrink1 <- shrinkem(estimates, covmatrix, type="horseshoe")
# posterior modes of middle three estimates are practically zero

# plot posterior densities
old.par.mfrow <- par(mfrow = c(1,1))
old.par.mar <- par(mar = c(0, 0, 0, 0))
par(mfrow = c(11,1))
par(mar = c(1,2,1,2))
for(p in 1:ncol(shrink1$draws$beta)){plot(density(shrink1$draws$beta[,p]),
  xlim=c(-10,10),main=colnames(shrink1$draws$beta)[p])}
par(mfrow = old.par.mfrow)
par(mar = old.par.mar)

# Bayesian horseshoe where first three and last three beta's have different
# global shrinkage parameter than other beta's
shrink2 <- shrinkem(estimates, covmatrix, type="horseshoe",
   group=c(rep(1,3),rep(2,5),rep(1,3)))
# posterior modes of middle five estimates are virtually zero

# plot posterior densities
par(mfrow = c(11,1))
par(mar = c(1,2,1,2))
for(p in 1:ncol(shrink2$draws$beta)){plot(density(shrink2$draws$beta[,p]),xlim=c(-10,10),
  main=colnames(shrink2$draws$beta)[p])}
par(mfrow = old.par.mfrow)
par(mar = old.par.mar)
}

}
\references{
Karimovo, van Erp, Leenders, and Mulder (2024). Honey, I Shrunk the Irrelevant Effects! Simple and Fast Approximate Bayesian
Regularization. <https://doi.org/10.31234/osf.io/2g8qm>
}
