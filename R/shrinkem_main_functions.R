#' @title Fast Bayesian regularization using Gaussian approximations
#'
#' @description The \code{shrinkem} function can be used for regularizing a vector
#' of estimates using Bayesian shrinkage methods where the uncertainty of the estimates
#' are assumed to follow a Gaussian distribution.
#'
#' @param x A vector of estimates.
#' @param Sigma A covariance matrix capturing the uncertainty of the estimates (e.g., error covariance matrix).
#' @param type A character string which specifies the type of regularization method is used. Currently, the
#' types "ridge", "lasso", and "horseshoe", are supported.
#' @param group A vector of integers denoting the group membership of the estimates, where each group receives
#' a different global shrinkage parameter which is adapted to the observed data.
#' @param cred.level The significance level that is used to check whether a parameter is nonzero depending on whether
#' 0 is contained in the credible interval. The default is \code{cred.level = 0.95}.
#' @param df1 First hyperparameter (degrees of freedom) of the prior for a shrinkage parameter lambda^2, which follows a F(df1,df2,scale2)
#' distribution. The default is \code{df1 = 1}. For \code{df1 = 1}, this corresponds to half-t distribution for lambda with degrees of freedom \code{df2}
#' and scale parameter \code{sqrt(scale2/df2)}.
#' @param df2 Second hyperparameter (degrees of freedom) of the prior for a shrinkage parameter lambda^2, which follows a F(df1,df2,scale2)
#' distribution. The default is \code{df2 = 1}.
#' @param scale2 Second hyperparameter (scale parameter) of the prior for a shrinkage parameter lambda^2, which follows a F(df1,df2,scale2)
#' distribution. The default is \code{scale2 = 1e3}.
#' @param iterations Number of posterior draws after burnin. Default = 5e4.
#' @param burnin Number of posterior draws in burnin. Default = 1e3.
#' @param store Store every store-th draw from posterior. Default = 1 (implying that every draw is stored).
#' @param lambda2.fixed Logical indicating whether the penalty parameters(s) is/are fixed. Default is FALSE.
#' @param lambda2 Positive scalars of length equal to the number of groups in 'group'. The argument is only
#' used if the argument 'lambda2.fixed' is 'TRUE'.
#' @param nugget A small positive value close to 0 which is used to avoid numerically singular matrices.
#' The default is \code{1e-8}.
#' @param ... Parameters passed to and from other functions.
#' @return The output is an object of class \code{shrinkem}. The object has elements:
#' \itemize{
#' \item \code{estimates}: A data frame with the input estimates, the shrunken posterior mean, median, and mode,
#' the lower and upperbound of the credbility interval based on the shrunken posterior, and a logical which indicates if
#' zero is contained in the credibility interval.
#' \item \code{draws}: List containing the posterior draws of the effects (\code{beta}), the prior parameters (\code{tau2}, \code{gamma2}),
#' and the penalty parameters (\code{psi2} and \code{lambda2}).
#' \item \code{dim.est}: The dimension of the input estimates of \code{beta}.
#' \item \code{input.est}: The input vector of the unshrunken estimates of \code{beta}.
#' \item \code{call}: Input call.
#' }
#' @rdname shrinkem
#' @references Karimovo, van Erp, Leenders, and Mulder (2024). Honey, I Shrunk the Irrelevant Effects! Simple and Fast Approximate Bayesian
#' Regularization. <https://doi.org/10.31234/osf.io/2g8qm>
#' @examples
#' \donttest{
#' # EXAMPLE
#' estimates <- -5:5
#' covmatrix <- diag(11)
#' # Bayesian horseshoe where all beta's have the same global shrinkage
#' # (using default 'group' argument)
#' shrink1 <- shrinkem(estimates, covmatrix, type="horseshoe")
#' # posterior modes of middle three estimates are practically zero
#'
#' # plot posterior densities
#' old.par.mfrow <- par(mfrow = c(1,1))
#' old.par.mar <- par(mar = c(0, 0, 0, 0))
#' par(mfrow = c(11,1))
#' par(mar = c(1,2,1,2))
#' for(p in 1:ncol(shrink1$draws$beta)){plot(density(shrink1$draws$beta[,p]),
#'   xlim=c(-10,10),main=colnames(shrink1$draws$beta)[p])}
#' par(mfrow = old.par.mfrow)
#' par(mar = old.par.mar)
#'
#' # Bayesian horseshoe where first three and last three beta's have different
#' # global shrinkage parameter than other beta's
#' shrink2 <- shrinkem(estimates, covmatrix, type="horseshoe",
#'    group=c(rep(1,3),rep(2,5),rep(1,3)))
#' # posterior modes of middle five estimates are virtually zero
#'
#' # plot posterior densities
#' par(mfrow = c(11,1))
#' par(mar = c(1,2,1,2))
#' for(p in 1:ncol(shrink2$draws$beta)){plot(density(shrink2$draws$beta[,p]),xlim=c(-10,10),
#'   main=colnames(shrink2$draws$beta)[p])}
#' par(mfrow = old.par.mfrow)
#' par(mar = old.par.mar)
#' }
#'
#' @export
shrinkem <- function(x, Sigma, type, group,
                     iterations, burnin, store,
                     cred.level,
                     df1, df2, scale2,
                     lambda2.fixed,
                     lambda2,
                     nugget,
                     ...) {
  UseMethod("shrinkem", x)
}


#' @method shrinkem default
#' @export
shrinkem.default <- function(x, Sigma, type="horseshoe", group=1,
                     iterations = 5e4, burnin = 1e3, store = 1,
                     cred.level = .95,
                     df1=1, df2=1, scale2=1e3,
                     lambda2.fixed = FALSE,
                     lambda2 = NA,
                     nugget = 1e-8,
                     ...){

  if(is.na(lambda2.fixed)){
    stop("'lambda2.fixed' must be 'TRUE' of 'FALSE'.")
  }
  if(is.null(lambda2.fixed)){
    stop("'lambda2.fixed' must be 'TRUE' of 'FALSE'.")
  }
  if(!is.logical(lambda2.fixed)){
    stop("'lambda2.fixed' must be 'TRUE' of 'FALSE'.")
  }
  if(lambda2.fixed){

    if(is.na(lambda2[1])){
      stop("If 'lambda2.fixed = TRUE', then 'lambda2' needs to be a positive scalar.")
    }
    if(is.null(lambda2[1])){
      stop("If 'lambda2.fixed = TRUE', then 'lambda2' needs to be a positive scalar.")
    }
    if(max(lambda2<=0)==1){
      stop("If 'lambda2.fixed = TRUE', then 'lambda2' needs to be positive and of equal length as the number of groups.")
    }
    if(length(lambda2)!=length(unique(group))){
      stop("If 'lambda2.fixed = TRUE', then 'lambda2' needs to be positive and of equal length as the number of groups.")
    }
  }

  if(type=="lasso"){
    Gibbsresults <- normal.lasso(estimate=x, covmatrix=Sigma, group=group,
                                 iterations = iterations, burnin = burnin, store = store,
                                 a1 = df1/2, a2 = df2/2, b1 = scale2,
                                 lambda2.fixed = lambda2.fixed,
                                 lambda2.input = lambda2,
                                 nugget = nugget)
  }
  if(type=="ridge"){
    Gibbsresults <- normal.ridge(estimate=x, covmatrix=Sigma, group=group,
                                 iterations = iterations, burnin = burnin, store = store,
                                 a1 = df1/2, a2 = df2/2, b1 = scale2,
                                 lambda2.fixed = lambda2.fixed,
                                 lambda2.input = lambda2,
                                 nugget = nugget)
  }
  if(type=="horseshoe"){
    Gibbsresults <- normal.horseshoe(estimate=x, covmatrix=Sigma, group=group,
                                     iterations = iterations, burnin = burnin, store = store,
                                     a1 = df1/2, a2 = df2/2, b1 = scale2,
                                     a3 = .5, a4 = 0.5, b2 = 1,
                                     lambda2.fixed = lambda2.fixed,
                                     lambda2.input = lambda2,
                                     nugget = nugget)
  }
  if(sum(type==c("lasso","ridge","horseshoe"))==0){
    stop("argument 'type' must be 'lasso', 'ridge', or 'horseshoe'.")
  }

  if(is.null(names(x))){
    names(x) <- paste0("beta",1:length(x))
  }else{
    colnames(Gibbsresults$beta) <- names(x)
  }

  # handle output
  estimates <- data.frame(
      input.est=x,
      shrunk.mean=apply(Gibbsresults$beta,2,mean),
      shrunk.median=apply(Gibbsresults$beta,2,stats::median),
      shrunk.mode=apply(Gibbsresults$beta,2,get.mode),
      shrunk.lower=apply(Gibbsresults$beta,2,quantile,(1-cred.level)/2),
      shrunk.upper=apply(Gibbsresults$beta,2,quantile,1-(1-cred.level)/2),
      nonzero=!((0 < apply(Gibbsresults$beta,2,quantile,1-(1-cred.level)/2)) &
               (apply(Gibbsresults$beta,2,quantile,(1-cred.level)/2) < 0))
  )

  shrinkem_out <- list(
    estimates=estimates,
    draws=Gibbsresults,
    dim.est=length(x),
    input.est=x,
    call=match.call()
  )

  class(shrinkem_out) <- "shrinkem"

  return(shrinkem_out)
}


# ---------------------------------------------------------------------------
# Gibbs samplers (C++ backends in src/shrinkem_samplers.cpp).
# The R functions below are thin wrappers: they handle group bookkeeping and
# re-attach the column names that downstream code expects, then call the
# compiled samplers.
# ---------------------------------------------------------------------------

#' @useDynLib shrinkem, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats density quantile rgamma
#' @importFrom extraDistr rinvgamma
NULL

# Gibbs sampler for Bayesian (group) horseshoe
normal.horseshoe <- function(estimate, covmatrix, group = 1,
                             iterations = 1e4, burnin = 1e3, store = 10,
                             a1 = .5, a2 = .5, b1 = 1,
                             a3 = .5, a4 = .5, b2 = 1,
                             lambda2.fixed = FALSE, lambda2.input = NA, nugget){

  P <- length(estimate)
  if(is.null(names(estimate))) names(estimate) <- paste0("beta", 1:P)
  if(is.na(group[1]) || length(group) == 1) group <- rep(1, P)
  if(length(group) != P) stop("length of 'group' differs from length of 'estimate'")

  group_idx <- match(group, unique(group)) - 1L
  numGroup  <- length(unique(group))
  li <- if(isTRUE(lambda2.fixed)) lambda2.input else rep(1, numGroup)

  out <- normal_horseshoe_cpp(estimate, covmatrix, as.integer(group_idx), numGroup,
                              iterations, burnin, store, a1, a2, b1, a3, a4, b2,
                              isTRUE(lambda2.fixed), as.numeric(li), nugget)

  colnames(out$beta)    <- names(estimate)
  colnames(out$tau2)    <- paste0("tau2_", 1:P)
  colnames(out$psi2)    <- paste0("psi2_", 1:P)
  colnames(out$lambda2) <- paste0("lambda2_", 1:numGroup)
  colnames(out$gamma2)  <- paste0("gamma2_", 1:numGroup)
  out
}

# Gibbs sampler for Bayesian (group) lasso (Laplace prior)
normal.lasso <- function(estimate, covmatrix, group = 1,
                         iterations = 1e4, burnin = 1e3, store = 10,
                         a1 = .5, a2 = .5, b1 = 1,
                         lambda2.fixed = FALSE, lambda2.input = NA, nugget){

  P <- length(estimate)
  if(is.null(names(estimate))) names(estimate) <- paste0("beta", 1:P)
  if(is.na(group[1]) || length(group) == 1) group <- rep(1, P)
  if(length(group) != P) stop("length of 'group' differs from length of 'estimate'")

  group_idx <- match(group, unique(group)) - 1L
  numGroup  <- length(unique(group))
  li <- if(isTRUE(lambda2.fixed)) lambda2.input else rep(1, numGroup)

  out <- normal_lasso_cpp(estimate, covmatrix, as.integer(group_idx), numGroup,
                          iterations, burnin, store, a1, a2, b1,
                          isTRUE(lambda2.fixed), as.numeric(li), nugget)

  colnames(out$beta)    <- names(estimate)
  colnames(out$tau2)    <- paste0("tau2_", 1:P)
  colnames(out$lambda2) <- paste0("lambda2_", 1:numGroup)
  colnames(out$gamma2)  <- paste0("gamma2_", 1:numGroup)
  out
}

# Gibbs sampler for Bayesian (group) ridge (normal prior)
normal.ridge <- function(estimate, covmatrix, group = 1,
                         iterations = 1e4, burnin = 1e3, store = 10,
                         a1 = .5, a2 = .5, b1 = 1,
                         lambda2.fixed = FALSE, lambda2.input = NA, nugget){

  P <- length(estimate)
  if(is.null(names(estimate))) names(estimate) <- paste0("beta", 1:P)
  if(is.na(group[1]) || length(group) == 1) group <- rep(1, P)
  if(length(group) != P) stop("length of 'group' differs from length of 'estimate'")

  group_idx <- match(group, unique(group)) - 1L
  numGroup  <- length(unique(group))
  li <- if(isTRUE(lambda2.fixed)) lambda2.input else rep(1, numGroup)

  out <- normal_ridge_cpp(estimate, covmatrix, as.integer(group_idx), numGroup,
                          iterations, burnin, store, a1, a2, b1,
                          isTRUE(lambda2.fixed), as.numeric(li), nugget)

  colnames(out$beta) <- names(estimate)
  if(!is.null(out$lambda2)) colnames(out$lambda2) <- paste0("lambda2_", 1:numGroup)
  if(!is.null(out$gamma2)) colnames(out$gamma2)  <- paste0("gamma2_", 1:numGroup)
  out
}



#' @title The (scaled) F Distribution
#' @description Density and random generation for the F distribution with first degrees of freedom \code{df1},
#' second degrees of freedom \code{df2}, and scale parameter \code{beta}.
#' @param x vector of quantities.
#' @param df1 First degrees of freedom
#' @param df2 Second degrees of freedom
#' @param beta Scale parameter
#' @param n number of draws
#' @param log logical; if TRUE, density is given as log(p).
#' @return \code{dF} gives the probability density of the F distribution. \code{rF} gives random draws from the F distribution.
#' @references Mulder and Pericchi (2018). The Matrix-F Prior for Estimating and Testing Covariance Matrices.
#' Bayesian Analysis, 13(4), 1193-1214. <https://doi.org/10.1214/17-BA1092>
#' @importFrom stats rgamma
#' @importFrom extraDistr rinvgamma
#' @name F
#' @examples
#'
#' draws_F <- rF(n=1e4, df1=2, df2=4, beta=1)
#' hist(draws_F,500,xlim=c(0,10),freq=FALSE)
#' seqx <- seq(0,10,length=1e5)
#' lines(seqx,dF(seqx, df1=2, df2=4, beta=1),col=2,lwd=2)
#'

#' @rdname F
#' @export
dF <- function(x, df1, df2, beta, log = FALSE){
  dF_log <- lgamma((df1 + df2)/2) - lgamma(df1/2) - lgamma(df2/2) - df1/2 * log(beta) + (df1/2 - 1) * log(x) -
    (df1 + df2)/2 * log(1 + x/beta)
  if(log == TRUE){
    dF_log
  }else{
    exp(dF_log)
  }
}

#' @rdname F
#' @export
rF <- function(n, df1, df2, beta){
  psi2 <- stats::rgamma(n,shape=df1/2,rate=1/beta)
  extraDistr::rinvgamma(n,alpha=df2/2,beta=psi2)
}

#' @title The matrix F Distribution
#' @description Density and random generation for the matrix variate F distribution with first degrees
#' of freedom \code{df1}, second degrees of freedom \code{df2}, and scale matrix \code{B}.
#' @param x Positive definite matrix of quantities.
#' @param df1 First degrees of freedom
#' @param df2 Second degrees of freedom
#' @param B Positive definite scale matrix
#' @param n Number of draws
#' @param log logical; if TRUE, density is given as log(p).
#' @return \code{dmvF} returns the probability density of the matrix F distribution.
#' \code{rmvF} returns a numeric array, say \code{R}, of dimension  \eqn{p \times p \times n}, where each element
#' \code{R[,,i]} is a positive definite matrix, a realization of the matrix F distribution.
#' @references Mulder and Pericchi (2018). The Matrix-F Prior for Estimating and Testing Covariance Matrices.
#' Bayesian Analysis, 13(4), 1193-1214. <https://doi.org/10.1214/17-BA1092>
#' @importFrom CholWishart lmvgamma
#' @importFrom matrixcalc is.positive.definite
#' @name mvF
#' @examples
#'
#' set.seed(20180222)
#' draws_F <- rmvF(n=1, df1=2, df2=4, B=diag(2))
#' dmvF(draws_F[,,1], df1=2, df2=4, B=diag(2))

#' @rdname mvF
#' @export
dmvF <- function(x,df1,df2,B,log=FALSE){
  if(!matrixcalc::is.positive.definite(B)){
    stop("B must be a positive definite square matrix")
  }
  k <- nrow(B)
  if(df1 <= k-1){
    stop("df1 must be larger than nrow(B)-1")
  }
  if(df2 <= 0){
    stop("df2 must be larger than 0")
  }
  dmvF_log <- CholWishart::lmvgamma((df1+df2+k-1)/2, k) - CholWishart::lmvgamma(df1/2, k) - CholWishart::lmvgamma((df2+k-1)/2, k) -
    df1/2 * c(determinant(diag(2)+1,logarithm=TRUE)$modulus) + (df1-k-1)/2*log(det(x)) -
    (df1+df2+k-1)/2*c(determinant(diag(k)+x%*%solve(B),logarithm=TRUE)$modulus)
  if(log == TRUE){
    dmvF_log
  }else{
    exp(dmvF_log)
  }
}

#' @rdname mvF
#' @export
rmvF <- function(n,df1,df2,B){
  if(!is.matrix(B)){
    stop("B must be a square matrix")
  }
  k <- nrow(B)
  if(df1 <= k-1){
    stop("df1 must be larger than nrow(B)-1")
  }
  if(df2 <= 0){
    stop("df2 must be larger than 0")
  }

  Psi <- stats::rWishart(n,df=df1,Sigma=B)
  array(unlist(lapply(1:n,function(i){
    solve(stats::rWishart(1,df=df2+k-1,Sigma=solve(Psi[,,i]))[,,1])
  })),dim=c(k,k,n))
}

get.mode <- function(draws){
  draws <- draws[is.finite(draws)]
  bw <- tryCatch(stats::bw.SJ(draws), error = function(e) stats::bw.nrd0(draws))
  dens <- stats::density(draws, bw = bw, n = 2048)
  dens$x[which.max(dens$y)]
}

#' @importFrom logspline logspline qlogspline dlogspline
get.mode.logspline <- function(x){
  fit <- logspline(x,error.action=1)
  u1 <- qlogspline(0.01, fit)
  u2 <- qlogspline(0.99, fit)
  u3 <- 1.1 * u1 - 0.1 * u2
  u4 <- 1.1 * u2 - 0.1 * u1
  xx <- (0:(100 - 1))/(100 - 1) * (u4 - u3) + u3
  yy <- dlogspline(xx, fit)
  xx[which(yy==max(yy))]
}
