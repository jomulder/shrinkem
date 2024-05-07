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
#' distribution. The default is \code{df2 = 1e3}.
#' @param iterations Number of posterior draws after burnin. Default = 5e4.
#' @param burnin Number of posterior draws in burnin. Default = 1e3.
#' @param store Store every store-th draw from posterior. Default = 10.
#' @param ... Parameters passed to and from other functions.
#' @rdname shrinkem
#' @references Karimovo, D., Leenders, R., van Erp, S., and Mulder, J. (in prep.).
#' @examples
#' \donttest{
#' # EXAMPLE
#' estimates <- -5:5
#' covmatrix <- diag(11)
#' # Bayesian horseshoe where all beta's have the same global shrinkage (using default 'group' argument)
#' shrink1 <- shrinkem(estimates, covmatrix, type="horseshoe")
#' # posterior modes of middle three estimates are practically zero
#' print(shrink1)
#' # how traceplots
#' bayesplot::mcmc_trace(shrink1$draws$beta)
#' # plot posterior densities
#' bayesplot::mcmc_areas_ridges(shrink1$draws$beta)
#' # Bayesian horseshoe where first three and last three beta's have different global shrinkage parameter
#' # than other beta's
#' shrink2 <- shrinkem(estimates, covmatrix, type="horseshoe", group=c(rep(1,3),rep(2,5),rep(1,3)))
#' # posterior modes of middle five estimates are practically zero
#' print(shrink2)
#' # show traceplots
#' bayesplot::mcmc_trace(shrink2$draws$beta)
#' # plot posterior densities
#' bayesplot::mcmc_areas_ridges(shrink2$draws$beta)
#' }
#' @export
shrinkem <- function(x, Sigma, type, group,
                     iterations, burnin, store,
                     cred.level,
                     df1, df2, scale2, ...) {
  UseMethod("shrinkem", x)
}


#' @method shrinkem default
#' @export
shrinkem.default <- function(x, Sigma, type="horseshoe", group=1,
                     iterations = 5e4, burnin = 1e3, store = 10,
                     cred.level = .95,
                     df1=1, df2=1, scale2=1e3, ...){

  if(type=="lasso"){
    Gibbsresults <- normal.lasso(estimate=x, covmatrix=Sigma, group=group,
                                 iterations = iterations, burnin = burnin, store = store,
                                 a1 = df1/2, a2 = df2/2, b1 = scale2)
  }
  if(type=="ridge"){
    Gibbsresults <- normal.ridge(estimate=x, covmatrix=Sigma, group=group,
                                     iterations = iterations, burnin = burnin, store = store,
                                     a1 = df1/2, a2 = df2/2, b1 = scale2)
  }
  if(type=="horseshoe"){
    Gibbsresults <- normal.horseshoe(estimate=x, covmatrix=Sigma, group=group,
                                 iterations = iterations, burnin = burnin, store = store,
                                 a1 = df1/2, a2 = df2/2, b1 = scale2,
                                 a3 = .5, a4 = 0.5, b2 = 1)
  }
  if(is.null(names(x))){
    names(x) <- paste0("beta",1:length(x))
  }else{
    colnames(Gibbsresults$beta) <- names(x)
  }

  # handle output
  estimates <- data.frame(
      est=x,
      shrunk.mean=apply(Gibbsresults$beta,2,mean),
      shrunk.median=apply(Gibbsresults$beta,2,mean),
      shrunk.mode=apply(Gibbsresults$beta,2,get.mode),
      shrunk.lower=apply(Gibbsresults$beta,2,quantile,(1-cred.level)/2),
      shrunk.upper=apply(Gibbsresults$beta,2,quantile,1-(1-cred.level)/2),
      nonzero=!((0 < apply(Gibbsresults$beta,2,quantile,1-(1-cred.level)/2)) &
               (apply(Gibbsresults$beta,2,quantile,(1-cred.level)/2) < 0))
  )

  shrinkem_out <- list(
    estimates=estimates,
    draws=Gibbsresults,
    dim.beta=length(x),
    model=x,
    call=match.call()
  )

  class(shrinkem_out) <- "shrinkem"

  return(shrinkem_out)
}

# Gibbs sampler for Bayesian (group) horseshoe
#' @importFrom stats density
#' @importFrom bayesplot mcmc_areas_ridges mcmc_trace
#' @importFrom stats quantile
#' @importFrom mvtnorm rmvnorm
#' @importFrom extraDistr rinvgamma
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats rgamma
#' @importFrom brms rinv_gaussian
normal.horseshoe <- function(estimate, covmatrix, group=1,
                                   iterations = 1e4, burnin = 1e3, store = 10,
                                   a1 = .5, a2 = 0.5, b1 = 1,
                                   a3 = .5, a4 = 0.5, b2 = 1){

  # Define dimensions
  P <- length(estimate) # number of covariates without intercept
  if(is.null(names(estimate))){
    names(estimate) <- paste0("beta",1:P)
  }
  if(is.na(group[1])){
    group <- rep(1,P)
  }
  if(length(group)==1){
    group <- rep(1,P)
  }
  if(length(group)!=P){
    stop("length of 'group' differs from length of 'estimate'")
  }
  groupIndices <- unique(group)
  numGroup <- length(groupIndices)
  whichGroup <- lapply(1:numGroup,function(gr){
    groupIndices[gr] == group
  })
  whichGroupMat <- do.call(rbind,whichGroup)
  lenthGroup <- unlist(lapply(whichGroup,sum))

  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = iterations/store, ncol = P)
  colnames(beta_STORE) <- names(estimate)
  tau2_STORE <- psi2_STORE <- matrix(0,nrow = iterations/store, ncol = P)
  lambda2_STORE <- gamma2_STORE <- matrix(0,nrow = iterations/store, ncol = numGroup)
  colnames(lambda2_STORE) <- paste0("lambda2_",1:numGroup)
  colnames(gamma2_STORE) <- paste0("gamma2_",1:numGroup)
  colnames(tau2_STORE) <- paste0("tau2_",1:P)
  colnames(psi2_STORE) <- paste0("psi2_",1:P)

  # initial values
  lambda2 <- gamma2 <- rep(1,numGroup)
  tau2 <- psi2 <- rep(1, P)
  lambda2vec <- c(lambda2 %*% whichGroupMat)
  D_inv <- diag(1 / (tau2 * lambda2vec))
  beta <- estimate

  covmatrixInv <- solve(covmatrix)

  #nugget for avoiding singular covariance matrix
  nugget <- 1e-8

  print("Start burnin ... ")
  # Sampler iterations ----
  pb = txtProgressBar(min = 0, max = burnin, initial = 0)
  for (t in 1:burnin){

    # Sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # Sample tau2
    tau2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*lambda2vec))
    tau2 <- tau2 + nugget

    # Sample psi2 (mixing parameter of tau2)
    psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 + 1/tau2)

    # Sample lambda2 & gamma2
    for(gr in 1:numGroup){
      lambda2[gr] <- extraDistr::rinvgamma(1, a1 + lenthGroup[gr]/2, gamma2[gr] + sum(beta[whichGroup[[gr]]]^2/(2*tau2[whichGroup[[gr]]])) )
      gamma2[gr] <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2[gr])
    }
    lambda2vec <- c(lambda2 %*% whichGroupMat)

    # update conditional prior covariance matrix
    D_inv <- diag(1/(tau2 * lambda2vec),nrow=P)

    setTxtProgressBar(pb,t)

  }

  cat("\n")
  print("Start posterior sampling ... ")
  # Sampler iterations ----
  pb = txtProgressBar(min = 0, max = iterations, initial = 0)
  storecount <- 0
  for (t in 1:iterations){

    # Sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # Sample tau2
    tau2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*lambda2vec))
    tau2 <- tau2 + nugget

    # Sample psi2 (mixing parameter of tau2)
    psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 + 1/tau2)

    # Sample lambda2 & gamma2
    for(gr in 1:numGroup){
      lambda2[gr] <- extraDistr::rinvgamma(1, a1 + lenthGroup[gr]/2, gamma2[gr] + sum(beta[whichGroup[[gr]]]^2/(2*tau2[whichGroup[[gr]]])) )
      gamma2[gr] <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2[gr])
    }
    lambda2vec <- c(lambda2 %*% whichGroupMat)

    # update conditional prior covariance matrix
    D_inv <- diag(1/(tau2 * lambda2vec),nrow=P)

    #print(c(tau2,lambda2))

    # save each 'store' iterations
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      tau2_STORE[storecount,] <- tau2
      lambda2_STORE[storecount,] <- lambda2
      gamma2_STORE[storecount,] <- gamma2
      psi2_STORE[storecount,] <- psi2

    }
    ###################

    setTxtProgressBar(pb,t)
  }

  return(list(beta = beta_STORE, tau2 = tau2_STORE, gamma2 = gamma2_STORE,
              psi2 = psi2_STORE, lambda2 = lambda2_STORE))
}

# Gibbs sampler for Bayesian (group) lasso (LaPlace prior)
normal.lasso <- function(estimate, covmatrix, group=1,
                               iterations = 1e4, burnin = 1e3, store = 10,
                               a1 = .5, a2 = 0.5, b1 = 1){

  P <- length(estimate) # number of covariates without intercept
  if(is.null(names(estimate))){
    names(estimate) <- paste0("beta",1:P)
  }
  if(is.na(group[1])){
    group <- rep(1,P)
  }
  if(length(group)==1){
    group <- rep(1,P)
  }
  if(length(group)!=P){
    stop("length of 'group' differs from length of 'estimate'")
  }
  groupIndices <- unique(group)
  numGroup <- length(groupIndices)
  whichGroup <- lapply(1:numGroup,function(gr){
    groupIndices[gr] == group
  })
  whichGroupMat <- do.call(rbind,whichGroup)
  lenthGroup <- unlist(lapply(whichGroup,sum))

  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = iterations/store, ncol = P)
  colnames(beta_STORE) <- names(estimate)
  lambda2_STORE <- gamma2_STORE <- matrix(0,nrow = iterations/store, ncol = numGroup)
  colnames(lambda2_STORE) <- paste0("lambda2_",1:numGroup)
  colnames(gamma2_STORE) <- paste0("gamma2_",1:numGroup)
  tau2_STORE <- matrix(0,nrow = iterations/store, ncol = P)
  colnames(tau2_STORE) <- paste0("tau2_",1:P)

  # initial values
  tau2 <- rep(1, P)
  lambda2 <- gamma2 <- rep(1,numGroup)
  lambda2vec <- c(lambda2 %*% whichGroupMat)
  D_inv <- diag(1/(tau2*lambda2vec))
  beta <- estimate

  covmatrixInv <- solve(covmatrix)

  print("Start burnin ... ")
  # Sampler iterations ----
  pb = txtProgressBar(min = 0, max = burnin, initial = 0)
  for (t in 1:burnin){
    # 1. sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 3. Sample 1/tau^2
    mu_tau <- sqrt(lambda2vec/beta^2)
    tau2 <- 1/brms::rinv_gaussian(P, mu = mu_tau, shape = 1)

    # 4. Sample lambda2 & gamma2
    for(gr in 1:numGroup){
      lambda2[gr] <- extraDistr::rinvgamma(1, a1 + lenthGroup[gr]/2, gamma2[gr] + sum(beta[whichGroup[[gr]]]^2/(2*tau2[whichGroup[[gr]]])) )
      gamma2[gr] <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2[gr])
    }
    lambda2vec <- c(lambda2 %*% whichGroupMat)

    # update conditional prior covariance matrix
    D_inv <- diag(1/(tau2 * lambda2vec),nrow=P)

    setTxtProgressBar(pb,t)
  }

  cat("\n")
  print("Start posterior sampling ... ")
  # Sampler iterations ----
  pb = txtProgressBar(min = 0, max = iterations, initial = 0)
  storecount <- 0
  for (t in 1:iterations){
    # 1. sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 3. Sample 1/tau^2
    mu_tau <- sqrt(lambda2vec/beta^2)
    tau2 <- 1/brms::rinv_gaussian(P, mu = mu_tau, shape = 1)

    # 4. Sample lambda2 & gamma2
    for(gr in 1:numGroup){
      lambda2[gr] <- extraDistr::rinvgamma(1, a1 + lenthGroup[gr]/2, gamma2[gr] + sum(beta[whichGroup[[gr]]]^2/(2*tau2[whichGroup[[gr]]])) )
      gamma2[gr] <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2[gr])
    }
    lambda2vec <- c(lambda2 %*% whichGroupMat)

    # update conditional prior covariance matrix
    D_inv <- diag(1/(tau2 * lambda2vec),nrow=P)

    # save each 'store' iterations
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      lambda2_STORE[storecount,] <- lambda2
      tau2_STORE[storecount,] <- tau2
      gamma2_STORE[storecount,] <- gamma2

    }
    ###################

    setTxtProgressBar(pb,t)
  }

  colnames(beta_STORE) <- names(estimate)

  return(list(beta = beta_STORE,
              lambda2 = lambda2_STORE,
              gamma2 = gamma2_STORE,
              tau2 = tau2_STORE))
}

# Gibbs sampler for Bayesian (group) ridge (normal prior)
normal.ridge <- function(estimate, covmatrix, group = 1,
                         iterations = 1e4, burnin = 1e3, store = 10,
                         a1 = .5, a2 = 0.5, b1 = 1){

  P <- length(estimate) # number of covariates without intercept
  if(is.null(names(estimate))){
    names(estimate) <- paste0("beta",1:P)
  }
  if(is.na(group[1])){
    group <- rep(1,P)
  }
  if(length(group)==1){
    group <- rep(1,P)
  }
  if(length(group)!=P){
    stop("length of 'group' differs from length of 'estimate'")
  }
  groupIndices <- unique(group)
  numGroup <- length(groupIndices)
  whichGroup <- lapply(1:numGroup,function(gr){
    groupIndices[gr] == group
  })
  whichGroupMat <- do.call(rbind,whichGroup)
  lenthGroup <- unlist(lapply(whichGroup,sum))

  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = iterations/store, ncol = P)
  lambda2_STORE <- gamma2_STORE <- matrix(0,nrow = iterations/store, ncol = numGroup)
  colnames(lambda2_STORE) <- paste0("lambda2_",1:numGroup)
  colnames(gamma2_STORE) <- paste0("gamma2_",1:numGroup)

  # initial values
  lambda2 <- gamma2 <- rep(1,numGroup)
  lambda2vec <- c(lambda2 %*% whichGroupMat)
  D_inv <- diag(1/lambda2vec)
  beta <- estimate

  covmatrixInv <- solve(covmatrix)

  print("Start burnin ... ")
  # Sampler iterations ----
  pb = txtProgressBar(min = 0, max = burnin, initial = 0)
  for (t in 1:burnin){

    # sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(mvtnorm::rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 4. Sample lambda2 & gamma2
    for(gr in 1:numGroup){
      lambda2[gr] <- extraDistr::rinvgamma(1, a1 + lenthGroup[gr]/2, gamma2[gr] + sum(beta[whichGroup[[gr]]]^2/2) )
      gamma2[gr] <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2[gr])
    }
    lambda2vec <- c(lambda2 %*% whichGroupMat)

    # update conditional prior covariance matrix
    D_inv <- diag(1/lambda2vec)

    setTxtProgressBar(pb,t)
  }

  cat("\n")
  print("Start posterior sampling ... ")
  # Sampler iterations ----
  pb = txtProgressBar(min = 0, max = iterations, initial = 0)
  storecount <- 0
  for (t in 1:iterations){

    # sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 4. Sample lambda2 & gamma2
    for(gr in 1:numGroup){
      lambda2[gr] <- extraDistr::rinvgamma(1, a1 + lenthGroup[gr]/2, gamma2[gr] + sum(beta[whichGroup[[gr]]]^2/2) )
      gamma2[gr] <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2[gr])
    }
    lambda2vec <- c(lambda2 %*% whichGroupMat)

    # update conditional prior covariance matrix
    D_inv <- diag(1/lambda2vec)

    # save each 'store' iterations
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      lambda2_STORE[storecount,] <- lambda2
      gamma2_STORE[storecount,] <- gamma2
    }
    ###################

    setTxtProgressBar(pb,t)
  }

  colnames(beta_STORE) <- names(estimate)

  return(list(beta = beta_STORE,
              lambda2 = lambda2_STORE,
              gamma2 = gamma2_STORE))
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
  dens <- density(draws)
  dens$x[which(max(dens$y)==dens$y)[1]]
}

