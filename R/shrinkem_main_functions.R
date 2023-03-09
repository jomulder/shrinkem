#' @title Fast Bayesian regularization using Gaussian approximations
#' @rdname shrinkem
#' @export
shrinkem <- function(x, covmatrix, group=1, type="lasso",
                     iterations = 1e4, burnin = 1e4, store = 10,
                     sign.level=.95,
                     df1=1, df2=1, scale2=1, ...) {
  UseMethod("BF", x)
}


#' @method shrinkem default
#' @export
shrinkem.default <- function(x, sigma, group=1, type="lasso",
                     iterations = 1e4, burnin = 1e4, store = 10,
                     df1=1, df2=1, scale2=1){

  if(type=="lasso"){
    Gibbsresults <- normal.lasso(estimate=x, covmatrix=sigma, group=group,
                                 iterations = iterations, burnin = burnin, store = store,
                                 a1 = df1/2, a2 = df2/2, b1 = scale2)
  }
  if(type=="ridge"){
    Gibbsresults <- normal.ridge(estimate=x, covmatrix=sigma, group=group,
                                     iterations = iterations, burnin = burnin, store = store,
                                     a1 = df1/2, a2 = df2/2, b1 = scale2)
  }
  if(type=="horseshoe"){
    Gibbsresults <- normal.horseshoe(estimate=x, covmatrix=sigma, group=group,
                                 iterations = iterations, burnin = burnin, store = store,
                                 a1 = df1/2, a2 = df2/2, b1 = scale2,
                                 a3 = .5, a4 = 0.5, b2 = 1)
  }

  # handle

}

# Gibbs sampler for Bayesian (group) horseshoe
#' @importFrom mvtnorm rmvnorm
#' @importFrom extraDistr rinvgamma
#' @importFrom stats rgamma
#' @importFrom brms rinv_gaussian
normal.horseshoe <- function(estimate, covmatrix, group=1,
                                   iterations = 1e4, burnin = 1e4, store = 10,
                                   a1 = .5, a2 = 0.5, b2 = 1,
                                   a3 = .5, a4 = 0.5, b1 = 1){

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

  print("Start burnin ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
  for (t in 1:burnin){

    # Sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # Sample tau2
    tau2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*lambda2vec))

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

    pb$tick()
    Sys.sleep(1/burnin)
  }

  print("Start posterior sampling ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = iterations, clear = F, width = 80)
  storecount <- 0
  for (t in 1:iterations){

    # Sample beta
    var_beta <- solve(covmatrixInv + D_inv)
    mu_beta <- c(var_beta %*% covmatrixInv %*% estimate)
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # Sample tau2
    tau2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*lambda2vec))

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

    pb$tick()
    Sys.sleep(1/iterations)
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
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
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

    pb$tick()
    Sys.sleep(1/burnin)
  }

  print("Start posterior sampling ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = iterations, clear = F, width = 80)
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

    pb$tick()
    Sys.sleep(1/iterations)
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
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
  for (t in 1:burnin){

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

    pb$tick()
    Sys.sleep(1/burnin)
  }

  print("Start posterior sampling ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = iterations, clear = F, width = 80)
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

    pb$tick()
    Sys.sleep(1/iterations)
  }

  colnames(beta_STORE) <- names(estimate)

  return(list(beta = beta_STORE,
              lambda2 = lambda2_STORE,
              gamma2 = gamma2_STORE))
}



