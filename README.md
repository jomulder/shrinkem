

The `R` package **shrinkem** allows approximate Bayesian regularization
(e.g., ridge, lasso, horseshoe) using Gaussian approximations of the errors of the estimates of the key parameters.

## Installation

The developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("jomulder/shrinkem")
```

## Example analyses

``` r
estimates <- -5:5
covmatrix <- diag(11)
# Bayesian horseshoe where all beta's have the same global shrinkage (using default 'group' argument)
shrink1 <- shrinkem(estimates, covmatrix, type="horseshoe")
# posterior modes of middle three estimates are practically zero
print(shrink1)
# how traceplots
par(mfrow=c(3,4))
par(mar=c(2,2,2,2))
for(p in 1:ncol(shrink1$draws$beta)){plot(shrink1$draws$beta[,p],type="l",
                                          main=colnames(shrink1$draws$beta)[p])}
# plot posterior densities
par(mfrow=c(11,1))
par(mar=c(1,2,1,2))
for(p in 1:ncol(shrink1$draws$beta)){plot(density(shrink1$draws$beta[,p]),xlim=c(-10,10),
                                          main=colnames(shrink1$draws$beta)[p])}
# Bayesian horseshoe where first three and last three beta's have different global shrinkage
# parameter than other beta's (using the 'group' argument)
shrink2 <- shrinkem(estimates, covmatrix, type="horseshoe", group=c(rep(1,3),rep(2,5),rep(1,3)))
# posterior modes of middle five estimates are practically zero
print(shrink2)
# show traceplots
par(mfrow=c(3,4))
par(mar=c(2,2,2,2))
for(p in 1:ncol(shrink2$draws$beta)){plot(shrink2$draws$beta[,p],type="l",
                                          main=colnames(shrink2$draws$beta)[p])}
# plot posterior densities
par(mfrow=c(11,1))
par(mar=c(1,2,1,2))
for(p in 1:ncol(shrink2$draws$beta)){plot(density(shrink2$draws$beta[,p]),xlim=c(-10,10),
                                          main=colnames(shrink2$draws$beta)[p])}
```


## Citing **shrinkem**

You can cite the package and the paper using the following reference

> Karimova, D., Leenders, R., van Erp, S., and Mulder, J. (preprint). Honey, I shunk the irrelevant
effects! Simple and Fast Approximate Bayesian Regularization. 

## Contributing and Contact Information

If you have suggestions, please get involved. You can contribute by opening an
issue on GitHub, or sending a pull request with proposed features.

  - File a GitHub issue [here](https://github.com/jomulder/shrinkem)
  - Make a pull request [here](https://github.com/jomulder/shrinkem/pulls)

By participating in this project, you agree to abide by the Contributor
Code of Conduct v2.0.
