

The `R` package **shrinkem** implements several Bayesian regularization algorithms
(e.g., adaptive ridge, adaptivelasso, adaptive horseshoe) and assuming a Gaussian approximation of the estimates and (error) covariance matrix. 

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
shrink1 <- shrinkem(estimates, covmatrix, type="horseshoe", group=rep(1,11))
# posterior modes of middle three estimates are practically zero
print(shrink1)
# how traceplots
bayesplot::mcmc_trace(shrink1$draws$beta)
# plot posterior densities
bayesplot::mcmc_areas_ridges(shrink1$draws$beta)
# Bayesian horseshoe where first three and last three beta's have different global shrinkage
# parameter than other beta's
shrink2 <- shrinkem(estimates, covmatrix, type="horseshoe", group=c(rep(1,3),rep(2,5),rep(1,3)))
# posterior modes of middle five estimates are practically zero
print(shrink2)
# show traceplots
bayesplot::mcmc_trace(shrink2$draws$beta)
# plot posterior densities
bayesplot::mcmc_areas_ridges(shrink2$draws$beta)
```


## Citing **shrinkem**

You can cite the package and the paper using the following reference

> Karimova, D., Leenders, R., van Erp, S., and Mulder, J. (in prep.)

## Contributing and Contact Information

If you have ideas, please get involved. You can contribute by opening an
issue on GitHub, or sending a pull request with proposed features.

  - File a GitHub issue [here](https://github.com/jomulder/shrinkem)
  - Make a pull request [here](https://github.com/jomulder/shrinkem/pulls)

By participating in this project, you agree to abide by the Contributor
Code of Conduct v2.0.
