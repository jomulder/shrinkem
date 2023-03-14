#' @importFrom bain get_estimates

#' @method shrinkem lm
#' @export
shrinkem.lm <- function(x, Sigma, type="horseshoe", group=1,
                        iterations = 5e4, burnin = 1e3, store = 10,
                        cred.level = .95,
                        df1=1, df2=1, scale2=1e3, ...){

  #Extract summary statistics
  Args <- as.list(match.call()[-1])
  get_est <- bain::get_estimates(x)
  Args$x <- get_est$estimate
  Args$Sigma <- get_est$Sigma[[1]]
  Args$type <- type
  Args$group <- group
  Args$iterations <- iterations
  Args$burnin <- burnin
  Args$store <- store
  Args$cred.level <- cred.level
  Args$df1 <- df1
  Args$df2 <- df2
  Args$scale2 <- scale2
  out <- do.call(shrinkem, Args)
  out$model <- x
  out$call <- match.call()
  out

}
