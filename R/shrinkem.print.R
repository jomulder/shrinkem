

#' @method print shrinkem
#' @export
print.shrinkem <- function(x,
                      ...){

  cat("Call:")
  cat("\n")
  print(x$call)

  cat("\n")

  digits <- 3

  print(round(x$estimates,digits))

}


#' @method summary shrinkem
#' @export
summary.shrinkem <- function(object,
                     ...){

  cat("Call:")
  cat("\n")
  print(object$call)

  cat("\n")

  digits <- 3

  print(round(object$estimates,digits))

}


