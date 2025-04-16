#' squeeze a proportion by a finite amount
#' 
#' @param p   a vector of proportions (for squeeze)
#' @param x   a vector of real numbers (for unsqueeze)
#' @param sqz a squeeze factor (1e-6, by default)
#' @param min minimum unsqueezed value possible (0)
#' @param max maximum unsqueezed value possible (1) 
#'
#' @return (0 + sqz) <= p <= (1 - sqz), or for unsqueeze, min <= p <= max.
#'
#' @details squeeze and unsqueeze avoid -Inf, +Inf, and NaNs in transforms.
#'
#' @examples
#'
#' expit <- function(x) exp(x) / (1 + exp(x))
#' logit <- function(p) log(p / (1 - p))
#'
#' set.seed(1234)
#' x <- rnorm(n=1e6)
#' p <- round(expit(x), 2)
#' summary(logit(p) - x)
#'
#' set.seed(12345)
#' x <- rnorm(n=1e6)
#' p <- round(expit(x), 2) 
#' summary(logit(p) - x)          # oopsie!
#' summary(logit(squeeze(p)) - x) # better!
#' 
#' @aliases unsqueeze
#' 
#' @export 
#'
squeeze <- function(p, sqz=0.000001, min=0, max=1) {

  p[ p > max ] <- max 
  p[ p < min ] <- min 
  p <- ((((p - min) / (max - min)) - 0.5) * (1 - sqz)) + 0.5
  return(p)

}


#' @export 
#'
unsqueeze <- function(p, sqz=0.000001, min=0, max=1) {

  p <- ((((p - 0.5) / (1 - sqz)) + 0.5) * (max - min)) + min
  p[ p > max ] <- max 
  p[ p < min ] <- min
  return(p) 

}
