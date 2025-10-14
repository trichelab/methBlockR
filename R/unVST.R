#' Inverse variance stabilizing transform (sin(x)**2) for a proportion
#' 
#' @param x     a vector of transformed proportions (as from VST(p)) 
#'
#' @return      values between 0 and 1
#'
#' @details     Using asin(sqrt(p)) for VST avoids negatives and infinities.
#'
#' @examples
#' 
#' dat <- rbeta(n=100, 1, 7)
#' vsted <- VST(dat)
#' all(abs(unVST(vsted) - dat) < .Machine$double.eps)
#' 
#' @seealso VST
#' @seealso extremality
#'
#' @export 
#'
#' @export 
#'
unVST <- function(x) return(sin(x)**2)
