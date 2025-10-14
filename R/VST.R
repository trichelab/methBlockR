#' Variance stabilizing transform (asin(sqrt(p))) for a proportion
#' 
#' @param p     a vector of proportions (VST) 
#' @param min   minimum valid input (0)
#' @param max   maximum valid input (1)
#'
#' @return      values between 0 and 1.570796327
#'
#' @details     Using asin(sqrt(p)) for VST avoids negatives and infinities.
#'
#' @examples
#' 
#' dat <- rbeta(n=100, 1, 7)
#' vsted <- VST(dat)
#' all(abs(unVST(vsted) - dat) < .Machine$double.eps)
#' 
#' @seealso unVST
#' @seealso extremality
#'
#' @export 
#'
VST <- function(p, min=0, max=1) return(asin(sqrt((p - min) / max)))
