#' Helper function: squeezed logit
#'
#' @param p       a vector of values between 0 and 1 inclusive
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between -Inf and +Inf
#'
#' @export 
flogit <- function(p, sqz=0.000001) {
  p <- squeeze(p, sqz=sqz)
  log(p/(1 - p))
}
