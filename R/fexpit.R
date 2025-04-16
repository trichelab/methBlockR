#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @export 
fexpit <- function(x, sqz=1e-6) {
  unsqueeze(exp(x)/(1+exp(x)), sqz=sqz)
}
