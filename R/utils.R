# Small internal helper functions

#' Clamp values to limits
#'
#' Clamp numeric vectors to be within specified upper and/or lower bounds.
#'
#' @param x Numeric vector to clamp.
#' @param upper,lower Limits to clamp to.
#' @return Numeric vector with the same values as x, except those outside the limits are changed to be at the limit.
#' @keywords internal
clamp <- function(x, lower=-Inf, upper=Inf) {
  x[x > upper] <- upper
  x[x < lower] <- lower
  return(x)
}
