# Small internal helper functions

#' Clamp values to limits
#'
#' @param x Numeric vector to clamp
#' @param upper,lower Limits to clamp to
clamp <- function(x, lower=-Inf, upper=Inf) {
  x[x > upper] <- upper
  x[x < lower] <- lower
  return(x)
}
