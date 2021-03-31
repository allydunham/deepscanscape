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

# TODO - Expand docs
#' Cosine distance matrix
#'
#' Calculate the cosine distance between rows of a matrix or matrices.
#'
#' Floating point imprecision in the implementation of \code{\link[base]{tcrossprod}} leads to values slightly outside
#' the domain of acos, so cosine similarities are rounded to 8 digits before being converted to distances. This is well
#' beyond the precision required for most applications and does not pose a problem, especially here where the
#' experimental uncertainty of ER scores are greater than this.
#'
#' @param x,y Numeric matrices. The distance between rows of x is calculated if y = NULL.
#' @return Numeric matrix of the cosine distance between the rows of x and y.
#' @export
cosine_distance_matrix <- function(x, y = NULL) {
  if (is.null(y)) {
    y <- x
  }
  cosine <- tcrossprod(x, y) / sqrt(tcrossprod(rowSums(x^2), rowSums(y^2)))
  cosine <- acos(round(cosine, digits = 8)) / pi
  return(cosine)
}

#' Identify permissive positions
#'
#' Identify rows of a matrix of ER scores where all scores are less than 0.4 in absolute value
#'
#' @param x Numeric matrix of ER scores
#' @return Logical vector corresponding to results on each row
#' @keywords internal
permissive_positions <- function(x) {
  apply(abs(x) < 0.4, 1, all)
}
