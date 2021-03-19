# Normalise and standardise DMS scores
# TODO - test functions
# TODO - Add additional transforms

#' Normalise deep mutational scanning fitness scores
#'
#' @param x Vector of fitness scores
#' @param q Quantile to normalise against. Values other than 0.1 differ from
#' those used in the combined dataset
#'
#' @export
normalise_er <- function(x, q=0.1) {
  q <- stats::quantile(x, q, na.rm = TRUE)
  return(x / -stats::median(x[x <= q], na.rm = TRUE))
}

#' Standardise DMS data
#'
#' @param x Vector of fitness scores
#' @param trans Transform to apply (see description). Accepts either a
#' string or function
#'
#' @export
transform_er <- function(x, trans = c("log2", "vamp-seq", "unit")) {
  if (is.character(trans)) {
    trans <- match.arg(trans)
    methods <- c("vamp-Seq" = transform_vamp, "unit" = transform_unit, "log2" = log2)
    f <- methods[[trans]]
  } else if (is.function(trans)) {
    f <- trans
  } else {
    stop("Unrecognised transformation. Pass a supported string or a function")
  }

  return(f(x))
}

#' Transform VAMP-seq data
#' @keywords internal
transform_vamp <- function(x, ...) {
  y <- 1 + (x - 1) / -min(x - 1, na.rm = TRUE)
  return(log2(y + min(y[y > 0], na.rm = TRUE)))
}

#' Transform unit scaled data
#' @keywords internal
transform_unit <- function(x, ...) {
  stop("Not implemented yet")
}
