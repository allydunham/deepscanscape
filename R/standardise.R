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
#' @param base Log base of current scores (used for log transform)
#'
#' @export
transform_er <- function(x, trans, base=2) {
  methods <- c("vamp-Seq" = transform_vamp,
              "log" = transform_log,
              "unit" = transform_unit)

  if (is.character(trans)) {
    ind <- pmatch(stringr::str_to_lower(trans), names(methods))
    if (is.na(ind)) {
      stop(paste0("Unrecognised transform \"", trans, "\"\n",
                  "Recognised options: ", paste(methods, collapse = ", ")))
    }
    trans <- methods[[ind]]
  }

  return(trans(x, base))
}

# Transform data pro
# @keywords internal
transform_vamp <- function(x, ...) {
  y <- 1 + (x - 1) / -min(x - 1, na.rm = TRUE)
  return(log2(y + min(y[y > 0], na.rm = TRUE)))
}

transform_log <- function(x, base=2, ...) {
  stop("Not implemented yet")
}

transform_unit <- function(x, ...) {
  stop("Not implemented yet")
}
