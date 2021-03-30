# Check a deep_mutational_scan for common anomalies
# TODO - Test these funcs

#' Check a deep_mutational_scan for common anomalies
#'
#' Identify common potential anomalies in deep mutational scanning data, for example an unusual distribution of ER
#' scores or unexpected cluster assignment. Some checks require the data to be annotated and will be skipped, with a
#' warning, if it is not.
#'
#' The following checks are performed:
#' \itemize{
#'   \item Median ER is positive, suggesting it has possibly been inverted.
#' }
#'
#' If \code{x} is annotated the following additional checks are performed:
#' \itemize{
#'   \item There's an unusually low number of common subtypes (1 and 2), suggesting unusual positions or an error like
#'   an incorrect sequence.
#' }
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @param warn Produce warnings when anomalies are detected.
#' @return Data frame of anomalies detected
#' @export
check_data <- function(x, warn = TRUE) {
  if (!is.deep_mutational_scan(x)) {
    stop("x is not a deep_mutational_scan")
  }

  x <- validate_deep_mutational_scan(x)

  ids <- c()

  # Check median ER score
  if (stats::median(as.matrix(x$data[amino_acids])) > 0) {
    ids <- c(ids, "positive_er")
  }


  if (x$annotated) {
    # Expected proportion of 1/2 subtypes
    common_subtype <- sum(stringr::str_ends(x$data$cluster, "[12]"))
    n <- sum(stringr::str_ends(x$data$cluster, "[^A]"))
    btest <- stats::binom.test(common_subtype, n, p = 0.6, alternative = "less")

    if (btest$p.value < 0.05) {
      ids <- c(ids, "low_x1/2")
    }

  } else {
    warning("x is not annotated, skipping some checks")
  }


  out <- check_msgs[check_msgs$id %in% ids, ]

  if (warn) {
    sapply(stringr::str_wrap(stringr::str_c(out$cause, ":", out$explanation), width = 100), warning)
  }

  return(out)
}
