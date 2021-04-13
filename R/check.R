# Check a deep_mutational_scan for common anomalies

#' Check a deep_mutational_scan for common anomalies
#'
#' Identify common potential anomalies in deep mutational scanning data, for example an unusual distribution of ER
#' scores or unexpected cluster assignment. Some checks require the data to be annotated and will be skipped, with a
#' warning, if it is not.
#'
#' The following checks are performed:
#' \itemize{
#'   \item Median ER is positive, suggesting it has possibly been inverted.
#'   \item Mean WT |ER| > 0.1, suggesting an incorrect sequence.
#' }
#'
#' If \code{x} is annotated the following additional checks are performed:
#' \itemize{
#'   \item Positions fall outside the mutational landscape, suggesting they may have unusual function or the experiment
#'   selects unusual properties.
#'   \item There are an unusually high number of permissive positions, based on a binomial test. This can suggest a
#'   weakly selective region or experiment.
#' }
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @param warn Produce warnings when anomalies are detected.
#' @return Data frame of anomalies detected
#' @examples
#' check_data(deep_scans$p53, warn = TRUE)
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

  # Check WT scores are low
  wt <- tidyr::pivot_longer(x$data[c("position", "wt", amino_acids)], cols = dplyr::all_of(amino_acids),
                            names_to = "mut", values_to = "score")
  wt <- dplyr::filter(wt, .data$wt == .data$mut)
  if (abs(mean(wt$score)) > 0.1) {
    ids <- c(ids, "wt_high")
  }


  if (x$annotated) {
    # Positions outside the explored landscape
    if (length(landscape_outliers(x)) > 0) {
      ids <- c(ids, "landscape_outliers")
    }

    # Proportion of permissive
    nperm <- sum(stringr::str_ends(x$data$cluster, "P"))
    btest <- stats::binom.test(nperm, dim(x)[1], p = 0.116, alternative = "greater")
    if (btest$p.value < 0.05) {
      ids <- c(ids, "permissive")
    }

  } else {
    warning("x is not annotated, skipping some checks")
  }


  out <- check_msgs[check_msgs$id %in% ids, ]

  if (warn & nrow(out) > 0) {
    warning(stringr::str_c(out$cause, sep = "\n"), call. = FALSE)
  }

  return(out)
}

#' Identify positions outside the mutational landscape
#'
#' Find positions that lie outside the studied deep mutational landscape.
#'
#' @param x Annotated \code{\link{deep_mutational_scan}}.
#' @param cutoff Minimum distance from any studied position to be labelled an outlier.
#' @return integer vector containing the rows of outlying positions
#' @examples
#' x <- annotate(deep_scans$p53)
#' landscape_outliers(x)
#' @export
landscape_outliers <- function(x, cutoff = 0.25) {
  if (!is.deep_mutational_scan(x)) {
    stop("x is not a deep_mutational_scan")
  }

  if (!x$annotated) {
    stop("x is not annotated. Annotate using annotate()")
  }

  d <- as.matrix(pdist::pdist(
    as.matrix(x$data[c("umap1", "umap2")]),
    as.matrix(deepscanscape::deep_landscape[c("umap1", "umap2")])
  ))
  min_d <- apply(d, 1, min)
  return(which(min_d > cutoff))
}
