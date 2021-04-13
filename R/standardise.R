# Normalise and standardise DMS scores

#' Normalise deep mutational scanning fitness scores
#'
#' Normalise log2(ER) scores by dividing my the median of the lowest 10% of score, under the assumption these will
#' represent complete protein knockout and will therefore be somewhat comparable between studies.
#'
#' @param x Vector of fitness scores
#' @param q Quantile to normalise against. Values other than 0.1 differ from
#' those used in the combined dataset
#' @return Numeric vector of normalised ER scores
#' @export
normalise_er <- function(x, q=0.1) {
  q <- stats::quantile(x, q, na.rm = TRUE)
  return(x / -stats::median(x[x <= q], na.rm = TRUE))
}

#' Transform deep_mutational scanning scores onto a standard scale
#'
#' Transform deep mutational scanning fitness scores onto a standard scale aligned to the log2 enrichment ratio (ER)
#' scale. This scale has neutral variants at 0, deleterious variants with negative scores and beneficial variants with
#' positive scores. Standard transforms are provided and common transforms can be passed for more unusual ones.
#'
#' The supported transforms are:
#' \itemize{
#'   \item log2: log2 transform to change raw read ratios onto the log2 ER scale.
#'   \item \code{\link[=transform_vamp]{vamp-seq}}: Transform scores produced by the VAMP-Seq method. These have the
#'   median null mutation at 0, neutral variants at 1 and positive variants greater than 1.
#'   \item invert: Multiply scores by -1, for example if the experimental fitness is inverted compared to natural
#'   fitness
#' }
#'
#' @param x Vector of fitness scores
#' @param trans Transform to apply (see description). Accepts either a string or a function.
#' @return Numeric vector of transformed scores
#' @export
transform_er <- function(x, trans = c("log2", "vamp-seq", "invert")) {
  if (is.character(trans)) {
    trans <- match.arg(trans)
    methods <- c(`vamp-seq` = transform_vamp, `log2` = log2, `invert` = function(x) -x)
    f <- methods[[trans]]
  } else if (is.function(trans)) {
    f <- trans
  } else {
    stop("Unrecognised transformation. Pass a supported string or a function")
  }

  return(f(x))
}

#' Transform VAMP-Seq fitness scores
#'
#' Transform VAMP-Seq fitness scores to the standard log2(ER) scale. The VAMP-Seq scores have the median null variant
#' at 0, neutral variants at 1 and beneficial variant >1.
#'
#' @param x Vector of fitness scores
#' @return Numeric vector of transformed scores
#' @export
transform_vamp <- function(x) {
  y <- 1 + (x - 1) / -min(x - 1, na.rm = TRUE)
  return(log2(y + min(y[y > 0], na.rm = TRUE)))
}
