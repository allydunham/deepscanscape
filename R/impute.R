# Impute missing DMS data
# TODO test
# TODO allow passing custom impute er scores

#' Impute missing deep mutational scan data
#'
#' @param x \code{\link{deep_mutational_scan}} to impute data from
#' @param na_value Value to set missense NA values to. Setting this to "impute" sets scores to the median score for that
#' substitution (e.g A -> C) in the combined dataset
#'
#' @export
impute <- function(x, na_value="impute") {
  if (!is.deep_mutational_scan(x)) {
    stop("Unrecognised data.\nCreate a standardised dataset using deep_mutational_scan()")
  }

  df <- tidyr::pivot_longer(x$data, dplyr::one_of(amino_acids), names_to = "mut", values_to = "score")

  # Calculate impute mask
  # TODO document the structure of this
  mask <- df
  mask$mask <- 0
  mask$mask[is.na(mask$score) & mask$wt == mask$mut] <- 1
  mask$mask[is.na(mask$score) & mask$wt != mask$mut] <- 2
  mask <- tidyr::pivot_wider(mask[c("position", "wt", "mut", "mask")], names_from = .data$mut,
                             values_from = .data$mask, names_prefix = "impute_")


  df$score[df$wt == df$mut & is.na(df$score)] <- 0

  if (na_value == "impute") {
    df$score[is.na(df$score)] <- median_scores[as.matrix(df[is.na(df$score), c("wt", "mut")])]
  } else {
    df$score[is.na(df$score)] <- na_value
  }

  df <- tidyr::pivot_wider(df, names_from = "mut", values_from = "score")
  df <- dplyr::bind_cols(df, dplyr::select(mask, dplyr::starts_with("impute_")))
  x$data <- df
  x$imputed <- TRUE
  return(x)
}
