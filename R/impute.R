# Impute missing DMS data

# TODO - track which values have been imputed
#' Impute missing deep mutational scan data
#'
#' @param x link{deep_mutational_scan} to impute data from
#' @param na_value Value to set missense NA values to. Setting this to "impute" sets scores to the median score for that
#' substitution (e.g A -> C) in the combined dataset
#'
#' @export
impute_dms <- function(x, na_value="impute") {
  if (!"deep_mutational_scan" %in% class(x)) {
    stop("Unrecognised data.\nCreate a standardised dataset using deep_mutational_scan()")
  }

  df <- tidyr::pivot_longer(x$data, dplyr::one_of(amino_acids), names_to = "mut", values_to = "score")
  df$score[df$wt == df$mut & is.na(df$score)] <- 0

  if (na_value == "impute") {
    df$score[is.na(df$score)] <- median_scores[as.matrix(df[is.na(df$score), c("wt", "mut")])]
  } else {
    df$score[is.na(df$score)] <- na_value
  }

  df <- tidyr::pivot_wider(df, names_from = "mut", values_from = "score")
  x$data <- dplyr::select(df, .data$position, .data$wt, amino_acids, dplyr::everything())
  return(x)
}
