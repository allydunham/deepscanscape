# Impute missing DMS data

#' Impute missing deep mutational scan data
#'
#' Impute NA values in a deep_mutational_scan objects ER scores. Values can be imputed as the median scores from the
#' combined landscape dataset or using custom values.
#'
#' If na_value == "impute" missense NA scores are imputed to be the median value of that substitution (e.g. A -> C) from
#' the \code{\link{deep_landscape}} dataset. If na_value == "average" missense NA scores are set to the average missense
#' score for that position. If na_value is a matrix it should have rows and column names corresponding
#' to single letter amino acid codes and have cell i,j correspond to the imputed score for substitutions from i to j.
#' Any other value of na_value is interpreted as the score to impute all NA values to.
#'
#' An impute mask is also generated and added to the data tibble of the deep_mutational_scan. This consists of one
#' column for each amino acid (impute_X) which contains '0' if the corresponding score in that row is not imputed, '1'
#' for synonymous substitution imputed as 0 and '2' for non-synonymous substitutions that have undergone imputation.
#'
#' @param x \code{\link{deep_mutational_scan}} to impute data from
#' @param na_value Value to set missense NA values to. Can be "impute", "average", a single value or a matrix of
#' substitution scores (see details).
#' @return An imputed \code{\link{deep_mutational_scan}}
#' @examples
#' # Load an unimputed DMS object
#' path <- system.file("extdata", "urn_mavedb_00000011_a_1_scores.csv",
#'                     package = "deepscanscape")
#' csv <- read.csv(path, skip = 4)
#' dms <- deep_mutational_scan(csv, name = "Hietpas Hsp90", scheme = "mave", trans = NULL,
#'                             na_value = NULL, annotate = FALSE)
#'
#' # Set NA to a constant
#' impute(dms, na_value = 1)
#'
#' # Use the built in imputed values
#' impute(dms, na_value = "impute")
#'
#' @export
impute <- function(x, na_value="impute") {
  if (!is.deep_mutational_scan(x)) {
    stop("Unrecognised data.\nCreate a standardised dataset using deep_mutational_scan()")
  }

  x <- validate_deep_mutational_scan(x)

  if (x$imputed) {
    stop("x has already been imputed")
  }

  df <- tidyr::pivot_longer(x$data, dplyr::one_of(amino_acids), names_to = "mut", values_to = "score")

  # Calculate impute mask
  mask <- df
  mask$mask <- 0
  mask$mask[is.na(mask$score) & mask$wt == mask$mut] <- 1
  mask$mask[is.na(mask$score) & mask$wt != mask$mut] <- 2
  mask <- tidyr::pivot_wider(mask[c("name", "position", "wt", "mut", "mask")], names_from = .data$mut,
                             values_from = .data$mask, names_prefix = "impute_")

  df$score[df$wt == df$mut & is.na(df$score)] <- 0

  if (na_value == "impute") {
    df$score[is.na(df$score)] <- median_scores[as.matrix(df[is.na(df$score), c("wt", "mut")])]

  } else if (na_value == "average") {
    df <- dplyr::group_by(df, .data$position)
    df <- dplyr::mutate(df, mean = mean(.data$score[.data$wt != .data$mut], na.rm = TRUE))
    df <- dplyr::ungroup(df)
    df$score[is.na(df$score)] <- df$mean[is.na(df$score)]
    df <- dplyr::select(df, -.data$mean)

  } else if (is.matrix(na_value)) {
    if (!(rownames(na_value) == amino_acids & colnames(na_value) == na_value)) {
      stop("na_value must have rows and column names corresponding to the amino acids (A, C, D, etc.)")
    }

    if (!is.numeric(na_value) | any(is.na(na_value))) {
      stop("The imputation matrix must be numeric and contain no NA values")
    }

    df$score[is.na(df$score)] <- na_value[as.matrix(df[is.na(df$score), c("wt", "mut")])]

  } else if (length(na_value) == 1 & is.numeric(na_value)) {
    df$score[is.na(df$score)] <- na_value

  } else {
    stop("na_value must be 'impute', a correctly formatted matrix or a single numeric value")
  }

  df <- tidyr::pivot_wider(df, names_from = "mut", values_from = "score")
  df <- dplyr::bind_cols(df, dplyr::select(mask, dplyr::starts_with("impute_")))
  x$data <- df
  x$imputed <- TRUE
  return(x)
}
