# Parse common input formats in preparation for analysis

#' Parse deep mutational scanning data
#'
#' Processes deep mutational scanning data from recognised formats and convert it into a standardised data frame, with
#' one row per variant and position, wt, mut and score columns. Duplicate scores can also be averaged to account for
#' multiple replicates and single variant scores averaged from multiply mutated sequences for some input formats.
#'
#' The following input formats are supported:
#' \itemize{
#'   \item \code{\link[=parse_mavedb]{mavedb}}: Data downloaded from \url{mavedb.org}. This format additionally supports
#'   averaging over multiply mutated sequences.
#'   \item long: Data in long format, which already has the required columns (position, wt, mut and score). This allows
#'   duplicates to be averaged conveniently if data is already in the standard format.
#'   \item \code{\link[=parse_wide]{wide}}: Process wide format data, with columns for position, wt and scores for each
#'   amino acid
#'   \item \code{\link[=parse_sequence]{sequence}}: Each score is associated with a sequence, from which variants are
#'   extracted. This method supports averaging over multiply mutated sequences. It is the most frequent use for the
#'   position_offset parameter.
#' }
#'
#' @param x A tibble or an object that can be converted to one
#' @param scheme Original data scheme (see description). Accepts any unambiguous substring.
#' @param position_offset Offset all position
#' @param duplicates How to handle duplicate scores for a single mutation
#' @param synonymous Character vector of strings to consider as synonymous variants. Amino acid letters are silently filtered
#' @param stop_codons Character vector of strings to consider as stop codons. Amino acid letters are silently filtered
#' @param ... Arguments passed to the chosen parsing function, including specifying parsing for multiply mutated
#' sequences (see individual methods for details).
#' @return A long format \code{\link[tibble]{tibble}} with columns for 'position', 'wt', 'mut' and 'score'.
#' @export
parse_deep_scan <- function(x, scheme = c("mavedb", "long", "wide", "sequence"), position_offset = 0,
                            duplicates = c("warn", "mean", "median", "error"), synonymous = c("=", "sym"),
                            stop_codons = c("*", "Ter", "ter", "stop"), ...) {
  scheme <- match.arg(scheme)
  methods <- c("mavedb" = parse_mavedb, "long" = identity,
               "wide" = parse_wide, "sequence" = parse_sequence)
  synonymous <- synonymous[!synonymous %in% amino_acids]
  stop_codons <- stop_codons[!stop_codons %in% amino_acids]
  
  x <- tibble::as_tibble(x)
  f <- methods[[scheme]]
  x <- f(x, ...)

  # Handle duplicate positions
  duplicates <- match.arg(duplicates)
  if (nrow(x) > nrow(dplyr::distinct(x, .data$position, .data$wt, .data$mut))) {
    if (duplicates == "error") {
      stop("Duplicate scores present - control behaviour with the 'duplicates' parameter")
    } else if (duplicates == "warn") {
      warning("Duplicate scores present - averaging scores for each variant using 'mean'", immediate. = TRUE)
      duplicates <- "mean"
    }

    f <- match.fun(duplicates)
    x <- dplyr::group_by(x, .data$position, .data$wt, .data$mut)
    x <- dplyr::summarise(x, score = f(.data$score, na.rm = TRUE), .groups = "drop")
  }

  x$position <- x$position + position_offset

  # Convert synonymous variants to the corresponding wt
  x$mut[x$mut %in% synonymous] <- x$wt[x$mut %in% synonymous]
  
  # Normalise stop codons to stop (avoiding * symbol column name)
  x$mut[x$mut %in% stop_codons] <- "stop"
  
  # Check amino acids - filter any that are unrecognised
  x <- dplyr::filter(x, .data$wt %in% amino_acids, .data$mut %in% c(amino_acids, "stop"))

  return(x)
}

#' Parse data in the MaveDB format
#'
#' Parse data downloaded from \href{mavedb.org}{MaveDB}, including averaging across multiply mutated sequences.
#' Usually called internally by\code{\link{parse_deep_scan}} but is exposed to help users get their input into the
#' correct format.
#'
#' @param x A data frame with a column 'hgvs_pro' gicing the HGVS protein mutation string describing the variant(s) and
#' a fitness score column ('score' by default).
#' @param score_col String. Column containing fitness scores, to conveniently use an additional data column as the score
#' where mutliple measurements are included in the data.
#' @param average_multi Average scores for variants included in multiply mutated sequences, where they have not been
#' measured individually. Care should be taken to check that the type of multiple mutation makes this appropriate.
#' @param ... Ignored.
#' @return A long format \code{\link[tibble]{tibble}} with columns specifying 'position', 'wt', 'mut' and 'score'.
#' @export
parse_mavedb <- function(x, score_col = "score", average_multi = FALSE, ...) {
  # Regexes
  re_any_sub <- "p\\.\\[?([A-Z][a-z]{2}[0-9]*([A-Z][a-z]{2}|=);?)+\\]?"
  re_multi_sub <- "p\\.\\[([A-Z][a-z]{2}[0-9]*([A-Z][a-z]{2}|=);?)*\\]"
  re_single_sub <- "p\\.([A-Z][a-z]{2})([0-9]*)([A-Z][a-z]{2}|=|\\*)"

  # Filter tibble to substitutions
  x <- tibble::as_tibble(x)
  x <- dplyr::filter(x[c("hgvs_pro", score_col)], stringr::str_detect(.data$hgvs_pro, re_any_sub))

  if (any(stringr::str_detect(x$hgvs_pro, re_multi_sub))) {
    if (average_multi) {
      x_multi <- dplyr::filter(x, stringr::str_detect(.data$hgvs_pro, re_multi_sub))
      x_multi$hgvs_pro <- stringr::str_split(stringr::str_remove_all(x_multi$hgvs_pro, "(p\\.\\[|\\])"), ";")
      x_multi <- tidyr::unnest(x_multi, .data$hgvs_pro)
      x_multi$hgvs_pro <- stringr::str_c("p.", x_multi$hgvs_pro)
      x_multi <- dplyr::summarise(dplyr::group_by(x_multi, .data$hgvs_pro), score = mean(.data$score))

      x <- dplyr::filter(x, stringr::str_detect(.data$hgvs_pro, re_single_sub))
      x <- dplyr::bind_rows(x, x_multi[!x_multi$hgvs_pro %in% x$hgvs_pro, ])

    } else {
      warning("Multiply mutated sequences detected. These are currently ignored, if you want to average these scores ",
              "to use for variants that have not been measured alone use average_multi = TRUE", immediate. = TRUE)
      x <- dplyr::filter(x, stringr::str_detect(.data$hgvs_pro, re_single_sub))
    }
  }

  x[c("wt", "position", "mut")] <- stringr::str_match(x$hgvs_pro, re_single_sub)[, 2:4]
  x$wt <- unname(aa_3_to_1[stringr::str_to_lower(x$wt)])
  x$mut <- unname(aa_3_to_1[stringr::str_to_lower(x$mut)])
  x$position <- as.numeric(x$position)
  x <- dplyr::rename(x, score = .data[[score_col]])
  return(x[c("position", "wt", "mut", "score")])
}

#' Parse data in the wide format
#'
#' Parse wide format data into the standard layout. The input data frame must have columns specifying 'position',
#' 'wt', and ER score ('score') for each amino acid (named 'A', 'C', ...). Usually called internally by
#' \code{\link{parse_deep_scan}} but is exposed to help users get their input into the correct format.
#'
#' @param x A data frame.
#' @param ... Ignored.
#' @return A long format \code{\link[tibble]{tibble}} with columns specifying each variants 'position', 'wt', 'mut' and
#' 'score'.
#' @export
parse_wide <- function(x, ...) {
  x <- tibble::as_tibble(x)[, c("position", "wt", amino_acids)]
  x <- tidyr::pivot_longer(x, dplyr::all_of(amino_acids), names_to = "mut", values_to = "score")
  return(x)
}

#' Parse data in the sequence format
#'
#' Extract variants from sequences and ER scores. The most frequent amino acid at each position in the sequence is
#' assumed to be the wild type if a wild type sequence is not provided. Usually called internally by
#' \code{\link{parse_deep_scan}} but is exposed to help users get their input into the correct format.
#'
#' @param x A data frame with columns 'sequence' and 'score'.
#' @param wt_seq Character vector. Wild type sequence to calculate variants against. If NULL it is assumed to be the
#' most common amino acid in each position.
#' @param average_multi Average scores for variants included in multiply mutated sequences, where they have not been
#' measured individually. Care should be taken to check that the type of multiple mutation makes this appropriate.
#' @param ... Ignored.
#' @return A long format data frame with columns for 'position', 'wt', 'mut' and 'score'.
#' @export
parse_sequence <- function(x, wt_seq = NULL, average_multi = FALSE, ...) {
  x <- tibble::as_tibble(x)[, c("sequence", "score")]
  seqs <- stringr::str_split(x$sequence, "", simplify = TRUE)

  if (is.null(wt_seq)) {
    wt_seq <- apply(seqs, 2, function(x) names(which.max(table(x))))
  } else if (length(wt_seq) == 1) {
    warning("length(wt_seq) == 1, splitting into a character vector with one amino acid per element")
    wt_seq <- stringr::str_split(wt_seq, "", simplify = TRUE)
  }

  pos <- seq_along(wt_seq)
  x$mut <- apply(seqs, 1, function(x) stringr::str_c(wt_seq, pos, x)[!x == wt_seq])
  x <- tidyr::unnest(x, .data$mut)

  x <- dplyr::add_count(x, .data$sequence)

  if (any(x$n > 1)) {
    if (average_multi) {
      x_multi <- x[x$n > 1, ]
      x_multi <- dplyr::summarise(dplyr::group_by(x_multi, .data$mut), score = mean(.data$score))
      x <- x[x$n == 1, c("mut", "score")]
      x <- dplyr::bind_rows(x, x_multi[!x_multi$mut %in% x$mut, ])

    } else {
      warning("Multiply mutated sequences detected. These are currently ignored, if you want to average these scores ",
              "to use for variants that have not been measured alone use average_multi = TRUE", immediate. = TRUE)
      x <- x[x$n == 1, c("mut", "score")]
    }
  }
  x <- tidyr::extract(x, .data$mut, into = c("wt", "position", "mut"),
                      regex = "([A-Z])([0-9]*)([A-Z\\*])", convert = TRUE)
  return(x[c("position", "wt", "mut", "score")])
}
