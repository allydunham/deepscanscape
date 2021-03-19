# Parse common input formats in preparation for analysis
# TODO - Test funcs
# TODO - Add additional formats
# TODO - support data with multiple scores for each variant
# TODO - support data with multiple substitutions - take average of occurrence?

#' Parse deep mutational scanning data
#'
#' @param x A tibble or an object that can be converted to one
#' @param scheme Original data scheme (see description). Accepts any unambiguous
#' substring.
#' @param position_offset Offset all positions by this amount
#' @param duplicates How to handle duplicate scores for a single mutation
#' @param ... Arguments passed to the chosen parsing function
#'
#' @export
parse_deep_scan <- function(x, scheme = c("mavedb", "long", "wide", "sequence"), position_offset = 0,
                            duplicates = c("warn", "mean", "median", "error"), ...) {
  scheme <- match.arg(scheme)
  methods <- c("mavedb" = parse_mavedb, "long" = parse_long,
               "wide" = parse_wide, "sequence" = parse_sequence)

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
  return(x)
}

#' Parse data in the MaveDB format (called from parse_dms_data)
#'
#' @param x A tibble.
#' @param score_col String. Column containing fitness scores.
#' @param ... Ignored.

#' @export
parse_mavedb <- function(x, score_col = "score", ...) {
  x <- x[c("hgvs_pro", score_col)]

  # Filter tibble to substitutions
  x <- dplyr::filter(x, stringr::str_detect(.data$hgvs_pro, "p\\.[A-Z][a-z]{2}[0-9]*([A-Z][a-z]{2}|=)"))

  x[c("wt", "position", "mut")] <- stringr::str_match(x$hgvs_pro, "p\\.([A-Za-z=\\*]*)([0-9]*)([A-Za-z=\\*]*)")[, 2:4]
  x$wt <- unname(aa_3_to_1[stringr::str_to_lower(x$wt)])
  x$mut <- unname(aa_3_to_1[stringr::str_to_lower(x$mut)])
  x$position <- as.numeric(x$position)
  x <- dplyr::rename(x, score = .data[[score_col]])
  return(x[c("position", "wt", "mut", "score")])
}

#' Parse data in the long format (called from parse_dms_data)
#'
#' @param x A tibble.
#' @param ... Ignored.
parse_long <- function(x, ...) {
  stop("Not implemented yet")
}

#' Parse data in the wide format (called from parse_dms_data)
#'
#' @param x A tibble.
#' @param ... Ignored.
parse_wide <- function(x, ...) {
  stop("Not implemented yet")
}

#' Parse data in the sequence format (called from parse_dms_data)
#'
#' @param x A tibble.
#' @param ... Ignored.
parse_sequence <- function(x, ...) {
  stop("Not implemented yet")
}
