# Make summary objects for deep mutational scan objects

#' Summarise Deep Mutational Scans
#'
#' @param object,x A \code{\link{deep_mutational_scan}} object.
#' @param ... Additional arguments
#' @return deep_mutational_scan_summary object containing:
#'   \describe{
#'     \item{studies}{The number of studies scanned}
#'     \item{annotated}{Logical showing whether the dataset has been annotated with deep mutational landscape data}
#'     \item{imputed}{Named vector summarising the number of imputed scores}
#'     \item{multi_study}{Logical showing whether the dataset contains multiple studies}
#'     \item{positions}{Number of positions with data}
#'     \item{er}{Named vector containing summary statistics to ER scores}
#'     \item{clusters}{Named vector summarising cluster annotations}
#'   }
#' @name dms_summary
#' @export
summary.deep_mutational_scan <- function(object, ...) {
  out <- as.list(object)[c("meta", "annotated", "imputed", "multi_study")]

  # Basic overview
  out$positions <- dim(object)[1]
  out$studies <- nrow(object$meta)

  # ER values
  er <- as.vector(as.matrix(object$data[amino_acids]))
  q <- unname(stats::quantile(er))
  out$er <- c(Min = q[1], `1st Q.` = q[2], Mean = mean(er), Median = stats::median(er), `3rd Q.` = q[4], Max = q[5])

  # Imputation
  impute_mask <- as.matrix(dplyr::select(object$data, dplyr::starts_with("impute")))
  out$imputed <- c(Mutant = sum(impute_mask == 2, na.rm = TRUE),
                   WildType = sum(impute_mask == 1, na.rm = TRUE),
                   `Mean per position` = mean(rowSums(impute_mask > 0)))

  if (out$annotated) {
    out$clusters <- list(ambiguous = sum(stringr::str_ends(object$data$cluster, "A")),
                         permissive = sum(stringr::str_ends(object$data$cluster, "P")),
                         outlier = sum(stringr::str_ends(object$data$cluster, "O")),
                         subtype = sum(!stringr::str_detect(object$data$cluster, "[A-Z][AOP]{1}")),
                         nclusters = length(unique(object$data$cluster[!stringr::str_ends(object$data$cluster, "O")])))
  }

  class(out) <- "deep_mutational_scan_summary"
  return(out)
}

#' @describeIn dms_summary S3 format method
#' @export
format.deep_mutational_scan_summary <- function(x, ...) { # nolint
  stusdy_str <- stringr::str_c(x$meta$study, " - ", x$meta$gene, " (", x$meta$positions,
                               " positions, name: ", x$meta$name, ")")
  out <- c(
    stringr::str_c("Deep mutational scan dataset of", x$positions, "positions from", x$studies,
                   ifelse(x$studies == 1, "study", "studies"), sep = " "),
    "",
    stusdy_str,
    "",
    "ER Scores:",
    utils::capture.output(print(x$er)),
    "",
    "Imputed Scores:",
    utils::capture.output(print(x$imputed))
  )

  if (x$annotated) {
    out <- c(
      out,
      "",
      "Annotated with amino acid subtypes, pricipal components and UMAP coordinates",
      stringr::str_c("  Permissive positions: ", x$clusters$permissive),
      stringr::str_c("  Ambiguous positions: ", x$clusters$ambiguous),
      stringr::str_c("  Outlier positions: ", x$clusters$outlier),
      stringr::str_c("  Subtype positions: ", x$clusters$subtype),
      stringr::str_c("  Subtypes represented (incl. permissive): ", x$clusters$nclusters)
    )
  }

  return(out)
}

#' @describeIn dms_summary S3 print method
#' @export
print.deep_mutational_scan_summary <- function(x, ...) { # nolint
  cat(format(x, ...), sep = "\n")
}
