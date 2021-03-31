# Perform a full analysis on user supplied data and generate

# TODO show total energy etc?
#' Generate an HTML report for new deep mutational scan data
#'
#' Perform a full analysis on new deep mutational scan data and summarise the results as an HTML report. The
#' report provides an overview of the new data and then compares the properties of the new data to the datasets
#' included in the \code{\link{deep_landscape}} dataset. This helps to identify unusual properties, map the
#' features of the protein on the common mutational landscape and assign each position to an amino acid subtype, which
#' informs on its properties.
#'
#' @param ... \code{\link{deep_landscape}} objects or lists of such objects to report on.
#' @param output_path Path to write report to.
#' @param highlight A named list of positions to highlight, where each element is named for a study in ... and contains
#' the positions in that study to highlight.
#' @returns Path to the generated report.
#' @examples
#' \dontrun{
#'   # For a single study
#'   generate_report(deep_scans$p53)
#'
#'   # For multiple studies
#'   highlight <- list(`Elazar GpA`=c(1,2,3), `Kotler p53`=c(50, 100, 110))
#'   generate_report(deep_scans, highlight = highlight)
#' }
#' @export
generate_report <- function(..., output_path = "deepscanscape_report.html", highlight = NULL) {
  scans <- rlang::flatten(list(...))

  if (length(scans) == 1) {
    dms <- scans[[1]]
  } else {
    dms <- bind_scans(scans)
  }

  if (!is.null(highlight) & !is.list(highlight)) {
    stop("highlight must be NULL or a list of position vectors to select")
  }

  rmarkdown::render(system.file("template", "dms_report_template.Rmd", package = "deepscanscape"),
                    output_dir = dirname(output_path),
                    output_file = basename(output_path),
                    clean = TRUE,
                    envir = new.env(),
                    params = list(dms = dms, highlight = highlight))
}
