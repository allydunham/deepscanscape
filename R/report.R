# Perform a full analysis on user supplied data and generate

#' Generate an HTML report for new deep mutational scan data
#'
#' Perform a full analysis on new deep mutational scan data and summarise the results as an HTML report. The
#' report provides an overview of the new data and then compares the properties of the new data to the datasets
#' included in the \code{\link{deep_landscape}} dataset. This helps to identify unusual properties, map the
#' features of the protein on the common mutational landscape and assign each position to an amino acid subtype, which
#' informs on its properties.
#'
#' @param ... \code{\link{deep_landscape}} objects to report.
#' @param output_path Path to write report to.
#' @returns Path to the generated report.
#' @export
generate_report <- function(..., output_path = "deepscanscape_report.html") {
  rmd_file <- system.file("template", "dms_report_template.Rmd", package = "deepscanscape")
  params <- list(studies = list(...))
  output_dir <- dirname(output_path)
  rmarkdown::render(rmd_file, output_dir = output_dir, output_file = output_path,
                    clean = TRUE, params = params, envir = new.env())
}
