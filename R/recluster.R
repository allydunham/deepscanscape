# Recluster the whole dataset with new data and report on the output results

#' Recluster DMS data with a combined dataset
#'
#' @param ... \link{deep_mutational_scan} objects to recluster
#'
#' @export
recluster_dms <- function(...) {
  if (!all(sapply(list(...), is.deep_mutational_scan))) {
    stop("Datasets to recluster must all be deep mutational scans.\nGenerate these using deep_mutational_scan()")
  }

  df <- dplyr::select(deepscanscape::deep_mutational_scans, .data$study, .data$gene,
                      .data$position, .data$wt, dplyr::one_of(amino_acids))
  df <- dplyr::bind_rows(rbind(...), df)

  # Perform PCA

  # Identify Permissive Positions

  # Hierarchical Clustering

  # Dynamic Tree Cut
  # TODO how to do deep split parameter?

  # TODO implement recluster
  stop("Not implemented yet")
}
