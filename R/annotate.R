# Annotate positions based on the overall dataset

#' Annotate positions based on the deep mutational landscape
#'
#' @param x \link{deep_mutational_scan} to annotate
#'
#' @export
annotate_dms <- function(x) {
  if (!"deep_mutational_scan" %in% class(x)) {
    stop("Unrecognised data.\nCreate a standardised dataset using deep_mutational_scan()")
  }

  if (any(is.na(x[amino_acids]))) {
    stop("NA fitness scores present\nUse impute_dms")
  }

  # Assign subtypes
  # TODO

  # Map onto PCAs
  pca <- stats::predict(dms_pca, newdata = as.matrix(x[amino_acids]))
  x$data <- dplyr::bind_cols(x$data, tibble::as_tibble(pca))

  # Map onto UMAP space
  model <- uwot::load_uwot(system.file("extdata", "dms_umap", package = "deepscanscape"))
  umap <- uwot::umap_transform(as.matrix(x[amino_acids]), model = model)
  x[c("umap1", "umap2")] <- umap

  x$annotated <- TRUE
  return(x)
}
