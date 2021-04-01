# Annotate positions based on the overall dataset

#' Annotate positions based on the deep mutational landscape
#'
#' Add PCA and UMAP transformed coordinates and amino acid subtypes based on analysis and clustering of the combined
#' \link[=deep_landscape]{deep mutational landscape dataset}.
#'
#' PCA and UMAP coordinates are assigned using the original models fit on the
#' \link[=deep_landscape]{deep mutational landscape dataset}. The PCA coordinates are then used to assign
#' each position to an amino acid cluster. These are marked X1-8 for the main described clusters of amino acid X, XP
#' for the permissive cluster (|ER| < 0.4 for all positions), XO for outliers and XA when assignment is ambiguous.
#'
#' Clusters are initially assigned based on cosine distance to cluster centroids. Positions that are a similar distance
#' from multiple centroids (difference < 0.03) are marked as ambiguous, those that are > 0.45 away from all clusters
#' marked as outliers and those with all |ER| scores < 0.4 marked as permissive. The first two thresholds were
#' determined by benchmarking assignment on positions in the original dataset and the latter is the the threshold used
#' for permissive clusters in the original dataset.
#'
#' @param x \code{\link{deep_mutational_scan}} to annotate
#' @return An annotated \code{\link{deep_mutational_scan}} whose data object contains the following added columns:
#' \itemize{
#'   \item cluster: the assigned amino acid subtype.
#'   \item PC1 - PC20: Principal Component coordinates.
#'   \item umap1/2: UMAP coordinates.
#'   \item base_cluster: The nearest primary cluster centroid (i.e. not outlier or permissive clusters).
#'   \item permissive: The position is identified as permissive (|ER| < 0.4 for all amino acids).
#'   \item ambiguous: The distance to two clusters is very similar, so the assignment is low confidence and marked as
#'   an ambiguous.
#'   \item high_distance: The position is distant from all cluster centroids and is marked an outlier.
#'   \item dist1-8: The distance to each cluster of the WT amino acid.
#'   \item cluster_notes: Notes on the cluster assignment.
#' }
#' @examples
#' dms <- annotate(deepscanscape::deep_scans$p53)
#' @export
annotate <- function(x) {
  if (!is.deep_mutational_scan(x)) {
    stop("x must be a deep_mutational_scan")
  }

  x <- validate_deep_mutational_scan(x)

  if (any(is.na(x$data[amino_acids]))) {
    stop("NA fitness scores present, which prevents annotation. Remove missing values by filtering or using impute().")
  }

  if (x$annotated) {
      warning("deep_mutational_scan already annotated. Clearing annotation and reapplying", immediate. = TRUE)
      x$data <- dplyr::select(x$data, -dplyr::starts_with("PC"), -.data$umap1, -.data$umap2, -.data$cluster,
                              -.data$base_cluster, .data$permissive, -.data$ambiguous, -.data$high_distance,
                              dplyr::starts_with("dist"), -.data$cluster_notes)
      x$annotated <- FALSE
  }

  er_mat <- as.matrix(x$data[amino_acids])

  # Map onto PCAs
  pca <- stats::predict(dms_pca, newdata = er_mat)
  x$data <- dplyr::bind_cols(x$data, tibble::as_tibble(pca))

  # Map onto UMAP space
  model <- uwot::load_uwot(system.file("extdata", "dms_umap", package = "deepscanscape"))
  umap <- uwot::umap_transform(er_mat, model = model)
  x$data[c("umap1", "umap2")] <- umap

  # Assign clusters
  distance <- calculate_cluster_distances(x$data$wt, pca[, -1])
  closest <- apply(distance, 1, which.min)

  x$data$cluster <- stringr::str_c(x$data$wt, closest)
  x$data$base_cluster <- x$data$cluster
  x$data$permissive <- permissive_positions(er_mat)
  x$data$ambiguous <- ambiguous_assignment(distance)
  x$data$high_distance <- outlier_positions(distance)

  x$data$cluster[x$data$ambiguous] <- stringr::str_c(x$data$wt[x$data$ambiguous], "A")
  x$data$cluster[x$data$high_distance] <- stringr::str_c(x$data$wt[x$data$high_distance], "O")
  x$data$cluster[x$data$permissive] <- stringr::str_c(x$data$wt[x$data$permissive], "P")

  x$data <- dplyr::bind_cols(x$data, as.data.frame(distance))
  x$data$cluster_notes <- cluster_notes(x$data$permissive, x$data$ambiguous, x$data$high_distance)

  x$annotated <- TRUE
  return(x)
}

#' Calculate distance to each cluster of the amino acid
#'
#' Calculate the distances to each cluster of that amino acid for each row of the input PCA matrix
#'
#' Internal function called by \code{\link{annotate}}.
#'
#' @param wt Wild type amino acid for each position
#' @param pca Numeric matrix of PCA values for PC2 to PC20
#' @return Matrix of cosine distances to the centroids of clusters 1 to 8 for the WT amino acid. When the wild type
#' amino acid has fewer than 8 clusters NA values occur in the additional columns
#' @keywords internal
calculate_cluster_distances <- function(wt, pca) {
  dists <- matrix(nrow = length(wt), ncol = 8)
  colnames(dists) <- stringr::str_c("dist", 1:8)

  for (aa in unique(wt)) {
    c <- cluster_centers[stringr::str_starts(rownames(cluster_centers), aa), , drop = FALSE]
    dists[wt == aa, seq_len(nrow(c))] <- cosine_distance_matrix(pca[wt == aa, , drop = FALSE], c)
  }

  return(dists)
}

#' Identify ambiguous cluster assignment
#'
#' Identify positions where the difference between the minimum distance to a cluster and the distance any other is
#' less than 0.03, which indicates ambiguous assignment.
#'
#' Internal function called by \code{\link{annotate}}.
#'
#' @param x Numeric matrix of distances
#' @return Logical vector corresponding to results on each row
#' @keywords internal
ambiguous_assignment <- function(x) {
  rowSums((x - apply(x, 1, min, na.rm = TRUE)) < 0.03, na.rm = TRUE) > 1
}

#' Identify outlier positions
#'
#' Identify positions where the distance to all clusters is > 0.45, indicating a likely outlier.
#'
#' Internal function called by \code{\link{annotate}}.
#'
#' @param x Numeric matrix of distances
#' @return Logical vector corresponding to results on each row
#' @keywords internal
outlier_positions <- function(x) {
  apply(x, 1, min, na.rm = TRUE) > 0.45
}

#' Generate cluster notes message
#'
#' Create a standard note for cluster assignments, explaining why O, P or A cluster status is assigned.
#'
#' Internal function called by \code{\link{annotate}}.
#'
#' @param permissive Logical showing which positions are permissive
#' @param ambiguous Logical showing which positions are ambiguous
#' @param high_distance Logical showing which positions are outliers due to distance from all centroids.
#' @return Character vector
#' @keywords internal
cluster_notes <- function(permissive, ambiguous, high_distance) {
  notes <- rep(NA, length(permissive))
  notes[ambiguous] <- "Top cluster is only nearer by a small margin, making assignment ambiguous"
  notes[high_distance] <- "Not close to any cluster center"
  notes[ambiguous & high_distance] <- "Both distant and only nearer one by a small margin"
  notes[permissive] <- "No mutation with |ER| > 0.4"
  return(notes)
}

#' Predict positions properties from their position in the deep mutational landscape
#'
#' Predict positions properties based on the average properties of \code{\link{deep_landscape}} positions in the same
#' hexagonal bin as them.
#'
#' This function is not externally exposed or used for other annotation because subtype based annotation is found to be
#' at least as accurate in almost all cases, and often better. It is included for illustration and potential further
#' development work.
#'
#' @param x \link{deep_mutational_scan}.
#' @param bins Number of bins to divide landscape into. The default value was chosen to give the best results on a
#' benchmark dataset.
#' @returns A \code{\link[tibble]{tibble}}, with each row detailing a row of the input data and columns matching
#' those in the \code{\link{subtypes}} dataset, describing the predicted properties of the position.
#' @keywords internal
landscape_properties <- function(x, bins = 20) {
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop("Install hexbin to perform hex based binning")
  }

  if (!is.deep_mutational_scan(x)) {
    stop("x must be a deep_mutational_scan")
  }

  x <- validate_deep_mutational_scan(x)

  if (!x$annotated) {
    warning("deep_mutational_scan is not annotated. Annotating using annotate_dms().", immediate. = TRUE)
    x <- annotate(x)
  }

  df <- x$data
  xbnds <- range(c(deepscanscape::deep_landscape$umap1, df$umap1))
  ybnds <- range(c(deepscanscape::deep_landscape$umap2, df$umap2))

  hex <- hexbin::hexbin(deepscanscape::deep_landscape$umap1, deepscanscape::deep_landscape$umap2,
                        xbins = bins, IDs = TRUE, xbnds = xbnds, ybnds = ybnds)

  hex_summary <- dplyr::group_by(dplyr::mutate(deepscanscape::deep_landscape, hex = hex@cID), .data$hex)
  hex_summary <- dplyr::summarise(hex_summary,
                                  dplyr::across(c(.data$mean_sift:.data$energy_ionisation, .data$all_atom_rel),
                                                mean, na.rm = TRUE),
                                  .groups = "drop")

  new_hexes <- hexbin::hexbin(df$umap1, df$umap2, xbins = bins, xbnds = xbnds, ybnds = ybnds, IDs = TRUE)
  df <- df[, c("name", "position", "wt")]
  df$hex <- new_hexes@cID
  df <- dplyr::select(dplyr::left_join(df, hex_summary, by = "hex"), -.data$hex)
  return(df)
}

#' Describe amino acid positional subtypes
#'
#' Add descriptions and frequencies from the cluster assigned to each position in an annotated deep_mutational_scan
#' dataset. These were determined through analysis of the larger \code{\link{deep_landscape}} dataset. Numerical
#' summary statistics based on the same dataset can also be added, which give the mean characteristics of the subtype.
#' They include results from FoldX and Naccess and mean ER fitness score profiles.
#'
#' The deep landscape properties summarised where chosen by benchmarking subtype based predictions against the
#' landscape data and a separate deep mutational scan of SARS-CoV-2 Spike by Starr et al. (2020). Metrics are only shown
#' when the average value of positions of that subtype is meaningfully related to the observed values.
#'
#' Properties are not shown for outlier subtypes, which do not have consistent properties, and ambiguously assigned
#' positions, since the properties could be that of either of the subtypes the position could have been assigned to.
#'
#' @param x \link{deep_mutational_scan}.
#' @param full Logical. Include average statistics from the \code{\link{deep_landscape}} dataset in
#' addition to summary descriptions and notes.
#' @return A \code{\link[tibble]{tibble}}, with each row detailing a row of the input data and columns matching
#' those in the \code{\link{subtypes}} dataset.
#' @examples
#' dms <- deepscanscape::deep_scans$p53
#'
#' # Basic description
#' describe_clusters(dms)
#'
#' # More details
#' describe_clusters(dms, full = TRUE)
#'
#' @export
describe_clusters <- function(x, full = FALSE) {
  if (!is.deep_mutational_scan(x)) {
    stop("x must be a deep_mutational_scan")
  }

  x <- validate_deep_mutational_scan(x)

  if (!x$annotated) {
    warning("deep_mutational_scan is not annotated. Annotating using annotate_dms().", immediate. = TRUE)
    x <- annotate(x)
  }

  df <- x$data

  cols <- c("wt", "cluster", "prop", "group", "description", "notes", amino_acids,
            "mean_score", "total_energy", "backbone_hbond", "sidechain_hbond",
            "van_der_waals", "electrostatics", "solvation_polar", "solvation_hydrophobic",
            "van_der_waals_clashes", "entropy_sidechain", "entropy_mainchain", "cis_bond",
            "torsional_clash", "backbone_clash", "disulfide", "partial_covalent_bonds",
            "energy_ionisation", "all_atom_rel")
  extra <- dplyr::select(deepscanscape::subtypes, dplyr::all_of(cols))
  extra <- dplyr::rename(extra, global_cluster_freq = .data$prop)
  if (!full) {
    extra <- extra[c("wt", "cluster", "global_cluster_freq", "group", "description", "notes")]
  }

  df <- dplyr::left_join(df[c("name", "position", "wt", "cluster")], extra, by = c("wt", "cluster"))
  return(df)
}
