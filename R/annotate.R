# Annotate positions based on the overall dataset
# TODO test annotation functions

#' Annotate positions based on the deep mutational landscape
#'
#' Add PCA and UMAP transformed coordinates and amino acid subtypes based on analysis and clustering of the combined
#' \link[=deep_mutational_scans]{deep mutational landscape dataset}.
#'
#' PCA and UMAP coordinates are assigned using the original models fit on the
#' \link[=deep_mutational_scans]{deep mutational landscape dataset}. The PCA coordinates are then used to assign
#' each position to an amino acid cluster. These are marked X1-8 for the main described clusters of amino acid X, XP
#' for the permissive cluster (|ER| < 0.4 for all positions), XO for outliers and XA when assignment is ambiguous.
#'
#' Clusters are initially assigned based on cosine distance to cluster centroids. Positions that are a similar distance
#' from multiple centroids (difference < 0.03) are marked as ambiguous, those that are > 0.45 away from all clusters
#' marked as outliers and those with all |ER| scores < 0.4 marked as permissive. The first two thresholds were
#' determined by benchmarking assignment on positions in the original dataset and the latter is the the threshold used
#' for permissive clusters in the original dataset.
#'
#' @param x \code{\link{deep_mutational_scan}} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}
#' to annotate
#' @return An annotated \code{\link{deep_mutational_scan}} if passed a \code{\link{deep_mutational_scan}} or a
#' \code{\link[tibble]{tibble}} when passed a \link[=rbind.deep_mutational_scan]{combined DMS data frame}.
#' This data frame (or the cluster element of the annotated scan) contains the following added columns:
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
#'   \item notes: Notes on the cluster assignment.
#' }
#' @examples
#' dms <- annotate(deepscanscape::deep_scan)
#' @export
annotate <- function(x) {
  if (is.deep_mutational_scan(x)) {
    if (any(is.na(x[amino_acids]))) {
      stop("NA fitness scores present\nUse impute to remove these")
    }

    if (x$annotated) {
      warning("deep_mutational_scan already annotated. Clearing annotation and reapplying")
      x$cluster <- NA
      x$data <- dplyr::select(x$data, -.data$cluster)
      x$annotated <- FALSE
    }

    df <- x$data[c("position", "wt", amino_acids)]
    cluster <- annotate_df(df)

    x$cluster <- cluster$cluster
    x$data$cluster <- cluster$cluster$cluster
    x$data <- dplyr::bind_cols(x$data, cluster$coords)
    x$data <- dplyr::select(x$data, .data$position, .data$wt, .data$cluster, dplyr::everything())
    x$annotated <- TRUE

  } else if (is.data.frame(x)) {
    df <- validate_combined_dms(x)[c("study", "gene", "position", "wt", amino_acids)]
    cluster <- annotate_df(df)
    x <- dplyr::bind_cols(x, cluster$coords)
    x <- dplyr::bind_cols(x, cluster$cluster)
    x <- dplyr::select(x, .data$study, .data$wt, .data$position, .data$wt, .data$cluster, dplyr::everything())

  } else {
    stop("Unrecognised input: x must be a deep_mutational_scan or data frame containing data from multiple scans")
  }

  return(x)
}

#' Assign PCA, UMAP and clusters to a data frame of positional ER scores
#'
#' Internal function called by \code{\link{annotate}} to perform the mechanics of annotations
#'
#' @param df Data frame
#' @return A list containing two data frames, the first covering PCA and UMAP assignment and the second
#' cluster assignment
#' @keywords internal
annotate_df <- function(df) {
  # Map onto PCAs
  pca <- stats::predict(dms_pca, newdata = as.matrix(df[amino_acids]))
  coords <- tibble::as_tibble(pca)

  # Map onto UMAP space
  model <- uwot::load_uwot(system.file("extdata", "dms_umap", package = "deepscanscape"))
  umap <- uwot::umap_transform(as.matrix(df[amino_acids]), model = model)
  coords[c("umap1", "umap2")] <- umap

  # Assign clusters
  distance <- calculate_cluster_distances(df$wt, pca[, -1])
  closest <- apply(distance, 1, which.min)

  cluster <- tibble::tibble(wt = df$wt, cluster = stringr::str_c(.data$wt, closest), base_cluster = .data$cluster)
  cluster$permissive <- permissive_positions(df[, amino_acids])
  cluster$ambiguous <- ambiguous_assignment(distance)
  cluster$high_distance <- outlier_positions(distance)

  cluster$cluster[cluster$ambiguous] <- stringr::str_c(cluster$wt[cluster$ambiguous], "A")
  cluster$cluster[cluster$high_distance] <- stringr::str_c(cluster$wt[cluster$high_distance], "O")
  cluster$cluster[cluster$permissive] <- stringr::str_c(cluster$wt[cluster$permissive], "P")

  cluster <- dplyr::bind_cols(df[c("position")], cluster, as.data.frame(distance))
  cluster$notes <- cluster_notes(cluster$permissive, cluster$ambiguous, cluster$high_distance)

  return(list(coords = coords, cluster = cluster))
}

#' Calculate distance to each cluster of the amino acid
#'
#' Calculate the distances to each cluster of that amino acid for each row of the input PCA matrix
#'
#' Internal function called by \code{\link{annotate_df}}.
#'
#' @param wt Wild type amino acid for each position
#' @param pca Numeric matrix of PCA values for PC2 to PC20
#' @return Matrix of cosine distances to the centroids of clusters 1 to 8 for the WT amino acid. When the wild type
#' amino acid has fewer than 8 clusters NA values occur in the additional columns
#' @keywords internal
calculate_cluster_distances <- function(wt, pca) {
  dists <- matrix(nrow = length(wt), ncol = 8)
  colnames(dists) <- stringr::str_c("dist", 1:8)

  for (aa in amino_acids) {
    c <- cluster_centers[stringr::str_starts(rownames(cluster_centers), aa), , drop = FALSE]
    dists[wt == aa, seq_len(nrow(c))] <- cosine_distance_matrix(pca[wt == aa, ], c)
  }

  return(dists)
}

#' Identify ambiguous cluster assignment
#'
#' Identify positions where the difference between the minimum distance to a cluster and the distance any other is
#' less than 0.03, which indicates ambiguous assignment.
#'
#' Internal function called by \code{\link{annotate_df}}.
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
#' Internal function called by \code{\link{annotate_df}}.
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
#' Internal function called by \code{\link{annotate_df}}.
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

#' Describe amino acid positional subtypes
#'
#' Add descriptions and frequencies from the cluster assigned to each position in an annotated deep_mutational_scan
#' dataset. These were determined through analysis of the larger \code{\link{deep_mutational_scans}} dataset. Numerical
#' summary statistics based on the same dataset can also be added, which give the mean characteristics of the subtype.
#' They include results from SIFT4G, FoldX, Naccess and mean ER fitness score profiles.
#'
#' @param x \link{deep_mutational_scan} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}
#' @param full Logical. Include average statistics from the \code{\link{deep_mutational_scans}} dataset in
#' addition to summary descriptions and notes.
#' @return A \code{\link[tibble]{tibble}}, with each row detailing a row of the input data and columns matching
#' those in the \code{\link{subtypes}} dataset.
#' @examples
#' dms <- deepscanscape::deep_scan
#'
#' # Basic description
#' describe_clusters(dms)
#'
#' # More details
#' describe_clusters(dms, full = TRUE)
#'
#' @export
describe_clusters <- function(x, full = FALSE) {
  if (is.deep_mutational_scan(x)) {
    if (!x$annotated) {
      warning("deep_mutational_scan is not annotated. Annotating using annotate_dms().")
      x <- annotate(x)
    }

    df <- x$data
    df$gene <- x$gene
    df$study <- x$study
  } else if (is.data.frame(x)) {
    df <- validate_combined_dms(x, annotated = TRUE)
  } else {
    stop("Unrecognised input: x must be a deep_mutational_scan or a data frame")
  }

  extra <- dplyr::rename(deepscanscape::subtypes, global_cluster_freq = .data$prop)
  if (!full) {
    extra <- extra[c("wt", "cluster", "global_cluster_freq", "group", "description", "notes")]
  }

  df <- dplyr::left_join(df[c("study", "gene", "position", "wt", "cluster")], extra, by = c("wt", "cluster"))
  return(df)
}
