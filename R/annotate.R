# Annotate positions based on the overall dataset
# TODO test annotation functions

# TODO add more parameters
#' Annotate positions based on the deep mutational landscape
#'
#' Add PCA and UMAP transformed coordinates and amino acid subtypes based on analysis and clustering of the combined
#' \link[=deep_mutational_scans]{deep mutational landscape dataset}.
#'
#' PCA and UMAP coordinates are assigned using the original models fit on the
#' \link[=deep_mutational_scans]{deep mutational landscape dataset}.
#'
#' @param x \code{\link{deep_mutational_scan}} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}
#' to annotate
#' @return An annotated \code{\link{deep_mutational_scan}} if passed a \code{\link{deep_mutational_scan}} or a
#' \code{\link[tibble]{tibble}} when passed a \link[=rbind.deep_mutational_scan]{combined DMS data frame}.
#' This data frame contains the following added columns:
#' \itemize{
#'   \item cluster: the assigned amino acid subtype
#'   \item PC1 - PC20: Principal Component coordinates
#'   \item umap1/2: UMAP coordinates
#'   \item base_cluster: The nearest primary cluster centroid (i.e. not outlier or permissive clusters)
#'   \item permissive: The position is identified as permissive (|ER| < 0.4 for all amino acids)
#'   \item small_margin: The distance to two clusters is very similar, so the assignment is low confidence and marked as
#'   an outlier
#'   \item high_distance: The position is distant from all cluster centroids and is marked an outlier.
#'   \item dist1-8: The distance to each cluster of the WT amino acid
#'   \item notes: Notes on the cluster assignment
#' }
#' @export
annotate_dms <- function(x) {
  if (is.deep_mutational_scan(x)) {
    return(annotate_dms_deep_scan(x))
  } else if (is.data.frame(x)) {
    x <- tibble::as_tibble(x)
    return(annotate_dms_tibble(x))
  } else {
    stop("Unrecognised input: x must be a deep_mutational_scan or data frame containing data from multiple scans")
  }
}

#' Annotate a Deep Mutational Scan object
#'
#' Internal function called by \code{\link{annotate_dms}} when passed a \code{\link{deep_mutational_scan}}
#'
#' @param x \code{\link{deep_mutational_scan}}
annotate_dms_deep_scan <- function(x) {
  if (!is.deep_mutational_scan(x)) {
    stop("Unrecognised data.\nCreate a standardised dataset using deep_mutational_scan()")
  }

  if (any(is.na(x[amino_acids]))) {
    stop("NA fitness scores present\nUse impute_dms to remove these")
  }

  if (x$annotated) {
    stop("dataset already annotated")
  }

  # Map onto PCAs
  pca <- stats::predict(dms_pca, newdata = as.matrix(x[amino_acids]))
  x$data <- dplyr::bind_cols(x$data, tibble::as_tibble(pca))

  # Map onto UMAP space
  model <- uwot::load_uwot(system.file("extdata", "dms_umap", package = "deepscanscape"))
  umap <- uwot::umap_transform(as.matrix(x[amino_acids]), model = model)
  x[c("umap1", "umap2")] <- umap

  # Assign subtypes
  calc_dist <- function(i) {
    calculate_cluster_distances(x$data$wt[i], pca[i, -1])
  }

  cluster_dists <- t(sapply(seq_len(nrow(x$data)), calc_dist))
  colnames(cluster_dists) <- stringr::str_c("dist", 1:8)
  closest <- apply(cluster_dists, 1, which.min)
  cluster <- tibble::tibble(wt = x$data$wt, cluster = stringr::str_c(.data$wt, closest), base_cluster = .data$cluster)

  # Permissive positions
  cluster$permissive <- apply(x$data[, amino_acids], 1, check_permissive)

  # Identify low confidence  positions
  margin <- t(apply(cluster_dists, 1, calculate_margin))
  cluster$small_margin <- apply(margin, 1, any_less_than, threshold = 0.03)

  check_distance <- function(x) {
    min(x, na.rm = TRUE) > 0.45
  }
  cluster$high_distance <- apply(cluster_dists, 1, check_distance)

  # Resolve cluster assignment
  outlier_ind <- cluster$small_margin | cluster$high_distance
  cluster$cluster[outlier_ind] <- stringr::str_c(cluster$wt[outlier_ind], "O")
  cluster$cluster[cluster$permissive] <- stringr::str_c(cluster$wt[cluster$permissive], "P")
  cluster <- dplyr::bind_cols(cluster, as.data.frame(cluster_dists))
  cluster <- dplyr::select(cluster, -.data$wt)

  # Add cluster notes
  cluster$notes <- rep(NA, nrow(cluster))
  cluster$notes[cluster$small_margin] <- "Top cluster is only nearer by a small margin"
  cluster$notes[cluster$high_distance] <- "Not close to any cluster center"
  cluster$notes[cluster$small_margin & cluster$high_distance] <- "Both distant and only nearer one by a small margin"
  cluster$notes[cluster$permissive] <- "No mutation with |ER| > 0.4"

  # Tidy object
  x$cluster <- cluster
  x$data$cluster <- cluster$cluster
  x$data <- dplyr::select(x$data, .data$cluster, .data$position, .data$wt, dplyr::everything())
  x$annotated <- TRUE
  return(x)
}

#' Annotate a Deep Mutational Scan object
#'
#' Internal function called by \code{\link{annotate_dms}} when passed a \code{\link[tibble]{tibble}}
#'
#' @param x \code{\link[tibble]{tibble}}
annotate_dms_tibble <- function(x) {
  x <- validate_combined_dms(x)[c("study", "gene", "position", "wt", amino_acids)]

  # Map onto PCAs
  pca <- stats::predict(dms_pca, newdata = as.matrix(x[amino_acids]))
  x <- dplyr::bind_cols(x, tibble::as_tibble(pca))

  # Map onto UMAP space
  model <- uwot::load_uwot(system.file("extdata", "dms_umap", package = "deepscanscape"))
  umap <- uwot::umap_transform(as.matrix(x[amino_acids]), model = model)
  x[c("umap1", "umap2")] <- umap

  # Assign subtypes
  calc_dist <- function(i) {
    calculate_cluster_distances(x$wt[i], pca[i, -1])
  }

  cluster_dists <- t(sapply(seq_len(nrow(x)), calc_dist))
  colnames(cluster_dists) <- stringr::str_c("dist", 1:8)
  closest <- apply(cluster_dists, 1, which.min)
  cluster <- tibble::tibble(wt = x$wt, cluster = stringr::str_c(.data$wt, closest), base_cluster = .data$cluster)

  # Permissive positions
  cluster$permissive <- apply(x[, amino_acids], 1, check_permissive)

  # Identify low confidence  positions
  margin <- t(apply(cluster_dists, 1, calculate_margin))
  cluster$small_margin <- apply(margin, 1, any_less_than, threshold = 0.03)

  check_distance <- function(x) {
    min(x, na.rm = TRUE) > 0.45
  }
  cluster$high_distance <- apply(cluster_dists, 1, check_distance)

  # Resolve cluster assignment
  outlier_ind <- cluster$small_margin | cluster$high_distance
  cluster$cluster[outlier_ind] <- stringr::str_c(cluster$wt[outlier_ind], "O")
  cluster$cluster[cluster$permissive] <- stringr::str_c(cluster$wt[cluster$permissive], "P")
  cluster <- dplyr::bind_cols(cluster, as.data.frame(cluster_dists))
  cluster <- dplyr::select(cluster, -.data$wt)

  # Add cluster notes
  cluster$notes <- rep(NA, nrow(cluster))
  cluster$notes[cluster$small_margin] <- "Top cluster is only nearer by a small margin"
  cluster$notes[cluster$high_distance] <- "Not close to any cluster center"
  cluster$notes[cluster$small_margin & cluster$high_distance] <- "Both distant and only nearer one by a small margin"
  cluster$notes[cluster$permissive] <- "No mutation with |ER| > 0.4"

  # Tidy table
  x <- dplyr::bind_cols(x, cluster)
  x <- dplyr::select(x, .data$study, .data$gene, .data$position, .data$wt, .data$cluster, dplyr::everything())
  return(x)
}

#' Calculate distance to each cluster of the amino acid
#'
#' @param wt Wild type amino acid
#' @param pca PCA values for PC2 to PC20
calculate_cluster_distances <- function(wt, pca) {
  clus <- cluster_centers[stringr::str_starts(rownames(cluster_centers), wt), , drop = FALSE]
  m <- matrix(pca, nrow = nrow(clus), ncol = 19, byrow = TRUE)
  d <- cosine_distance(m, clus)
  names(d) <- NULL
  length(d) <- 8 # Maximum number of clusters for one AA in the dataset
  return(d)
}

#' Calculate Cosine distance between rows of two matrices
#'
#' @param x,y Numeric matrices
cosine_distance <- function(x, y) {
  return(acos(rowSums(x * y) / (sqrt(rowSums(x^2) * rowSums(y^2)))) / pi)
}

#' Check if ER scores qualifies as permissive
#'
#' @param x Vector of ER scores
check_permissive <- function(x) {
  return(all(abs(x) < 0.4))
}

#' Calculate the margin by which each element is largest than the smallest
#'
#' @param x Numeric vector
calculate_margin <- function(x) {
  i <- which.min(x)
  return(x - x[i])
}

#' Check if any non-zero values are below a threshold
#'
#' @param x Numeric vector
#' @param threshold Threshold below which TRUE is returned
any_less_than <- function(x, threshold) {
  if (any(x > 0, na.rm = TRUE)) {
    return(min(x[x > 0], na.rm = TRUE) < threshold)
  } else {
    return(FALSE)
  }
}

#' Describe amino acid positional subtypes
#'
#' @param x \link{deep_mutational_scan} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}
describe_clusters <- function(x) {
  if (is.deep_mutational_scan(x)) {
    if (!x$annotated) {
      warning("deep_mutational_scan is not annotated. Annotating using annotate_dms().")
      x <- annotate_dms(x)
    }

    df <- x$data
    df$gene <- x$gene
    df$study <- x$study
  } else if (is.data.frame(x)) {
    df <- validate_combined_dms(x, annotated = TRUE)
  } else {
    stop("Unrecognised input: x must be a deep_mutational_scan or a data frame")
  }
  df <- dplyr::left_join(df[c("study", "gene", "position", "wt", "cluster")],
                         dplyr::rename(deepscanscape::subtypes, global_cluster_freq = .data$prop),
                         by = c("wt", "cluster"))
  return(df)
}
