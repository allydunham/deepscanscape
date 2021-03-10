# Annotate positions based on the overall dataset
# TODO test annotation functions

# TODO Store distance measurements in class
# TODO Note secondary clusters when margin is very small?
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
  closest <- apply(cluster_dists, 1, which.min)
  cluster <- stringr::str_c(x$data$wt, closest)

  # Permissive positions
  permissive <- apply(x$data[, amino_acids], 1, check_permissive)
  cluster[permissive] <- stringr::str_c(x$data$wt[permissive], "P")

  x$data$cluster <- cluster

  # Identify low confidence  positions
  margin <- t(apply(cluster_dists, 1, calculate_margin))
  margin_small <- apply(margin, 1, any_less_than, threshold = 0.03)

  check_distance <- function(x) {
    min(x, na.rm = TRUE) > 0.45
  }
  dist_high <- apply(cluster_dists, 1, check_distance)

  notes <- rep(NA, length(cluster))
  # TODO Note which clusters are ambiguous
  # TODO better naming here?
  notes[!permissive & margin_small] <- "margin"
  notes[!permissive & dist_high] <- "distance"
  notes[!permissive & margin_small & dist_high] <- "both"
  x$data$cluster_notes <- notes

  x$data <- dplyr::select(x$data, .data$position, .data$wt, .data$cluster, .data$cluster_notes, dplyr::everything())
  x$annotated <- TRUE
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
