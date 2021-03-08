# S3 class / attributes thing to store a dms dataset
# TODO - Test these funcs

#' Deep mutational scan data
#'
#' Store deep mutational scanning data in a standardised format, alongside
#' metadata describing the state of the data. The primary data is contained in
#' a tibble
#'
#' @param df A data frame to parse
#' @param scheme Original data scheme (see \code{\link{parse_dms_data}})
#' @param trans Function to transform scores onto the standard scale. Accepts a string corresponding to known transforms
#' or a custom function (\code{\link{transform_dms}}).
#' @param na_value Value to set missense NA scores to (see \code{\link{impute_dms}})
#' @param annotate Annotate the dataset with mutational landscape data (PCA, UMAP and amino acid subtypes)
#' @param study Source of the deep mutational scan
#' @param gene Gene scanned
#' @return A deep mutational scan S3 object
#' @export
deep_mutational_scan <- function(df, scheme=NULL, trans=NULL, na_value="impute",
                                 annotate=TRUE, study=NULL, gene=NULL) {

  # Initialise object
  out <- list(study = study, gene = gene, annotated = FALSE)
  class(out) <- c("deep_mutational_scan", class(out))

  # Parse scheme
  if (!is.null(scheme)) {
    df <- parse_dms_data(df, scheme)
  } else {
    df <- tibble::as_tibble(df)
    if (!all(c("position", "wt", "mut", "score") %in% names(df))) {
      stop("Incorrect input df, must include position, wt, mut and score columns")
    }
    df <- dplyr::select(df, .data$position, .data$wt, .data$mut, .data$score, dplyr::everything())
  }

  # Transform
  if (!is.null(trans)) {
    df$score <- transform_dms(df$score, trans)
  }

  # Normalise
  df$score <- normalise_dms(df$score)

  # Format tibble
  df <- df[df$mut %in% amino_acids, ]
  df <- dplyr::arrange(df, .data$position, .data$mut)
  df <- tidyr::pivot_wider(df, names_from = "mut", values_from = "score")
  df <- dplyr::select(df, .data$position, .data$wt, amino_acids, dplyr::everything())

  # Set attributes
  out$data <- df

  # Impute
  if (!(is.na(na_value) | is.null(na_value))) {
    out <- impute_dms(out)
  }

  # Annotate
  if (annotate) {
    out <- annotate_dms(out)
  }

  return(out)
}

# TODO Document this fully - [ and [[ go to data tbl, $ goes to list
#' Extracting and replacing deep mutational scanning data
#'
#' @param x \link{deep_mutational_scan} object
#' @param i,j,... Indeces to access
#' @param drop Coerce result to lowest possible dimension
#' @param exact Ignored for tibbles (which is the underlying data structure)
#' @param value Value to set
#'
#' @name dms_extract
NULL
#> NULL

#' @describeIn dms_extract Extract
#' @export
`[.deep_mutational_scan` <- function(x, i, j, drop = FALSE, ...) {  # nolint
  if (missing(j)) {
    return(x$data[i])
  }
  return(x$data[i, j, ..., drop = drop])
}

#' @describeIn dms_extract Assign
#' @export
`[<-.deep_mutational_scan` <- function(x, i, j, ..., value) {  # nolint
  if (missing(j)) {
    x$data[i] <- value
  } else {
    x$data[i, j, ...] <- value
  }
  return(x)
}

#' @describeIn dms_extract Double Extract
#' @export
`[[.deep_mutational_scan` <- function(x, i, j, exact = FALSE, ...) {  # nolint
  if (missing(j)) {
    return(x$data[[i]])
  }
  return(x$data[[i, j, ..., exact = exact]])
}

#' @describeIn dms_extract Double Assign
#' @export
`[[<-.deep_mutational_scan` <- function(x, i, j, ..., value) {  # nolint
  if (missing(j)) {
    x$data[[i]] <- value
  } else {
    x$data[[i, j, ...]] <- value
  }
  return(x)
}

# TODO - implement other generics (length, dim, etc.)
#' General S3 Methods for Deep Mutational Scans
#'
#' @param x,object \link{deep_mutational_scan} object
#' @param ... Additional arguments
#'
#' @name dms_s3
NULL
#> NULL

# TODO - Pretty print when extra columns added
#' @describeIn dms_s3 S3 print method
#' @export
format.deep_mutational_scan <- function(x, ...) {
  out <- c(paste("# A deep_mutational_scan"),
           paste("# Study:", ifelse(is.null(x$study), "Unknown", x$study)),
           paste("# Gene:", ifelse(is.null(x$gene), "Unknown", x$gene)),
           paste("#", nrow(x$data), "positions"),
           "# Positional data:",
           format(x$data)[-1])

  if (requireNamespace("crayon", quietly = TRUE)) {
    out[1:5] <- crayon::make_style("darkgrey")(out[1:5])
  }

  return(out)
}

#' @describeIn dms_s3 S3 print method
#' @export
print.deep_mutational_scan <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
}

#' @describeIn dms_s3 S3 str method
#' @export
str.deep_mutational_scan <- function(object, ...) {
  stop("Not implemented yet")
}


#' @describeIn dms_s3 S3 as_tibble method
#' @importFrom tibble as_tibble
#' @export
as_tibble.deep_mutational_scan <- function(x, ...) { # nolint
  return(x$data)
}

#' @describeIn dms_s3  S3 as.data.frame method
#' @export
as.data.frame.deep_mutational_scan <- function(x, ...) {
  return(as.data.frame(as_tibble(x, ...)))
}

# Create tibble with study and gene columns
to_full_tibble <- function(x) {
  out <- tibble::as_tibble(x)
  out$study <- x$study
  out$gene <- x$gene
  return(dplyr::select(out, .data$study, .data$gene, dplyr::everything()))
}

#' Combine deep mutational scan data
#'
#' @param ... \link{deep_mutational_scan} objects to combine
#'
#' @export
rbind.deep_mutational_scan <- function(...) {
  dplyr::bind_rows(lapply(list(...), to_full_tibble))
}

# TODO - Need output object for this?
#' Summarise Deep Mutational Scans
#'
#' @param object A \link{deep_mutational_scan} object.
#' @param ... Additional arguments
#'
#' @export
summary.deep_mutational_scan <- function(object, ...) {
  stop("Not implemented yet")
}

#' Summary plot for Deep Mutational Scans
#'
#' @param object,x \link{deep_mutational_scan} to produce a
#' summary plot for.
#' @param ... Additional parameters
#'
#' @export
#' @name dms_plot
#' @importFrom ggplot2 autoplot
autoplot.deep_mutational_scan <- function(object, ...){ # nolint
  stop("Not implemented yet")
}

#' @describeIn dms_plot S3 plot method
#' @export
#' @importFrom graphics plot
plot.deep_mutational_scan <- function(x, ...) {
  print(autoplot(x, ...))
}
