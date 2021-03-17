#' deepscanscape: Check and Annotate Deep Mutational Scanning Data
#'
#' The deepscanscape package provides functions to standardise deep
#' mutational scanning data and compare it to a deep mutational landscape
#' calculated from a dataset of 28 previous studies. This allows the whole
#' dataset to be checked for unusual properties and individual positions
#' to be annotated with positional subtypes, indicating the positions in other
#' proteins they are most similar to.
#'
#' @section Data Processing Functions:
#' \itemize{
#'   \item \code{\link{deep_mutational_scan}} - Construct standardised deep mutational scan datasets. This creates an
#'   S3 object with various generics.
#'   \item \code{\link{rbind.deep_mutational_scan}} - Combined deep mutational scan datasets
#'   \item \code{\link{parse_deep_scan}} - Parse common deep scan data formats
#'   \item \code{\link{transform_er}} - Transform ER scores from common score types to a standardised scale
#'   \item \code{\link{normalise_er}} - Normalise ER scores
#'   \item \code{\link{impute}} - Impute missing data from deep mutational scans
#' }
#'
#' @section Annotation Functions:
#' \itemize{
#'   \item \code{\link{annotate}} - Add annotations from the combined landscape to deep mutational scan data
#'   \item \code{\link{describe_clusters}} - Add details on positions assigned clusters
#'   \item \code{\link{recluster}} - Perform the original clustering procedure on a new deep mutational scan dataset
#' }
#'
#' @section Visualisation Functions:
#' \itemize{
#'   \item \code{\link{plot_er_heatmap}} - Plot heatmaps show fitness scores across a protein
#'   \item \code{\link{plot_landscape}} - Project a new dataset onto the deep mutational landscape, including
#'   visualising various biophysical properties.
#'   \item \code{\link{plot_cluster_frequencies}} - Plot the frequencies of amino acid subtypes in a new dataset
#'   \item \code{\link{plot_recluster}} - Summarise the profiles of a new clustered dataset
#' }
#'
#' @importFrom rlang .data
#' @docType package
#' @name deepscanscape
NULL
#>NULL

# TODO - Add examples to all functions
# TODO - Add vignettes
# TODO - Where to process counts? if at all
# TODO - Use generics internally?
# TODO - Function names - dms bit?
# TODO - Add to documentation - values, description, details, authour etc.
# TODO - Add progress bars to longer opperations?
# TODO - consistency of cluster/subtype
# TODO - better example datasets - create a list of several DMS datasets
# TODO - update citations over whole package
# TODO - note data frame much more fragile - possible to abuse functions by constructing your own
