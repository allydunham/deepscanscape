# Parse common input formats in preparation for analysis
# TODO - Test funcs
# TODO - Add additional formats

aa_3_to_1 <- c(ala = "A", arg = "R", asn = "N", asp = "D", cys = "C", gln = "Q",
               glu = "E", gly = "G", his = "H", ile = "I", leu = "L", lys = "K",
               met = "M", phe = "F", pro = "P", ser = "S", thr = "T", trp = "W",
               tyr = "Y", val = "V", sec = "U", pyl = "O", asx = "B", xle = "J",
               glx = "Z", xaa = "X", ter = "*")

# aa_1_to_3 <- c(A = "ala", R = "arg", N = "asn", D = "asp", C = "cys", Q = "gln",
#                E = "glu", G = "gly", H = "his", I = "ile", L = "leu", K = "lys",
#                M = "met", F = "phe", P = "pro", S = "ser", T = "thr", W = "trp",
#                Y = "tyr", V = "val", U = "sec", O = "pyl", B = "asx", J = "xle",
#                Z = "glx", X = "xaa", `*` = "ter")

#' Parse deep mutational scanning data
#'
#' @param x A tibble or an object that can be converted to one
#' @param scheme Original data scheme (see description). Accepts any unambiguous
#' substring.
#' @param position_offset Offset all positions by this amount
#' @param score_col Name of column containing fitness scores (used for when
#' parsing MaveDB data)
#'
#' @export
parse_dms_data <- function(x, scheme, position_offset=0, score_col=NULL) {
  methods <- c("mavedb" = parse_mavedb,
               "long" = parse_long,
               "wide" = parse_wide,
               "sequence" = parse_sequence)

  scheme <- names(methods)[pmatch(stringr::str_to_lower(scheme), names(methods))]
  if (is.na(scheme)) {
    stop(paste0("Unrecognised scheme \"", scheme, "\"\n",
                "Recognised options: ", paste(methods, collapse = ", ")))
  }

  x <- tibble::as_tibble(x)
  func <- methods[[scheme]]
  x <- func(x, score_col)
  x$position <- x$position + position_offset
  return(x)
}

# Parse data in the MaveDB format (called from parse_dms_data)
parse_mavedb <- function(x, score_col, ...) {
  if (is.null(score_col)) {
    score_col <- "score"
  }
  x <- x[c("hgvs_pro", score_col)]
  x[c("wt", "position", "mut")] <- stringr::str_match(x$hgvs_pro, "p\\.([A-Za-z=\\*]*)([0-9]*)([A-Za-z=\\*]*)")[, 2:4]
  x$wt <- unname(aa_3_to_1[stringr::str_to_lower(x$wt)])
  x$mut <- unname(aa_3_to_1[stringr::str_to_lower(x$mut)])
  x$position <- as.numeric(x$position)
  return(x[c("position", "wt", "mut", "score")])
}

# Parse data in the long format (called from parse_dms_data)
parse_long <- function(x, ...) {
  stop("Not implemented yet")
}

# Parse data in the wide format (called from parse_dms_data)
parse_wide <- function(x, ...) {
  stop("Not implemented yet")
}

# Parse data in the sequence format (called from parse_dms_data)
parse_sequence <- function(x, ...) {
  stop("Not implemented yet")
}
