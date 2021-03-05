# Generate deep mutational scan data object
# This is based on the supplementary data from CITATION
# URL
# TODO - Add proper citation here
library(readr)
library(dplyr)
x <- read_tsv('data_raw/combined_mutational_scans.tsv')

x <- select(x, -starts_with('log10_sift'), -all_atom_abs, -side_chain_abs, -side_chain_rel,
            -backbone_abs, -backbone_rel, -non_polar_abs, -non_polar_rel, -polar_abs,
            -polar_rel, -starts_with('within_10_0'), -starts_with('angstroms_to_'),
            -starts_with('ss_'))

deep_mutational_scans <- x
usethis::use_data(deep_mutational_scans, overwrite = TRUE, compress = 'xz', version = 3)
