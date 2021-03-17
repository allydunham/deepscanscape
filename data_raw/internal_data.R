# Generate internal package datasets
library(tidyverse)

# PCA
dms_pca <- readRDS('data_raw/dms_pca.rds')

# Amino Acids
amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

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

# FoldX column names
foldx_terms <- c("total_energy", "backbone_hbond", "sidechain_hbond", "van_der_waals", "electrostatics",
                 "solvation_polar", "solvation_hydrophobic", "van_der_waals_clashes", "entropy_sidechain",
                 "entropy_mainchain", "cis_bond", "torsional_clash", "backbone_clash", "helix_dipole",
                 "disulfide", "electrostatic_kon", "partial_covalent_bonds", "energy_ionisation")

# Median DMS Scores
median_scores <- read_tsv('data_raw/raw_mutational_scans.tsv') %>%
  filter(!mut == '*', !is.na(score)) %>%
  group_by(wt, mut) %>%
  summarise(score = median(score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = 'mut', values_from = 'score') %>%
  tblhelpr::tibble_to_matrix(-wt, row_names = 'wt')

# Cluster centroids
cluster_centers <- read_tsv('data_raw/combined_mutational_scans.tsv') %>%
  filter(!str_detect(cluster, "[A-Y]O"), !str_detect(cluster, "[A-Y]P")) %>%
  group_by(cluster) %>%
  summarise(across(.cols = c(PC1:PC20), .fns = mean))
cluster_centers <- select(cluster_centers, PC2:PC20) %>%
  as.matrix() %>%
  magrittr::set_rownames(cluster_centers$cluster)

usethis::use_data(dms_pca, median_scores, cluster_centers, aa_3_to_1, amino_acids, foldx_terms,
                  internal = TRUE, overwrite = TRUE, compress = "xz", version = 3)
