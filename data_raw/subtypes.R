# Generate subtype descriptions object
# This is based on the supplementary data from
# Dunham & Beltrao (2021) - https://www.embopress.org/doi/full/10.15252/msb.202110305
library(readr)
library(dplyr)
library(stringr)

dms <- read_tsv("data_raw/combined_mutational_scans.tsv")

prop <- group_by(dms, wt, cluster) %>% summarise(n = n()) %>% mutate(p = n / sum(n)) %>% select(-n)

features <- group_by(dms, wt, cluster) %>%
  summarise(across(c(A:Y, mean_score, mean_sift, total_energy:energy_ionisation, all_atom_rel), mean, na.rm = TRUE),.groups = "drop") %>%
  mutate(across(c(A:Y, mean_score, mean_sift, total_energy:energy_ionisation, all_atom_rel), ~ifelse(str_ends(cluster, "O"), NA, .x)))

subtypes <- read_tsv("data_raw/subtype_descriptions.tsv") %>%
  left_join(prop, by = c("wt", "cluster")) %>%
  select(wt, cluster, prop = p, group, description, notes) %>%
  left_join(features, by = c("wt", "cluster"))

usethis::use_data(subtypes, overwrite = TRUE, compress = 'xz', version = 3)
