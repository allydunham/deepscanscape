# Generate subtype descriptions object
# This is based on the supplementary data from CITATION
# URL
# TODO - Add proper citation here
library(readr)
library(dplyr)

dms <- read_tsv("data_raw/combined_mutational_scans.tsv")

prop <- group_by(dms, wt, cluster) %>% summarise(n = n()) %>% mutate(p = n / sum(n)) %>% select(-n)

subtypes <- read_tsv("data_raw/subtype_descriptions.tsv") %>%
  left_join(prop, by = c("wt", "cluster")) %>%
  select(wt, cluster, prop = p, group, description, notes)

usethis::use_data(subtypes, overwrite = TRUE, compress = 'xz', version = 3)
