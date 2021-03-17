# Example deep mutational scan, primarily used for demonstrating/testing functions
# TODO Choose a better scan
library(devtools)
load_all()

csv <- read.csv('inst/extdata/urn_mavedb_00000036-a-1_scores.csv', skip = 4)
deep_scan <- deep_mutational_scan(csv, scheme = 'mave', trans = 'vamp', na_value = 'impute', annotate = FALSE,
                                  gene = 'LDLRAP1',  study = 'urn:mavedb:00000036')

usethis::use_data(deep_scan, overwrite = TRUE, compress = 'xz', version = 3)
