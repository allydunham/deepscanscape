# Example deep mutational scan, primarily used for demonstrating/testing functions
library(devtools)
load_all()

deep_scans <- list()

# HSP90 - https://www.mavedb.org/scoreset/urn:mavedb:00000011-a-1/
csv <- read.csv('inst/extdata/urn_mavedb_00000011_a_1_scores.csv', skip = 4)
deep_scans$hsp90 <- deep_mutational_scan(csv, scheme = 'mave', trans = NULL, na_value = 'impute',
                                         annotate = FALSE, gene = 'Hsp90',  study = 'Hietpas et al. (2011)')

# GpA - https://www.mavedb.org/scoreset/urn:mavedb:00000051-c-1/
csv <- read.csv('inst/extdata/urn_mavedb_00000051_c_1_scores.csv', skip = 4)
df <- parse_deep_scan(csv, scheme = "mavedb", score_col = "ratio")
deep_scans$gpa <- deep_mutational_scan(df, scheme = NULL, trans = "log2", na_value = 'impute',
                                       annotate = FALSE, gene = 'GpA',  study = 'Elazar et al. (2016)')

# p53 - https://www.mavedb.org/scoreset/urn:mavedb:00000059-a-1/
csv <- read.csv('inst/extdata/urn_mavedb_00000059_a_1_scores.csv', skip = 4)
deep_scans$p53 <- deep_mutational_scan(csv, scheme = "mave", trans = NULL, na_value = 'impute',
                                       annotate = FALSE, gene = 'p53',  study = 'Kotler et al. (2018)')

usethis::use_data(deep_scans, overwrite = TRUE, compress = 'xz', version = 3)
