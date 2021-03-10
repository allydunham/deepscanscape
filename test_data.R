# Reload script for interactive testing
library(devtools)
document()
load_all()

csv <- read.csv('inst/extdata/urn_mavedb_00000036-a-1_scores.csv', skip = 4)
dms <- deep_mutational_scan(csv, scheme = 'mave', trans = 'vamp', annotate = TRUE, gene = 'LDLRAP1',
                            study = 'urn:mavedb:00000036', na_value = 'impute')
