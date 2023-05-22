# DeepScanScape 1.0.0

Initial package release, accompanying the publication [Dunham & Beltrao (2021)](https://www.embopress.org/doi/full/10.15252/msb.202110305)

# DeepScanScape 1.0.1

- Fixed plain landscape plots

# DeepScanScape 1.1.0

- Added support for importing synonymous variants marked with "=", "sym" or a user defined list of symbols
- Bugfixed custom imputation matrices
- Added basic support for stop codons. Currently imputing early stop values isn't directly supported, although the user can easily apply their own algorithm. Frameshifts and indels are unsupported since they do not fit the package data model. Stop scores are not used in many downstream applications and are mainly made available for use by the user, although they are plotted on the summary heatmap.