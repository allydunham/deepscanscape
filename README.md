# DeepScanScape - Annotate deep mutational scans using the deep landscape
<!-- badges: start -->
[![R-CMD-check](https://github.com/allydunham/deepscanscape/workflows/R-CMD-check/badge.svg)](https://github.com/allydunham/deepscanscape/actions)
<!-- badges: end -->

Process and annotate deep mutational scanning data using a deep mutational landscape derived from 28 previous studies [(Dunham & Beltrao 2021)](https://doi.org/10.15252/msb.202110305), highlighting unusual properties, visualising data and assigning positions to amino acid subtypes that predict their likely features.

## Installation

Install the latest version from GitHub 

```R
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("allydunham/deepscanscape")
```

## Usage

We recommend new users start with the package [vignette](https://allydunham.github.io/deepscanscape/articles/deepscanscape.html), which shows an example analysis importing and assessing three typical deep scans.
The examples in the [documentation](https://allydunham.github.io/deepscanscape/reference/deepscanscape.html) also illustrate the usage of each function.
Reading the publication describing the theory behind the deep landscape may also help explain package usage.

## Citation

Alistair Dunham & Pedro Beltrao (2020). Exploring amino acid functions in a deep mutational landscape. bioRxiv [2020.05.26.116756](https://www.biorxiv.org/content/10.1101/2020.05.26.116756v1); doi:
<https://doi.org/10.1101/2020.05.26.116756>
