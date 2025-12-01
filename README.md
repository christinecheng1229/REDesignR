
<!-- README.md is generated from README.Rmd. Please edit that file -->

# REDesignR

<!-- badges: start -->

<!-- badges: end -->

## Description

The goal of REDesignR is to help scientists determine the best
restriction enzyme(s) to use for their sequence digestion experiment.
The tool will be an R package with a Shiny interface that simulates,
optimizes and visualizes restriction enzyme digestion experiments. While
web interfaces like REBASE (Roberts, 2015) exist to help with optimal
restriction enzyme choice and packages like DECIPHER (Wright, 2024) can
simulate the fragments that results from digesting a given sequence with
the specified restriction enzyme, these existing tools lack in
automating optimal enzyme choice, providing intuitive
visualization/interactive interfaces, and simulating multi-enzyme
co-digestion (using multiple restriction enzymes to digest a single
sequence) experiments that many scientists choose to do. ‘REDesignR’ was
developed using R version ‘4.5.1’, platform ‘x86_64-w64-mingw32/x64’,
and running under ‘Windows 11 x64 (build 26100)’.

TODO: expand on novelty e.g., - Visualize digest results (and
compatibility matrices) in one workflow. - certain R packages(REDseq)
and web tools (RestrictionMapper) help visualize digests as restriction
maps, however none/not many visualize results as agarose gel, which can
help scientists validate their wet lab results with the simulation. This
is a critical step in ensuring the digestion performed as expected and
won’t confound downstream experiments/results.

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install github("christinecheng1229/REDesignR",build vignettes = TRUE)
library("REDesignR")
```

To run the shinyApp: TODO Under construction

``` r
# runREDesignR()
```

## Overview

TODO ls() doesn’t work as expected below

``` r
ls("package:REDesignR") # View package functions and datasets
data(package = "REDesignR") # View package dataset(s)
browseVignettes("REDesignR")  # View package vignette(s)
```

`REDesignR` contains 3 functions.

1.  ***simulateCoDigest*** for simulating a co-digestion experiment and
    producing a table of resulting digests.

2.  ***plotFragments*** for visualizing the distribution of fragment
    sizes after a digestion.

3.  ***optimalRE*** for determining the optimal combination of
    restriction enzyme to use for digestion, based on experimental
    environments and desired average fragment length.

The package also contains a dataset, called ‘Enzymes’ that can be used
as a source of restriction enzymes and their recognition sites.

<!-- An overview of the package is illustrated below.-->

<!-- ![](./inst/extdata/ExampleImage.png) -->

## Contributions

The author of the package is Christine Cheng. OpenAI ChatGPT-5.1 Auto
model was used during the package development process to aid in
maintaining consistent formatting across R sripts, identifying gaps in
unit test coverage, TODO. TODO Outline contributions from other
packages/sources for each function, etc. Remember your individual
contributions to the package are (also) important.

## References

- Müller K, Wickham H (2025). *tibble: Simple Data Frames*.
  <doi:10.32614/CRAN.package.tibble>
  <https://doi.org/10.32614/CRAN.package.tibble>, R package version
  3.3.0, <https://CRAN.R-project.org/package=tibble>.
- OpenAI. (2025). ChatGPT (GPT-5.1) \[Large language model\].
  <https://openai.com/chatgpt>
- Pagès H, Aboyoun P, Gentleman R, DebRoy S (2025). *Biostrings:
  Efficient manipulation of biological strings*.
  <doi:10.18129/B9.bioc.Biostrings>
  <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version
  2.78.0, <https://bioconductor.org/packages/Biostrings>.
- Pagès H, Lawrence M, Aboyoun P (2025). *S4Vectors: Foundation of
  vector-like and list-like containers in Bioconductor*.
  <doi:10.18129/B9.bioc.S4Vectors>
  <https://doi.org/10.18129/B9.bioc.S4Vectors>, R package version
  0.48.0, <https://bioconductor.org/packages/S4Vectors>.
- R Core Team (2025). *R: A Language and Environment for Statistical
  Computing*. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.
- Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. REBASE-a database
  for DNA restriction and modification: enzymes, genes and genomes.
  Nucleic Acids Res. 43: D298-D299 (2015).
- Wickham H (2019). *assertthat: Easy Pre and Post Assertions*.
  <doi:10.32614/CRAN.package.assertthat>
  <https://doi.org/10.32614/CRAN.package.assertthat>, R package version
  0.2.1, <https://CRAN.R-project.org/package=assertthat>.
- Wright ES (2024). “Fast and Flexible Search for Homologous Biological
  Sequences with DECIPHER v3.” *The R Journal*, *16*(2), 191-200.

## Acknowledgements

This package was developed as part of an assessment for 2025 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. REDesignR welcomes issues, enhancement requests, and other
contributions. To submit an issue,use the GitHub issues.

<!-- ## Package Tree Structure -->

<!-- TODO -->
