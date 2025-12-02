# Purpose: Tests for plotRestrictionMap function
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.0
# Bugs and Issues: None known.

test_that("plotRestrictionMap fails on mis-formatted input data", {
  testDigestDF <- simulateCoDigest(Biostrings::DNAString("ATCGATGGATCCATCGATATCGATGGATCC"),
                                   Enzymes[c(4, 27),])
  expect_error(REDesignR::plotRestrictionMap(t(testDigestDF)))
})

test_that("plotRestrictionMap rejects non-tibble/input missing columns", {
  badInput <- tibble::tibble(A = 1)
  expect_error(plotRestrictionMap(badInput), "valid digestion dataframe")
})

test_that("plotRestrictionMap returns a plot for single-digest input", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")
simulateCoDigest(Biostrings::DNAString("CGGCCGATCGATCGGCCG"), Enzymes[c(1, 6), ])
  single <- simulateCoDigest(dna, Enzymes[c(1, 6), ])

  testPlot <- plotRestrictionMap(single$digestDf, multiDigest = FALSE)

  expect_s3_class(testPlot, "gg")
})

test_that("plotRestrictionMap returns patchwork for co-digest", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")

  co <- simulateCoDigest(dna, Enzymes[c(1, 4), ])

  testPlot <- plotRestrictionMap(co$digestDf, multiDigest = TRUE)

  # patchwork objects inherit class "patchwork"
  expect_true("patchwork" %in% class(testPlot))
})

# [END]
