# Purpose: Tests for simulateGel function
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.0
# Bugs and Issues: None known.

test_that("simulateGel fails on mis-formatted input data", {
  testDigestDF <- simulateCoDigest(Biostrings::DNAString("ATCGATGGATCCATCGATATCGATGGATCC"),
                                   Enzymes[c(4, 27),])
  expect_error(REDesignR::simulateGel(t(testDigestDF)))
})

test_that("simulateGel rejects incomplete digestion data", {
  badInput <- tibble::tibble(A = 1, B = 2)
  expect_error(simulateGel(badInput), "valid digestion dataframe")
})

test_that("simulateGel warns when multiDigest = TRUE but data is single-enzyme", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")

  single <- simulateCoDigest(dna, Enzymes[c(1, 6), ])

  expect_warning(simulateGel(single$digestDf))
})

test_that("simulateGel warns when multiDigest = FALSE but data is co-digest", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")

  co <- simulateCoDigest(dna, Enzymes[c(1, 4), ])

  expect_warning(simulateGel(co$digestDf, multiDigest = FALSE))
})

test_that("simulateGel returns a ggplot for single-digest", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")

  single <- simulateCoDigest(dna, Enzymes[c(1, 6), ])

  testPlot <- simulateGel(single$digestDf, multiDigest = FALSE)

  expect_s3_class(testPlot, "gg")
})

test_that("simulateGel returns a ggplot for co-digest", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")

  co <- simulateCoDigest(dna, Enzymes[c(1, 4), ])

  testPlot <- simulateGel(co$digestDf, multiDigest = TRUE)

  expect_s3_class(testPlot, "gg")
})

# [END]
