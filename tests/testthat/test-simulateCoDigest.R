# Purpose: Tests for simulateCoDigest function
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.0
# Bugs and Issues: None known.

test_that("simulateCoDigest works as intended with valid inputs", {
  expect_no_error(
    testDigestDF <- simulateCoDigest(
      Biostrings::DNAString("CGGCCGATCGATCGGCCG"),
      Enzymes[c(1,4),])
    )
  expectedResult <- tibble::tibble(Enzymes = rep(c("AaaI + AagI"), times = 4),
                                   FragmentID = c(1, 2, 3, 4),
                                   Start = c(1, 2, 9, 14),
                                   End = c(1, 8, 13, 17),
                                   Length = c(1, 7, 5, 4))
  expect_equal(testDigestDF$digestDf, expectedResult)
  expect_equal(testDigestDF$isCoDigest, TRUE)
})

test_that("simulateCoDigest informs user when no cleavage site is found", {
  expect_message(
    result <- simulateCoDigest(Biostrings::DNAString("AGGATAAACAA"),
                               Enzymes[c(1, 4), ]
    ),
    regexp = "No digestion site\\(s\\) found in sequence"
  )

  expect_null(result)
})


# [END]
