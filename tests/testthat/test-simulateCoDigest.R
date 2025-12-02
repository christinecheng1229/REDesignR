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
  expectedResult <- tibble::tibble(Enzymes = c(rep(c("AaaI + AagI"), times = 4),
                                               rep(c("AaaI"), times = 3),
                                               rep(c("AagI"), times = 2)),
                                   FragmentID = c(1, 2, 3, 4, 1, 2, 3, 1, 2),
                                   Start = c(1, 2, 9, 14, 1, 2 , 14, 1, 9),
                                   End = c(1, 8, 13, 17, 1, 13, 17, 8, 17),
                                   Length = c(1, 7, 5, 4, 1, 12, 4, 8, 9))
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

test_that("simulateCoDigest rejects invalid dnaSeq input", {
  data(Enzymes, package = "REDesignR")

  expect_error(
    simulateCoDigest("NOT_A_DNASTRING", Enzymes[c(1, 4), ]),
    "must be a DNAString object"
  )
})

test_that("simulateCoDigest rejects invalid enzyme table size", {
  data(Enzymes, package = "REDesignR")

  one_enzyme <- Enzymes[1, ]
  three_enzymes <- Enzymes[c(1, 2, 3), ]

  expect_error(simulateCoDigest(Biostrings::DNAString("ATCG"), one_enzyme))
  expect_error(simulateCoDigest(Biostrings::DNAString("ATCG"), three_enzymes))
})

test_that("simulateCoDigest rejects enzyme table missing required columns", {
  bad_tbl <- tibble::tibble(NotName = c("a", "b"), Wrong = c("x", "y"))

  expect_error(
    simulateCoDigest(Biostrings::DNAString("ATCG"), bad_tbl),
    "missing required column"
  )
})

test_that("simulateCoDigest produces expected tibble output for single-digest case", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("ATGCGGCCGTTT")

  result <- simulateCoDigest(dna, Enzymes[c(1, 6), ])

  expect_type(result, "list")

  expect_false(result$isCoDigest)

  expect_s3_class(result$digestDf, "tbl_df")
  expect_true(all(c("Enzymes", "FragmentID", "Start", "End", "Length") %in%
                    names(result$digestDf)))
})

test_that("simulateCoDigest produces expected tibble output for co-digest case", {
  data(Enzymes, package = "REDesignR")
  dna <- Biostrings::DNAString("CGGCCGATCGATCGGCCG")

  result <- simulateCoDigest(dna, Enzymes[c(1, 4), ])

  expect_type(result, "list")

  expect_true(result$isCoDigest)

  expect_s3_class(result$digestDf, "tbl_df")
  expect_true(all(c("Enzymes", "FragmentID", "Start", "End", "Length") %in%
                    names(result$digestDf)))
})

# [END]
