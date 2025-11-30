# Purpose: Tests for simulateCoDigest function
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.0
# Bugs and Issues: None known.

test_that("Input validation: codigestDf format must be as documented", {
  testDigestDF <- simulateCoDigest(DNAString("ATCGATGGATCCATCGATATCGATGGATCC"), Enzymes[c(4, 27),])
  testthat::expect_error(REDesignR::plotFragments(t(testDigestDF)))
})

# [END]
