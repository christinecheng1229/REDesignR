# Purpose: Code and documentation for RE digestion simulation function(s).
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.1
# Bugs and Issues: None known.
#' Simulates multi-restriction enzyme digestion
#'
#' A function that simulates a co-digestion experiment given a DNA sequence and the recognition sequences of 2 restriction enzymes.
#'
#' @param DNA A DNAString representing the DNA sequence to be digested in the simulation.
#' @param enzymes  A table of 2 restriction enzymes containing the following at minimum:
#' \itemize {
#'  \item Name A string indicating the name of the restriction enzymes.
#'  \item RecognitionSeq A string representing the recognition sequences of the corresponding restriction enzyme, including the cleavage site.
#' }
#'
#' @returns Returns a dataframe of digestion results if successful:
#' \itemize{
#'   \item Enzymes - A string indicating the name of the restriction enzyme used for digestion.
#'   \item FragmentID - An integer representing the order of resulting digested DNA fragments.
#'   \item Start - An integer representing the starting position of the DNA fragment, in reference to the pre-digested DNA sequence.
#'   \item End - An integer representing the ending position of the DNA fragment, in reference to the pre-digested DNA sequence.
#'   \item Length - An integer representing the length of the resulting DNA fragment after digestion.
#' }
#' @returns otherwise informs user of unsuccessful digestion.
#'
#' @examples
#' # Using restriction enzymes available with package.
#'
#' # Example 1:
#' # Simulate co-digestion with RE AaaI and AagI.
#' simulateCoDigest(Biostrings::DNAString("CGGCCGATCGATCGGCCG"), Enzymes[c(1,4),])
#'
#' # Example 2:
#' # Simulate digestion when no cleavage sites are found.
#' simulateCoDigest(Biostrings::DNAString("AGGATAAACAA"), Enzymes[c(1,4),])
#'
#' @references
#' Müller K, Wickham H (2025). _tibble: Simple Data Frames_. doi:10.32614/CRAN.package.tibble <https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0, <https://CRAN.R-project.org/package=tibble>.
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2025). _Biostrings: Efficient manipulation of biological strings_. doi:10.18129/B9.bioc.Biostrings <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.78.0, <https://bioconductor.org/packages/Biostrings>.
#' Pagès H, Lawrence M, Aboyoun P (2025). _S4Vectors: Foundation of vector-like and list-like containers in Bioconductor_. doi:10.18129/B9.bioc.S4Vectors <https://doi.org/10.18129/B9.bioc.S4Vectors>, R package version 0.48.0, <https://bioconductor.org/packages/S4Vectors>.
#' Wickham H (2019). _assertthat: Easy Pre and Post Assertions_. doi:10.32614/CRAN.package.assertthat <https://doi.org/10.32614/CRAN.package.assertthat>, R package version 0.2.1, <https://CRAN.R-project.org/package=assertthat>.
#' Wright ES (2024). “Fast and Flexible Search for Homologous Biological Sequences with DECIPHER v3.” _The R Journal_, *16*(2), 191-200.
#'
#' @import assertthat
#' @import Biostrings
#' @import DECIPHER
#' @import S4Vectors
#' @import tibble
#' @export
simulateCoDigest <- function(DNA, enzymes) {  # TODO: allow user to input their own dataset, otherwise use provided dataset as default
  # check for valid non-empty inputs
  assertthat::assert_that(not_empty(DNA), not_empty(enzymes))

  seqSet <- Biostrings::DNAStringSet(DNA)  # process input to valid form expected by DECIPHER::DigestDNA

  sequenceLen <- length(DNA)

  # Run the digestion experiments using the two enzymes
  firstDigest <- safeDigest(seqSet, enzymes[1, ]$RecognitionSeq)
  secondDigest <- safeDigest(seqSet, enzymes[2, ]$RecognitionSeq)

  if (length(firstDigest) == 0 && length(secondDigest) == 0) {

    # Case: both enzymes don't digest given DNA sequence
    message("No digestion site(s) found in sequence")
    return(invisible(NULL)) # to stop function execution & prevent returning 'NULL' to user

  } else if (S4Vectors::isEmpty(firstDigest) == TRUE) {  # Case: only second enzyme digests

    # Add sequence boundaries
    cutPositions <- unique(c(1L, secondDigest, sequenceLen))
    cutPositions <- sort(cutPositions)
    fragLen <- diff(cutPositions) # Compute fragment lengths

    # Update results table
    enzymeNames <- enzymes$Name[2]  # TODO: 'Enzymes' should be 'enzymes' (input argument), right?
    # result$Enzymes = Enzymes$Name[2]
    # result$FragmentID <- seq_along(fragLen)
    # result$Start <- head(cutPositions, -1)
    # result$End <- tail(cutPositions, -1)
    # result$Length <- fragLen

  } else if (S4Vectors::isEmpty(secondDigest) == TRUE) { # Case: only first enzyme digests

    # Add sequence boundaries
    cutPositions <- unique(c(1L, firstDigest, sequenceLen))
    cutPositions <- sort(cutPositions)
    fragLen <- diff(cutPositions) # Compute fragment lengths

    # Update results table
    enzymeNames <- enzymes$Name[1] # TODO: 'Enzymes' should be 'enzymes' (input argument), right?
    # result$Enzymes = Enzymes$Name[1]
    # result$FragmentID <- seq_along(fragLen)
    # result$Start <- head(cutPositions, -1)
    # result$End <- tail(cutPositions, -1)
    # result$Length <- fragLen

  } else {  # Case: both enzymes digest

    cutPositions <- sort(c(firstDigest, secondDigest)) # order digested positions

    # Add sequence boundaries
    cutPositions <- unique(c(1L, cutPositions, sequenceLen))
    cutPositions <- sort(cutPositions)
    fragLen <- diff(cutPositions) # Compute fragment lengths

    # Update results table
    enzymeNames <- paste(enzymes$Name, collapse = " + ")
    # result$FragmentID <- seq_along(fragLen)
    # result$Start <- head(cutPositions, -1)
    # result$End <- tail(cutPositions, -1)
    # result$Length <- fragLen

  }
  result <- tibble::tibble(Enzymes = enzymeNames,
                           FragmentID = seq_along(fragLen),
                           Start = head(cutPositions, -1),
                           End = tail(cutPositions, -1) - 1,
                           Length = fragLen)
  return(result)
}

#' Helper function of simulateCoDigest to handle errors thrown by DECIPHER::DigestDNA
#' @keywords internal
safeDigest <- function(dnaSeq, recognitionSite) {

  tryCatch(
    {
      DECIPHER::DigestDNA(
        myDNAStringSet = dnaSeq,
        sites = recognitionSite,
        type = "positions",
        strand = "top"
      )[[1]]$top
    },

    error = function(e) {
      if (grepl("No cut location\\(s\\) found in site:", conditionMessage(e))) {

        # Specific handling for NO CUT SITE error: return no cut site result so
        # that safeDigest always returns an integer vector for cleaner
        # downstream processing
        return(integer(0))

      } else {
        # Re-throw unexpected errors
        stop(e)
      }
    }
  )
}

# [END]
