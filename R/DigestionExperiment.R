# Purpose: Restriction enzyme co-digestion simulation.
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.2
# Bugs and Issues: None known.

# ---------------------- Main function: simulateCoDigest ----------------------
#' Simulates co-restriction enzyme digestion
#'
#' A function that simulates a co-digestion experiment given a DNA sequence and the recognition sequences of 2 restriction enzymes.
#'
#' @param dnaSeq A DNAString representing the DNA sequence to be digested in the simulation.
#' @param enzymes  A table of 2 restriction enzymes containing the following at minimum:
#'   - Name - A string indicating the name of the restriction enzymes.
#'   - RecognitionSeq - A string representing the recognition sequences of the corresponding restriction enzyme, including the cleavage site.
#'
#' @returns Returns `invisible(NULL)` and a message if no digestion sites are found. Otherwise, returns a [tibble()] of fragments with columns:
#' \itemize{
#'   \item Enzymes - A string indicating the name of the restriction enzyme used for digestion.
#'   \item FragmentID - An integer representing the order of resulting digested DNA fragments.
#'   \item Start - An integer representing the starting position of the DNA fragment, in reference to the pre-digested DNA sequence.
#'   \item End - An integer representing the ending position of the DNA fragment, in reference to the pre-digested DNA sequence.
#'   \item Length - An integer representing the length of the resulting DNA fragment after digestion.
#' }
#'
#' @examples
#' # Using restriction enzymes available with package.
#' data(Enzymes, package = "REDesignR")  # load data
#'
#' # Example 1: ---------------------------------------------------------------
#' # Simulate co-digestion with RE AaaI and AagI.
#' simulateCoDigest(Biostrings::DNAString("CGGCCGATCGATCGGCCG"), REDesignR::Enzymes[c(1,4),])
#'
#' # Example 2: ---------------------------------------------------------------
#' # Simulate digestion when no cleavage sites are found.
#' simulateCoDigest(Biostrings::DNAString("AGGATAAACAA"), REDesignR::Enzymes[c(1,4),])
#'
#' @references
#' Müller K, Wickham H (2025). _tibble: Simple Data Frames_. doi:10.32614/CRAN.package.tibble <https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0, <https://CRAN.R-project.org/package=tibble>.
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2025). _Biostrings: Efficient manipulation of biological strings_. doi:10.18129/B9.bioc.Biostrings <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.78.0, <https://bioconductor.org/packages/Biostrings>.
#' Wickham H (2019). _assertthat: Easy Pre and Post Assertions_. doi:10.32614/CRAN.package.assertthat <https://doi.org/10.32614/CRAN.package.assertthat>, R package version 0.2.1, <https://CRAN.R-project.org/package=assertthat>.
#' Wright ES (2024). “Fast and Flexible Search for Homologous Biological Sequences with DECIPHER v3.” _The R Journal_, *16*(2), 191-200.
#'
#' @import Biostrings
#' @import DECIPHER
#' @import tibble
#' @export
simulateCoDigest <- function(dnaSeq, enzymes) {  # TODO: allow user to input their own dataset, otherwise use provided dataset as default

  # ------------------------ Input validation ---------------------------------
  # inform user with message describing cause of error
  # call. argument in stop() set to FALSE for a cleaner error view and easier
  # error message checks in testing

  # Argument: dnaSeq
  if (!inherits(dnaSeq, "DNAString")) {
    stop("Argument 'dnaSeq' must not be empty and must be a DNAString object.", call. = FALSE)
  }

  # Argument: enzymes
  if (!is.data.frame(enzymes) || nrow(enzymes) != 2) {
    stop("Argument 'enzymes' must not be empty and must be a data frame with exactly 2 rows.", call. = FALSE)
  }

  requiredCols <- c("Name", "RecognitionSeq")
  missingCols  <- setdiff(requiredCols, names(enzymes))

  if (length(missingCols) > 0) {
    stop(
      "Argument 'enzymes' is missing required column(s): ",
      paste(missingCols, collapse = ", "), call. = FALSE
    )
  }

  # # Recognition sequences of enzymes list
  # if (any(!nzchar(enzymes$RecognitionSeq))) {
  #   stop("Each enzyme must have a non-empty recognition sequence.", call. = FALSE)
  # }


  # process input to valid form expected by DECIPHER::DigestDNA
  dnaSet <- Biostrings::DNAStringSet(dnaSeq)

  # ------------------------ Run digestion -----------------------------

  firstDigest <- safeDigest(dnaSet, enzymes$RecognitionSeq[1])
  secondDigest <- safeDigest(dnaSet, enzymes$RecognitionSeq[2])

  # ------------------------ Case: no digestion at all ------------------------

  if (length(firstDigest) == 0 && length(secondDigest) == 0) {

    # Case: both enzymes don't digest given DNA sequence
    message("No digestion site(s) found in sequence")
    return(invisible(NULL)) # to stop function execution & prevent returning 'NULL' to user
  }

  # ------------------------ Successful digestion cases ----------------------------

  sequenceLen <- length(dnaSeq)

  if (length(firstDigest) == 0) {

    # Case: only second enzyme digests
    result <- buildDigestDf(
      fragments = secondDigest,
      enzymeName = enzymes$Name[2],
      sequenceLen = sequenceLen
    )

  } else if (length(secondDigest) == 0) {

    # Case: only first enzyme digests
    result <- buildDigestDf(
      fragments = firstDigest,
      enzymeName = enzymes$Name[1],
      sequenceLen = sequenceLen
    )

  } else {

    # Case: both enzymes digest
    mergedCuts <- sort(c(firstDigest, secondDigest)) # order digested positions

    result <- buildDigestDf(
      fragments = mergedCuts,
      enzymeName = paste(enzymes$Name, collapse = " + "),
      sequenceLen = sequenceLen
    )

  }

  return(result)

}

# ----------------------------- Helper: safeDigest -----------------------------
# Handles DECIPHER::DigestDNA errors when no cut locations exist.
# Returns integer vector of cut positions (possibly length 0).

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

# --------------------------- Helper: buildDigestDf ---------------------------
# Helper function that builds digest dataframe to reduce repeated code chunks.
# Returns tibble of digest results to be returned in simulateCoDigest

buildDigestDf <- function(fragments, enzymeName, sequenceLen) {

  # Add sequence boundaries
  cutPositions <- sort(unique(c(1L, fragments, sequenceLen)))
  fragLen      <- diff(cutPositions)

  # Create digest dataframe
  digestDf <- tibble::tibble(
    Enzymes    = enzymeName,
    FragmentID = seq_along(fragLen),
    Start      = head(cutPositions, -1),
    End        = tail(cutPositions, -1) - 1,
    Length     = fragLen
  )

  return(digestDf)

}

# [END]
