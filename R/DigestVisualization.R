# Purpose: Code and documentation for visualization function(s).
# Author: Christine Cheng
# Date: November 18, 2025
# Version: 1.0
# Bugs and Issues: None known.
#' Visualize Restriction Enzyme co-digestion results
#'
#' TODO description
#'
#' @param codigestDf A dataframe of codigestion results in the format outputted by [simulateCoDigest()].
#' @param seq_length  An integer representing the length (# of nucleotides) in DNA sequence; inferred by the function if not provided.
#' @param show_lengths A boolean dictating whether numerical fragment lengths are shown as overlay on the plot; TRUE by default.
#' @param gel_style A boolean dictating how digested fragments are visualized; FALSE by default:
#' \itemize {
#'    \item TRUE = simulated agarose gel (bands by decreasing fragment length) view
#'    \item FALSE = linear restriction map view
#' }
#'
#' @returns A Plot of the results of a co-digestion experiment, either as an agarose gel or linear restriction map view and/or overlap of fragment lengths.
#'
#' TODO examples
#'
#' @seealso [simulateCoDigest()]
#'
#' @import ggplot2
#' @import dplyr
#' @export
plotFragments <- function(codigestDf, seq_length = NULL, show_lengths = TRUE, gel_style = FALSE) {
  # TODO: check this code in entirety

  # Basic input checks
  if (!"Enzymes" %in% names(codigestDf) ||
      !"FragmentID" %in% names(codigestDf) ||
      !"Start" %in% names(codigestDf) ||
      !"End" %in% names(codigestDf) ||
      !"Length" %in% names(codigestDf)) {
    stop("Please provide a valid codigestion dataframe (refer to sample output of simulateCoDigest()).")
  }

  # Infer sequence length if not provided
  if (is.null(seq_length)) {
    seq_length <- max(codigestDf$End)
  }

  # Prepare dataframe for plotting (in descending fragment length)
  codigestDf <- codigestDf %>%
    dplyr::mutate(
      y_pos = max(FragmentID) - FragmentID + 1, # TODO not needed if breaks on lines #* can use FragmentID instead of y_pos
      midpoint = (Start + End) / 2
    ) %>%
    dplyr::arrange(dplyr::desc(Length))

  if (gel_style) {
    # Simulated agarose gel view (bands proportional to fragment length)
    plot <- ggplot2::ggplot(codigestDf, aes(y = FragmentID, x = Length)) + #*
      ggplot2::geom_tile(aes(width = Length / seq_length * 0.9, height = 0.8),
                         fill = "steelblue",
                         alpha = 0.8,
                         color = "black") +
      ggplot2::scale_y_continuous(breaks = codigestDf$FragmentID, #*
                                  labels = codigestDf$FragmentID) +  #*
      # ggplot2::scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      ggplot2::labs(
        title = paste("Simulated Co-Digestion Pattern:", codigestDf$Enzymes[1]), # TODO
        x = "Fragment length (bp)",
        y = "Fragment ID"
      ) +
      ggplot2::theme_minimal(base_size = 13)

    if (show_lengths) {
      # add fragment size overlay labels
      plot <- plot + ggplot2::geom_text(aes(label = Length),
                                        color = "white",
                                        size = 4,
                                        vjust = 0.4)
    }
  } else {
    # Linear map view (genomic order along the sequence)
    plot <- ggplot2::ggplot(codigestDf, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1)) +
      ggplot2::geom_rect(fill = "steelblue", color = "black", alpha = 0.8) +
      ggplot2::labs(
        title = paste("Restriction Fragment Map:", codigestDf$Enzymes[1]), # TODO
        x = "Sequence position (bp)",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

    if (show_lengths) {
      # add fragment size overlay labels
      plot <- plot + ggplot2::geom_text(
        aes(x = midpoint, y = 0.5, label = paste0(Length, " bp")),
        color = "white", size = 3
      )
    }
  }

  return(plot)
}

#' Visualize distribution of digest fragments' characteristic(s)
#'
#' TODO description
#'
#' @param fragmentStats
#'
#' @returns
#'
#' TODO examples
#'
#' @seealso [simulateCoDigest()]
#'
#' @export
plotFragmentStats <- function(fragmentStats) {
  # idea: output a collection of mini plots that visualize all sorts of metadata
  # about the fragments (e.g., length distribution, distribution between 2 enzymes' cuts, etc.)
}

#' TODO Helper function: title
#'
#' @param codigestDF A dataframe of co-digestion results outputted by [simulateCoDigest()]
#'
#' @returns A dataframe of co-digestion results' metadata/stats to be used by [plotFragmentStats()] and the Shiny app.
getFragmentStats <- function(codigestDF) {
  fragmentStats <- data.frame()
  # include metadata of digestion to be plotted and simplier ones (e.g., average frag len) to be included in shiny app table
  return(fragmentStats)
}


# [END]
