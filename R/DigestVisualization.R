# Purpose: Code and documentation for visualization function(s).
# Author: Christine Cheng
# Date: November 4, 2025
#' Visualize distribution of digest fragments' characteristic(s)
#'
#' TODO description
#'
#' TODO examples
#'
#'
#' @import ggplot2
#' @import dplyr
#' @export
plotFragments <- function(codigestDf, seq_length = NULL, show_lengths = TRUE, gel_style = TRUE) {
  # TODO: check this code in entirety

  # Basic input checks
  if (missing(codigestDf) || !"FragmentID" %in% names(codigestDf)) {
    stop("Please provide a valid codigestion dataframe (output of simulateCoDigest()).")
  }

  # Infer sequence length if not provided
  if (is.null(seq_length)) {
    seq_length <- max(codigestDf$End)
  }

  # Prepare dataframe for plotting (shortest fragment at the top)
  codigestDf <- codigestDf %>%
    dplyr::mutate(
      y_pos = max(FragmentID) - FragmentID + 1,
      midpoint = (Start + End) / 2
    ) %>%
    dplyr::arrange(Length)

  if (gel_style) {
    # Simulated gel view (bands proportional to fragment length)
    plot <- ggplot(codigestDf, aes(y = y_pos, x = Length)) +
      geom_tile(aes(width = Length / seq_length * 0.9, height = 0.8),
                fill = "steelblue", alpha = 0.8, color = "black") +
      scale_y_continuous(breaks = codigestDf$y_pos, labels = codigestDf$FragmentID) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      labs(
        title = paste("Simulated Co-Digestion Pattern:", codigestDf$Enzyme[1]),
        x = "Fragment length (bp)",
        y = "Fragment ID"
      ) +
      theme_minimal(base_size = 13)

    if (show_lengths) {
      plot <- plot + geom_text(aes(label = Length), color = "white", size = 4, vjust = 0.4)
    }
  } else {
    # Linear map view (genomic order along the sequence)
    plot <- ggplot(codigestDf, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1)) +
      geom_rect(fill = "steelblue", color = "black", alpha = 0.8) +
      labs(
        title = paste("Restriction Fragment Map:", codigestDf$Enzyme[1]),
        x = "Sequence position (bp)",
        y = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

    if (show_lengths) {
      plot <- plot + geom_text(
        aes(x = (Start + End) / 2, y = 0.5, label = paste0(Length, " bp")),
        color = "white", size = 3
      )
    }
  }

  return(plot)
}


# [END]
