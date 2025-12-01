# Purpose: Visualization of restriction enzyme co-digestion simulation.
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.1
# Bugs and Issues: None known.
#' Visualize Restriction Enzyme co-digestion results
#'
#' TODO description
#'
#' @param codigestDf A dataframe of codigestion results in the format outputted by [simulateCoDigest()].
#' @param show_lengths A boolean dictating whether numerical fragment lengths are shown as overlay on the plot; TRUE by default.
#'
#' @returns A Plot of the results of a co-digestion experiment, either as an agarose gel or linear restriction map view and/or overlap of fragment lengths.
#'
#' TODO examples
#'
#' @references
#' Pedersen T (2025). _patchwork: The Composer of Plots._ R package version 1.3.2.9000, https://patchwork.data-imaginist.com.
#'
#' @seealso [simulateCoDigest()]
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
plotRestrictionMap <- function(codigestDf, multiDigest = TRUE, showLengths = TRUE) {
  # TODO: check this code in entirety

  # --------- Input Validation ----------------------------------

  # Ensure codigestDf is not NULL and is non-empty
  if (is.null(codigestDf) | nrow(codigestDf) == 0) {
    stop("No digestion detected. Please provide a valid digestion dataframe
         (refer to sample output of simulateCoDigest()).",
         call. = FALSE
    )
  }

  # Ensure codigestDf is of correct format
  requiredCols <- c("Enzymes", "FragmentID", "Start", "End", "Length")
  if (!is.data.frame(codigestDf) | !all(requiredCols %in% names(codigestDf))) {
    stop("Please provide a valid digestion dataframe (refer to sample output of
         simulateCoDigest()).",
         call. = FALSE)
  }
  # if (!"Enzymes" %in% names(codigestDf) ||
  #     !"FragmentID" %in% names(codigestDf) ||
  #     !"Start" %in% names(codigestDf) ||
  #     !"End" %in% names(codigestDf) ||
  #     !"Length" %in% names(codigestDf)) {
  #   stop("Please provide a valid digestion dataframe (refer to sample output of
  #        simulateCoDigest()).",
  #        call. = FALSE)
  # }

  # ------------------------ Case: single-enzyme digest ------------------------
  if (multiDigest == FALSE) {
    digestDf <- codigestDf %>%
      dplyr::mutate(
        Condition = Enzymes,
        midpoint = (Start + End) / 2
      )

    plot <- singleRestrictionPlot(codigestDf, showLengths)
  } else {  # multiDigest == TRUE
    # ---------------------------- Case: co-digest ----------------------------
    # format: c("enzyme1 + enzyme2", "enzyme1", "enzyme2")
    enzymeNames <- unique(codigestDf$Enzymes)
    codigestDf <- codigestDf %>%
      dplyr::mutate(midpoint = (Start + End) / 2)

    digestBoth <- codigestDf[codigestDf$Enzymes == enzymeNames[1], ]
    digestBoth$Condition <- enzymeNames[1]
    mapBoth <- singleRestrictionPlot(digestBoth, showLengths)

    digestOne <- codigestDf[codigestDf$Enzymes == enzymeNames[2], ]
    digestOne$Condition <- paste(enzymeNames[2], "only")
    mapOne <- singleRestrictionPlot(digestOne, showLengths)

    digestTwo <- codigestDf[codigestDf$Enzymes == enzymeNames[3], ]
    digestTwo$Condition <- paste(enzymeNames[3], " only")
    mapTwo <- singleRestrictionPlot(digestTwo, showLengths)

    # combine the 3 restriction maps into 1 plot using patchwork
    plot <- mapBoth / mapOne / mapTwo

  }

  return(plot)
}

# ----------------------- Helper: singleRestrictionPlot -----------------------
# Contains refactored code to create a single restriction map plot for to
# reduce repetition in plotRestrictionMap's function body.
singleRestrictionPlot <- function(codigestDf, showLengths) {
  seq_length <- max(codigestDf$End)

  # Prepare dataframe for plotting (in descending fragment length)
  codigestDf <- codigestDf %>%
    dplyr::arrange(dplyr::desc(Length))

  # Linear map view (genomic order along the sequence)
  restrictionMap <- ggplot2::ggplot(codigestDf, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1)) +
    ggplot2::geom_rect(fill = "steelblue", color = "black", alpha = 0.8) +
    ggplot2::labs(
      title = paste("Restriction Fragment Map:", codigestDf$Condition[1]), #***
      x = "Sequence position (bp)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())

  if (showLengths == TRUE) {
    # add fragment size as overlay labels
    restrictionMap <- restrictionMap + ggplot2::geom_text(
      aes(x = midpoint, y = 0.5, label = paste0(Length, " bp")),
      color = "white", size = 3
    )
  } else {
    # don't add overlay labels (don't do anything)
  }

  return(restrictionMap)
}

#' Agarose Gel Simulation of RE Digests
#'
#' TODO description
#'
#' @param codigestDf A dataframe of codigestion results in the format outputted by [simulateCoDigest()].
#' @param multiDigest A boolean indicating whether `codigestDf` resulted from a single-digest (`FALSE`) or a co-digest (`TRUE`), the latter being the default.
#' @param labelFragments  A boolean indicating whether plotted fragments should be labelled with the _FragmentID_ (`TRUE`) or not (`FALSE`).
#'
#' @examples
#' \dontrun{
#' TODO
#' data(Enzymes, package = "REDesignR")  # load data
#' # Simulate co-digestion with RE AaaI and AagI.
#' TODO
#' coDigestDf <- simulateCoDigest(Biostrings::DNAString("CGGCCGATCGATCGGCCG"), REDesignR::Enzymes[c(1,4),])
#'
#' # Example 1: ---------------------------------------------------------------
#' simulateGel(coDigestDf) # by default: multiDigest = TRUE, labelFragments = FALSE
#'
#' # Example 2: ---------------------------------------------------------------
#' simulateGel(coDigestDf, labelFragments = TRUE)  # by default: multiDigest = TRUE
#'
#' # Example 3: ---------------------------------------------------------------
#' # Simulate single-digestion with TODO
#' TODO
#' digestDF <- TODO
#' simulateGel(digestDF, multiDigest = FALSE)  # by default: labelFragments = FALSE
#' # > REDesignR::simulateGel(testDigestDF)
#' # > REDesignR::simulateGel(testDigestDF, multiDigest = TRUE)
#' # > REDesignR::simulateGel(testDigestDF, labelFragments = TRUE)
#' # > REDesignR::simulateGel(testDigestDF2)
#' # > REDesignR::simulateGel(testDigestDF2, multiDigest = FALSE)
#' # > REDesignR::simulateGel(testDigestDF2)
#' # > REDesignR::simulateGel(testDigestDF2, multiDigest = FALSE)
#' # > REDesignR::simulateGel(testDigestDF2, multiDigest = FALSE, TRUE)
#' }
#'
#' @returns A plot representing an agarose gel of the fragments resulting from [simulateCoDigest()].
#'
#' @references
#' Wickham H, Pedersen T, Seidel D (2025). _scales: Scale Functions for Visualization._ R package version 1.4.0, https://scales.r-lib.org.
#'
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import scales
simulateGel <- function(codigestDf,
                        multiDigest = TRUE,
                        labelFragments = FALSE) {
  # ----------------- Input Validation -------------------

  # Ensure codigestDf is of correct format
  requiredCols <- c("Enzymes", "FragmentID", "Start", "End", "Length")
  if (!is.data.frame(codigestDf) | !all(requiredCols %in% names(codigestDf))) {
    stop("Please provide a valid digestion dataframe (refer to sample output of
         simulateCoDigest()).",
         call. = FALSE)
  }

  # Ensure codigestDf is not NULL and is non-empty
  if (is.null(codigestDf) | nrow(codigestDf) == 0) {
    stop("No digestion detected. Please provide a valid digestion dataframe
         (refer to sample output of simulateCoDigest()).",
         call. = FALSE
    )
  }

  # Ensure digest type from dataframe matches multiDigest boolean
  if (length(unique(codigestDf$Enzymes)) == 1 & multiDigest == TRUE) {
    warning("Single-enzyme digestion detected. Please double-check inputted
            codigestDf if co-digestion was intended.",
            call. = FALSE)
  }
  if (length(unique(codigestDf$Enzymes)) > 1 & multiDigest == FALSE) {
    warning("Multi-enzyme digestion detected. Please double-check inputted
            codigestDf if single-digestion was intended.",
            call. = FALSE)
  }

  # ------------------------ Pre-processing ------------------------
  # Fragment grouping: all fragments of a single digest go in a single column
  gel_df <- codigestDf |>
    dplyr::mutate(
      lane = as.numeric(factor(Enzymes))
    ) |>
    dplyr::group_by(Enzymes, Length) |>
    dplyr::mutate(
      # band intensity: count of identical fragment lengths to simulate
      # realistic gel intensity variance
      band_count = dplyr::n(),
      band_alpha = scales::rescale(band_count, to = c(0.2, 1))
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(Length))

  # ------------------------- Plot ------------------

  plot <- ggplot2::ggplot(
    gel_df,
    ggplot2::aes(x = lane, y = Length)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(
        width = 0.8,
        height = band_count / length(unique(Length)),
        alpha = band_alpha
      ),
      fill = "pink",
      color = "red"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_continuous(
      name = "Fragment length (bp)",
      breaks = sort(unique(gel_df$Length), decreasing = TRUE),
      labels = sort(unique(gel_df$Length), decreasing = TRUE)
    ) +
    ggplot2::scale_x_continuous(
      name = "Enzyme(s)",
      breaks = unique(gel_df$lane),
      labels = unique(gel_df$Enzymes)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks.x = ggplot2::element_blank()
    )

  if (multiDigest == TRUE) {
    plot <- plot +
      ggplot2::labs(
        title = "Simulated Agarose Gel (Multi-Lane)"
      )
  } else {  # multiDigest == FALSE
    plot <- plot +
      ggplot2::labs(
        title = paste("Simulated Agarose Gel:", gel_df$Enzymes[1])
      )
  }

  # ------------------- Fragment labels -----------------------

  if (labelFragments == TRUE) {
    label_df <- gel_df |>
      dplyr::group_by(Enzymes, Length, lane) |>
      dplyr::summarise(
        label_text = paste(FragmentID, collapse = ", "),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        x = lane,
        y = Length
      )

    plot <- plot +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(
          x = x,
          y = y,
          label = paste("frgmt", label_text)
        ),
        color = "black",
        size = 4,
        hjust = 0.5
      )
  } else {  # labelFragments == FALSE
    # don't label fragments (do nothing)
  }

  return(plot)
}


# ---------- rough work -------------------------

# gel_df <- codigestDf |>
#   dplyr::group_by(Length) |>
#   dplyr::mutate(
#     band_count = dplyr::n(),
#     # Scaling for visual thickness
#     band_alpha = scales::rescale(band_count, to = c(0.2, 1))
#   ) |>
#   dplyr::ungroup() |>
#   dplyr::arrange(dplyr::desc(Length))
#
# # --------------------------------------------------------------------
# # Plot: tile placed at x = 1; identical lengths overlay -> thicker band
# # --------------------------------------------------------------------
# plot <- ggplot2::ggplot(gel_df, ggplot2::aes(x = 1, y = Length)) +
#   ggplot2::geom_tile(
#     ggplot2::aes(
#       width = 0.8,
#       height = band_count / length(unique(Length)),
#       alpha = band_alpha
#     ),
#     fill = "pink",
#     color = "red"
#   ) +
#   ggplot2::scale_alpha_identity() +
#   ggplot2::scale_y_continuous(
#     name = "Fragment length (bp)",
#     breaks = sort(unique(gel_df$Length), decreasing = TRUE),
#     labels = sort(unique(gel_df$Length), decreasing = TRUE)
#   ) +
#   ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.3, 0.3))) +
#   ggplot2::labs(
#     title = paste("Simulated Agarose Gel:", gel_df$Enzymes[1]),
#     x = NULL
#   ) +
#   ggplot2::theme_minimal(base_size = 14) +
#   ggplot2::theme(
#     axis.text.x = ggplot2::element_blank(),
#     axis.ticks.x = ggplot2::element_blank()
#   )
#
# # ------------------------------------------------
# if (labelFragments == TRUE) {
#   # Create one text label per unique fragment length
#   label_df <- gel_df |>
#     dplyr::group_by(Length) |>
#     dplyr::summarise(
#       label_text = paste(FragmentID, collapse = ", "),
#       .groups = "drop"
#     ) |>
#     dplyr::mutate(
#       x = 1.3,      # horizontal placement beside the band
#       y = Length    # same y as the band
#     )
#
#   # Add to plot
#   plot <- plot +
#     ggplot2::geom_text(
#       data = label_df,
#       ggplot2::aes(
#         x = x,
#         y = y,
#         label = label_text
#       ),
#       color = "black",
#       size = 4,
#       hjust = 1     # left aligned
#     )
# }

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
