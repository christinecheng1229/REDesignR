# Purpose: Visualization of restriction enzyme co-digestion simulation.
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.1
# Bugs and Issues: None known.

#' Visualize Restriction Enzyme co-digestion results
#'
#' TODO description
#'
#' @param codigestDf A tibble or data frame of co-digestion results in the format outputted by [simulateCoDigest()].
#' @param multiDigest A boolean indicating whether `codigestDf` resulted from a single-digest (`FALSE`) or a co-digest (`TRUE`), the latter being the default.
#' @param showLengths A boolean dictating whether numerical fragment lengths are shown as overlay on the plot; TRUE by default.
#'
#' @returns A plot of the results of a co-digestion experiment, either as an agarose gel or linear restriction map view and/or overlap of fragment lengths.
#'
#' @example path.R
#'
#' @references
#' Müller K, Wickham H (2025). _tibble: Simple Data Frames_. doi:10.32614/CRAN.package.tibble <https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0, <https://CRAN.R-project.org/package=tibble>.
#' Pedersen T (2025). _patchwork: The Composer of Plots._ R package version 1.3.2.9000, https://patchwork.data-imaginist.com.
#' Wickham, H. (2016). _ggplot2: Elegant Graphics for Data Analysis._ Springer-Verlag New York. ISBN: 978-3-319-24277-4.
#' Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. doi:10.32614/CRAN.package.dplyr <https://doi.org/10.32614/CRAN.package.dplyr>, R package version 1.1.4, <https://CRAN.R-project.org/package=dplyr>.
#'
#' @seealso [simulateCoDigest()]
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import tibble
#' @export
plotRestrictionMap <- function(codigestDf, multiDigest = TRUE, showLengths = TRUE) {
  # TODO: check this code in entirety

  # Standardize to tibble for downstream work
  codigestDf <- tibble::as_tibble(codigestDf)

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
  if (!all(requiredCols %in% names(codigestDf))) {
    stop("Please provide a valid digestion dataframe (refer to sample output of
         simulateCoDigest()).",
         call. = FALSE)
  }

  # ------------------------ Case: single-enzyme digest ------------------------
  if (multiDigest == FALSE) {
    codigestDf <- codigestDf %>%
      dplyr::mutate(
        Condition = Enzymes
      )

    plot <- singleRestrictionPlot(codigestDf, showLengths)

  } else {  # multiDigest == TRUE
    # ---------------------------- Case: co-digest ----------------------------
    # format: c("enzyme1 + enzyme2", "enzyme1", "enzyme2")
    enzymeNames <- unique(codigestDf$Enzymes)

    digestBoth <- codigestDf[codigestDf$Enzymes == enzymeNames[1], ]
    digestBoth$Condition <- rep(enzymeNames[1], nrow(digestBoth))
    mapBoth <- singleRestrictionPlot(digestBoth, showLengths)

    digestOne <- codigestDf[codigestDf$Enzymes == enzymeNames[2], ]
    digestOne$Condition <- rep(paste(enzymeNames[2], "only"), nrow(digestOne))
    mapOne <- singleRestrictionPlot(digestOne, showLengths)

    digestTwo <- codigestDf[codigestDf$Enzymes == enzymeNames[3], ]
    digestTwo$Condition <- rep(paste(enzymeNames[3], " only"), nrow(digestTwo))
    mapTwo <- singleRestrictionPlot(digestTwo, showLengths)

    # combine the 3 restriction maps into 1 plot using patchwork
    plot <- mapBoth / mapOne / mapTwo
  }

  return(plot)
}

# ----------------------- Helper: singleRestrictionPlot -----------------------
# Contains refactored code to create a single restriction map plot to
# reduce repetition in plotRestrictionMap's function body.
singleRestrictionPlot <- function(codigestDf, showLengths) {

  seq_length <- max(codigestDf$End)

  # Prepare dataframe for plotting (in descending fragment length)
  codigestDf <- codigestDf %>%
    dplyr::mutate(
      midpoint = (Start + End) / 2
    ) %>%
    dplyr::arrange(dplyr::desc(Length))

  # Linear map view (genomic order along the sequence)
  restrictionMap <- ggplot2::ggplot(
    codigestDf,
    ggplot2::aes(xmin = Start, xmax = End, ymin = 0, ymax = 1)
  ) +
    ggplot2::geom_rect(fill = "steelblue", color = "black", alpha = 0.8) +
    ggplot2::labs(
      title = paste("Restriction Fragment Map:", codigestDf$Condition[1]),
      x     = "Sequence position (bp)",
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  if (showLengths == TRUE) {
    # add fragment size as overlay labels
    restrictionMap <- restrictionMap +
      ggplot2::geom_text(
        ggplot2::aes(x = midpoint, y = 0.5, label = paste0(Length, " bp")),
        color = "white",
        size  = 3
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
#' @param codigestDf A tibble or data frame of co-digestion results in the format outputted by [simulateCoDigest()].
#' @param multiDigest A boolean indicating whether `codigestDf` resulted from a single-digest (`FALSE`) or a co-digest (`TRUE`), the latter being the default.
#' @param labelFragments  A boolean indicating whether plotted fragments should be labelled with the _FragmentID_ (`TRUE`) or not (`FALSE`).
#'
#' @examples
#' \dontrun{
#' TODO
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

  # Standardize to tibble
  codigestDf <- tibble::as_tibble(codigestDf)

  # ----------------- Input Validation -------------------

  # Ensure codigestDf is of correct format
  requiredCols <- c("Enzymes", "FragmentID", "Start", "End", "Length")
  if (!all(requiredCols %in% names(codigestDf))) {
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
        width  = 0.8,
        height = band_count / length(unique(Length)),
        alpha  = band_alpha
      ),
      fill  = "pink",
      color = "red"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_continuous(
      name   = "Fragment length (bp)",
      breaks = sort(unique(gel_df$Length), decreasing = TRUE),
      labels = sort(unique(gel_df$Length), decreasing = TRUE)
    ) +
    ggplot2::scale_x_continuous(
      name   = "Enzyme(s)",
      breaks = unique(gel_df$lane),
      labels = unique(gel_df$Enzymes)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
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
        .groups    = "drop"
      ) |>
      dplyr::mutate(
        x = lane,
        y = Length
      )

    plot <- plot +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(
          x     = x,
          y     = y,
          label = paste("frgmt", label_text)
        ),
        color = "black",
        size  = 4,
        hjust = 0.5
      )
  } else {  # labelFragments == FALSE
    # don't label fragments (do nothing)
  }

  return(plot)
}

# [END]
