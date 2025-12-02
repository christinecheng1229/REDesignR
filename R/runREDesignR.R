# Purpose: Shiny App launch function.
# Author: Christine Cheng
# Date: November 30, 2025
# Version: 1.0
# Bugs and Issues: None known.

#' Launch Shiny App for REDesignR
#'
#' A function that launches the Shiny app for REDesignR.
#' The purpose of this app is to simulate restriction enzyme co-digestion and
#' generate restriction maps and gel simulations to aid in digestion experiment
#' design optimization. The code of the app has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but opens up a Shiny app.
#'
#' @examples
#' \dontrun{
#'
#' REDesignR::runREDesignR()
#' }
#'
#' @references
#' Grolemund, G. (2015). _Learn Shiny - Video Tutorials._
#'   https://shiny.rstudio.com/tutorial/
#'
#' _Shiny - File Download_. (2014, July 29). Shiny.
#'   https://shiny.posit.co/r/gallery/widgets/file-download/‌
#'
#' _Shiny - File Upload._ (2014, July 29). Shiny.
#'   https://shiny.posit.co/r/gallery/widgets/file-upload/‌
#'
#' @export
#' @import shiny
runREDesignR <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "REDesignR")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")

  return(actionShiny)
}
# [END]
