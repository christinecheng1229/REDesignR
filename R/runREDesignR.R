#' Launch Shiny App for REDesignR
#'
#' A function that launches the Shiny app for REDesignR.
#' The purpose of this app is to [...] The code of the app has been placed in \code{./inst/shiny-scripts}.
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
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
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
