# viromeDashboard.R

#' @title viromeDashboard
#'
#' @description Load the OpenViRome shiny dashboard.
#'
#' @examples
#' viromeDashboard()
#'
#' @import shiny shinydashboard palmid plotly dplyr
#' @import RColorBrewer tidyr ggplot2 bs4Dash circlize 
#' @import grDevices fresh ComplexHeatmap
#'
#' @export
viromeDashboard <- function() {
  appDir <- system.file("shinyapp/app.R", package = "OpenViRome")
  if (appDir == "") {
    stop("Error: Could not find the app directory!")
  }
    shiny::runApp(appDir)
  }

# [END]
