# viromeDashboard.R

#' @title viromeDashboard
#'
#' @description Load the OpenViRome shiny dashboard.
#'
#' @examples
#' viromeDashboard()
#'
#' @import shiny
#' @import shinydashboard
#' @import palmid
#' @import plotly
#' @import dplyr
#' @import RColorBrewer
#' @import tidyr
#' @import ggplot2
#' @import bs4Dash
#' @import circlize
#' @import grDevices
#' @import fresh
#' @import ComplexHeatmap
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
