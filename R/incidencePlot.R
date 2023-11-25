#' @title plotVirusPositive
#' @description Plots the proportion of virus positive runs for each taxon
#' @param abundanceDF
#' @return A plotly bar chart
#' @importFrom dplyr group_by %>% summarise collect left_join
#' @importFrom tidyr pivot_longer
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @export
plotVirusPositive <- function(abundanceDF) {
  # Add a column for virus negative counts
  abundanceDF <- abundanceDF %>%
    dplyr::mutate(virus_negative = total - virus_positive)

  # Sort the dataframe by total count in descending order
  abundanceDF <- abundanceDF %>%
    dplyr::arrange(desc(total))

  # Create a long format dataframe for plotting
  longData <- abundanceDF %>%
    dplyr::select(scientific_name, virus_positive, virus_negative, total) %>%
    tidyr::pivot_longer(
      cols = c(virus_positive, virus_negative),
      names_to = "virus_status",
      values_to = "count"
    )

  # Adjust the levels of the virus_status factor to control the stack order
  longData$virus_status <- factor(longData$virus_status, levels = c("virus_positive", "virus_negative"))

  # Create the horizontal bar plot
  p <- plot_ly(data = longData, y = ~scientific_name, x = ~count, type = 'bar', color = ~virus_status,
               colors = c("virus_positive" = "cornflowerblue", "virus_negative" = "grey"),
               orientation = 'h') %>%
    layout(barmode = "stack", xaxis = list(title = "Count"), yaxis = list(title = "Run Source Species", categoryorder = "total ascending")) %>%
    plotly::config(displayModeBar = FALSE)

  return(p)
}





