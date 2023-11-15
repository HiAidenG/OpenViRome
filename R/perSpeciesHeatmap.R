#' @title perSpeciesHeatmap
#' @description creates a heatmap of the number of palmids aligned to a given
#' family for each species in the genus-level virome.
#' @param virome A virome object
#' @return A plotly heatmap
#' @importFrom dplyr group_by %>% summarise spread
#' @importFrom plotly plot_ly layout colorbar
#' @importFrom RColorBrewer brewer.pal colorRampPalette
#' @export
perSpeciesHeatmap <- function(virome) {
  countsDF <- virome %>%
    group_by(scientific_name, tax_family) %>%
    summarise(count = n()) %>%
    spread(tax_family, count, fill = 0)

  countsM <- t(as.matrix(countsDF[,2:ncol(countsDF)]))
  x <- countsDF$scientific_name
  y <- colnames(countsDF)[2:ncol(countsDF)]

  colors <- RColorBrewer::brewer.pal(9, "RdYlBu")
  colors <- colorRampPalette(colors)(100)
  colors <- rev(colors)
  colors <- c("white", colors)

  fig <- plotly::plot_ly(y=y, x=x, z=countsM, type="heatmap", colors=colors) %>%
   plotly::layout(xaxis=list(title="Source Species"), yaxis=list(autorange=TRUE)) %>%
   plotly::colorbar(title="# of aligned palmids")


  return(fig)
}
