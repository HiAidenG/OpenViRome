#' @title plotVirome
#' @description creates a pie chart of the number of palmids aligned to a given
#' family for each species in the genus-level virome.
#' @param virome A virome object
#' @return A plotly pie chart
#' @importFrom dplyr group_by %>% summarise arrange
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @export
plotVirome <- function(virome) {

    colors <- RColorBrewer::brewer.pal(9, "RdYlBu")

    virome <- virome %>%
    group_by(tax_family) %>%
    summarise(count=n()) %>%
    arrange(desc(count)) %>%
    plotly::plot_ly(labels=~tax_family, values=~count, type="pie",
                    textposition="inside", textinfo="label+percent",
                    hoverinfo="text", text=~paste(tax_family, count, sep=": "),
                    marker=list(colors=colors)) %>%
    plotly::layout(title="Virome Composition by Family")

}









