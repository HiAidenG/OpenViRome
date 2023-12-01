#' @title plotVirome
#' @description creates a pie chart of the number of palmids aligned to a given
#' family for each species in the genus-level virome.
#' @param virome A virome object
#' @param sourceSpecies A character string specifying the source species to
#' filter the virome by. Default is NULL and will not filter the data.
#' @return A plotly pie chart
#' @importFrom dplyr group_by %>% summarise arrange
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @export
plotVirome <- function(virome = NULL, sourceSpecies = NULL) {

    virome <- virome[[1]]

    # Check if source species is in virome
    if (!is.null(sourceSpecies)) {
      if (!sourceSpecies %in% virome$scientific_name) {
        stop("sourceSpecies not found in virome")
      }
      # Filter virome by source species
      virome <- virome %>%
        dplyr::filter(scientific_name == sourceSpecies)
    }


    # Get distinct sOTUs
    sotus <- virome %>%
      dplyr::select(sotu) %>%
      dplyr::distinct()

    sotus <- sotus %>%
      dplyr::left_join(virome, by = "sotu") %>%
      dplyr::group_by(tax_phylum) %>%
      dplyr::summarise(count = n())

    # Rename null to unknown
    sotus$tax_phylum[is.na(sotus$tax_phylum)] <- "Unknown"

    colors <- brewer.pal(length(unique(virome$tax_phylum)), "Set3")

    p <- plotly::plot_ly(sotus, labels=~tax_phylum, values=~count, type="pie",
                    textposition="inside", textinfo="label+percent",
                    hoverinfo="text", text=~paste(tax_phylum, count, sep=": "),
                    marker=list(colors=colors)) %>%
    plotly::layout(title="Unique sOTU Count by Predicted Phylum")


    return(p)

}









