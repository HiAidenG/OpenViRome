# plotVirome.R

#' @title plotViromePie
#'
#' @description Returns a plotly pie chart of predicted phyla for unique sOTUs
#' in the virome object.
#'
#' @param virome A virome object
#' @param phylumFilter Optional argument for filtering the virome by phylum.
#' Default is NULL and will not filter the data. Must be a character vector.
#' @param colors Optional argument for specifying colors for the pie chart.
#' Default is NULL and will use the Okabe-Ito palette if there are 9 or fewer
#' phyla. Otherwise, uses RColorBrewer Set3. Must be a character vector of hex
#' colors. length(colors) == length(unique(virome$tax_phylum))
#'
#' @return A plotly pie chart
#'
#' @importFrom dplyr group_by %>% summarise arrange
#' @importFrom plotly plot_ly layout
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Meloidogyne", con = con)
#' plotViromePie(virome = virome,
#'            phylumFilter = c("Pisuviricota", "Kitrinoviricota"),
#'            colors = c("#E69F00", "#56B4E9"))
#'
#' @export
plotViromePie <- function(virome = NULL, phylumFilter = NULL, colors = NULL) {

    if (is.null(virome)) {
      stop("Please provide a virome object (see OpenViRome::getVirome)")
    }

    virome <- plottingHelper(virome = virome, phylumFilter = phylumFilter,
                           colors = colors)


    # Get distinct sOTUs, match them to phyla
    sotus <- virome %>%
      dplyr::select(sotu) %>%
      dplyr::distinct() %>%
      dplyr::left_join(virome, by = "sotu", multiple = "first") %>%
      dplyr::group_by(tax_phylum, color) %>%
      dplyr::summarise(count = n())


    p <- plotly::plot_ly(sotus, labels=~tax_phylum, values=~count, type="pie",
                         textposition="inside", sort=FALSE,
                         textinfo="label+percent", hoverinfo="text",
                         text=~paste(tax_phylum, count, sep=": "),
                         marker=list(colors=~color))

    if (!is.null(phylumFilter)) {
      p <- p %>%
        plotly::layout(title="Unique sOTU Count by Predicted Filtered Phylum")
    }
    else {
      p <- p %>% plotly::layout(title="Unique sOTU Count by Predicted Phylum")
    }

    return(p)

}


#' @title drawVirusSankey
#'
#' @description Draws a Sankey diagram of the virome object with the following
#' columns (from left to right): predicted phylum, closest aligned species in
#' GenBank, Serratus sOTU, BioSample accession, metadata-annotated BioSample
#' source species. The width of the lines is proportional to the number of
#' runs for each category. For instance, if an sOTU has 10 outgoing
#' lines, then it appears in 10 different runs. Line colors are determined by
#' predicted phylum.
#' NOTE: for large viromes with many sOTUs, this will take a long time to
#' render and will likely be completely uninterepretable. It is recommended
#' that you filter by a specific viral phylum.
#'
#' @param virome A virome object
#' @param phylumFilter A character vector of viral phyla to filter by. Default
#' is NULL and will not filter the data.
#' @param colors A character vector of colors to use for the Sankey diagram.
#' Default is NULL and will use the Okabe-Ito palette if there are 9 or fewer
#' phyla. Otherwise, uses RColorBrewer Set3. Must be a character vector of hex
#' colors. length(colors) == length(unique(virome$tax_phylum))
#'
#' @return A plotly Sankey diagram of the virome object
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Meloidogyne", con = con)
#' drawVirusSankey(virome = virome,
#'                 phylumFilter = c("Pisuviricota", "Kitrinoviricota"),
#'                 colors = c("#E69F00", "#56B4E9"))
#'
#' @import dplyr
#' @import plotly
#' @import tidyr
#'
#' @export
drawVirusSankey <- function(virome, phylumFilter = NULL, colors = NULL) {

  if (is.null(virome)) {
    stop("Please provide a virome object (see OpenViRome::getVirome)")
  }

  virome <- plottingHelper(virome = virome, phylumFilter = phylumFilter,
                                 colors = colors)

  virome <- virome %>%
    dplyr::select(sotu, tax_phylum, tax_species,
                  bio_project, scientific_name, color)

  # Get counts of each sotu
  sotuCounts <- virome %>%
    dplyr::group_by(sotu) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::arrange(desc(count))

  links <- virome %>%
    dplyr::group_by(tax_phylum, tax_species) %>%
    dplyr::summarise(value = n(), color = first(color)) %>%
    dplyr::rename(source = tax_phylum, target = tax_species) %>%
    dplyr::bind_rows(virome %>%
                       dplyr::group_by(tax_species, sotu, color) %>%
                       dplyr::summarise(value = n(), color = first(color)) %>%
                       dplyr::rename(source = tax_species, target = sotu)) %>%
    dplyr::bind_rows(virome %>%
                       dplyr::group_by(sotu, bio_project, color) %>%
                       dplyr::summarise(value = n(), color = first(color)) %>%
                       dplyr::rename(source = sotu, target = bio_project)) %>%
    dplyr::bind_rows(virome %>%
                       dplyr::group_by(bio_project, scientific_name) %>%
                       dplyr::summarise(value = n(), color = '#dadee3') %>%
                       dplyr::rename(source = bio_project,
                                     target = scientific_name))

  # Create a list of unique nodes
  nodes <- data.frame(name = unique(c(links$source, links$target)))

  # Convert node names to indices
  links$source <- match(links$source, nodes$name) - 1
  links$target <- match(links$target, nodes$name) - 1

  p <- plot_ly(
    type = 'sankey',
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    node = list(
      pad = 10,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      ),
      label = nodes$name,
      color = rep("grey", length(nodes$name)) # Node Coloring
    ),
    link = list(
      source = links$source,
      target = links$target,
      value = links$value,
      color = links$color # Link Coloring
    )
  )


  return(p)
 }

#' @title plottingHelper
#'
#' @description Helper function for the pie chart and sankey diagram functions.
#' Does some sanity checks on user's input filters the data based on the
#' provided phylumFilter.
#' NOT EXPORTED
#'
#' @returns Filtered virome object with an added column for colors assigned to
#' phyla.
#'
#' @import dplyr
#' @import RColorBrewer
#' @import grDevices
#'
#' @tags interal
plottingHelper <- function(virome = NULL, phylumFilter = NULL, colors = NULL) {
  virome <- virome[[1]]

  # Check if phylumFilter is a character vector
  if (!is.null(phylumFilter)) {
    if (!is.character(phylumFilter)) {
      stop("phylumFilter must be a character vector")
    }
    # Filter the virome by phylum
    virome <- virome %>%
      dplyr::filter(tax_phylum %in% phylumFilter)
  }

  uniquePhyla <- unique(virome$tax_phylum)

  # Color sanity checks
  if (!is.null(colors)) {
    if (!is.character(colors)) {
      stop("colors must be a character vector!")
    }
    if (length(colors) != length(uniquePhyla)) {
      stop("colors must be the same length as the number of unique phyla: ",
           length(uniquePhyla))
    }
    tryCatch({
      col2rgb(colors)
    }, error = function(e) {
      stop("Invalid color(s)!")
    })
  } else {
    if (length(uniquePhyla) > 9) {
      print("Warning: more than 9 phyla detected. Discerning colors will be
              difficult. Consider filtering the virome by phylum.")
      colors <- RColorBrewer::brewer.pal(n = length(uniquePhyla),
                                         name = "Set3")
    }
    else {
      colors <- grDevices::palette.colors(n = length(uniquePhyla),
                                          palette = "Okabe-Ito",
                                          alpha = 0.7)
    }
  }
  # Assign colors to phyla
  virome <- virome %>%
    dplyr::mutate(color = colors[match(tax_phylum, uniquePhyla)])

  return(virome)

}

# [END]
