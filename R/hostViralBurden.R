#' @title getHostViralBurden
#' @description Takes a virome object as input and re-queries Serratus to get
#' the incidence of viral RNA positive SRA runs out of the total number queried
#' for each source species.
#' @param virome A virome object
#' @return A data frame with the total number of runs queried and the number
#' of runs that were virus positive for each source species.
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Meloidogyne", con = con)
#' getHostViralBurden(virome = virome, con = con)
#' @import dplyr
getHostViralBurden <- function(virome = NULL) {

  if (is.null(virome)) {
    stop("Please provide a virome object (see OpenViRome::getVirome)")
  }

  vbDF <- virome[[2]]
  vbDF <- vbDF %>%
    dplyr::group_by(scientific_name) %>%
    dplyr::summarise(total = n(),
                     virus_positive = sum(virus_positive)) %>%
    dplyr::mutate(percent_virus_positive = virus_positive / total)

  return(vbDF)
}

#' @title plotVirusPositive
#' @description Plots the proportion of virus positive runs for each taxon
#' @param vbDF
#' @return A plotly bar chart
#' @importFrom dplyr group_by %>% summarise collect left_join
#' @importFrom tidyr pivot_longer
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' vbDF <- getHostViralBurden(virome = virome)
#' plotVirusPositive(vbDF)
plotVirusPositive <- function(vbDF = NULL) {
  # Add a column for virus negative counts
  vbDF <- vbDF %>%
    dplyr::mutate(virus_negative = total - virus_positive)

  # Sort the dataframe by virus_positive count in descending order and convert scientific_name to factor
  vbDF <- vbDF %>%
    dplyr::arrange(desc(virus_positive)) %>%
    dplyr::mutate(scientific_name = factor(scientific_name,
                                           levels = unique(scientific_name)))

  # Plot a stacked horizontal bar chart
  p <- plot_ly(data = vbDF, x = ~virus_positive, y = ~scientific_name,
               type = 'bar',
               name = 'Virus Positive', marker = list(color = 'cornflowerblue'),
               orientation = 'h') %>%
    add_trace(x = ~virus_negative, name = 'Virus Negative',
              marker = list(color = 'grey'), orientation = 'h') %>%
    layout(barmode = 'stack', yaxis = list(title = 'Run Source Species'),
           xaxis = list(title = 'Runs'), showlegend = TRUE) %>%
    plotly::config(displayModeBar = FALSE)

  return(p)
}


#' @title drawVirusSankey
#' @description Draws a Sankey diagram of the virome object. NOTE: for
#' large viromes (i.e. >200 rows) this will take a long time to render and be
#' completely uninterpretable. It is recommended that you filter by a specific
#' viral phylum.
#' @param virome A virome object
#' @param phylumFilter A character vector of viral phyla to filter by. Default
#' is NULL and will not filter the data.
#' @param abundanceFilter Filter sOTUs with abundance across all BioSamples
#' less than this value. Default is 0.
#' @return A plotly Sankey diagram
#' @examples
#' virome <- getVirome(tax = "Meloidogyne", con = con)
#' drawVirusSankey(virome = virome)
#' @import dplyr
#' @import plotly
#' @import tidyr
#' @import RColorBrewer
#' @export
drawVirusSankey <- function(virome, phylumFilter = NULL, abundanceFilter = 0) {
  # Extract the data from the virome list
  data <- virome[[1]]

  # Assign 'unknown' to missing tax_phylum
  data$tax_phylum[is.na(data$tax_phylum)] <- "Unknown"

  # Filter the data based on the provided tax_phylum, if specified
  if (!is.null(phylumFilter)) {
    data <- data[data$tax_phylum %in% phylumFilter, ]
  }

  # Filter the data based on the provided abundanceFilter, if specified
  sotuCounts <- data %>%
    group_by(sotu) %>%
    summarise(abundance = n()) %>%
    filter(abundance >= abundanceFilter) %>%
    select(sotu)

  data <- data[data$sotu %in% sotuCounts$sotu, ]

  # Create the color mapping for each tax_phylum
  unique_phylum <- unique(data$tax_phylum)
  colors <- grDevices::palette.colors(n = length(unique_phylum), alpha = 0.7)
  color_mapping <- setNames(colors, unique_phylum)

  # Apply color mapping to each tax_phylum in the data
  data$color <- color_mapping[data$tax_phylum]

  # Create links for the Sankey diagram
  links <- data %>%
    group_by(tax_phylum, tax_species) %>%
    summarise(value = n(), color = first(color)) %>%
    rename(source = tax_phylum, target = tax_species) %>%
    bind_rows(data %>%
                group_by(tax_species, sotu) %>%
                summarise(value = n(), color = first(color)) %>%
                rename(source = tax_species, target = sotu)) %>%
    bind_rows(data %>%
                group_by(sotu, bio_project) %>%
                summarise(value = n(), color = first(color)) %>%
                rename(source = sotu, target = bio_project)) %>%
    bind_rows(data %>%
                group_by(bio_project, scientific_name) %>%
                summarise(value = n(), color = "#ebeced") %>%
                rename(source = bio_project, target = scientific_name))

  # Create a list of unique nodes
  nodes <- data.frame(name = unique(c(links$source, links$target)))


  # Convert node names to indices
  links$source <- match(links$source, nodes$name) - 1
  links$target <- match(links$target, nodes$name) - 1

  # Create Sankey diagram
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
      color = rep("grey", length(nodes$name))  # All nodes are grey
    ),
    link = list(
      source = links$source,
      target = links$target,
      value = links$value,
      color = links$color  # Link colors based on originating tax_phylum
    )
  )

  return(p)
}

















