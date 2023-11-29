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
#' @import
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


createSankeyDiagram <- function(virome) {
  # Create connections between the nodes
  data <- virome[[1]]
  # Create connections between the nodes
  links <- data %>%
    count(tax_phylum, tax_species) %>%
    rename(source = tax_phylum, target = tax_species, value = n) %>%
    bind_rows(data %>%
                count(tax_species, sotu) %>%
                rename(source = tax_species, target = sotu, value = n)) %>%
    bind_rows(data %>%
                count(sotu, bio_project) %>%
                rename(source = sotu, target = bio_project, value = n)) %>%
    bind_rows(data %>%
                count(bio_project, scientific_name) %>%
                rename(source = bio_project, target = scientific_name, value = n))

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
      label = nodes$name
    ),
    link = list(
      source = links$source,
      target = links$target,
      value = links$value
    )
  )

  return(p)
}









