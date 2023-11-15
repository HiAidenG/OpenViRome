#incidencePlot
#' @title incidencePlot
#' @description Plot the incidence of virus positive sra runs against all
#' runs analyzed
#' @param virome A virome object
#' @param con A database connection
#' @return A plotly bar chart
#' @importFrom dplyr group_by %>% summarise collect left_join pivot_longer
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @export
incidencePlot <- function(virome, con) {

  tax <- virome$scientific_name %>% unique()
  searchString <- paste(tax, collapse = "|")

  virome <- virome %>%
    dplyr::group_by(scientific_name) %>%
    dplyr::summarise(virus_positive=n()) %>%
    collect()

  sra <- tbl(con, 'srarun') %>%
    dplyr::filter(grepl(pattern = searchString, x = scientific_name,
                        ignore.case = TRUE)) %>%
    dplyr::group_by(scientific_name) %>%
    dplyr::summarise(total=n()) %>%
    collect()

  # join the two tables
  result <- dplyr::left_join(sra, virome,
                             by=c('scientific_name'='scientific_name'))

  long_result <- result %>%
    pivot_longer(cols = c(virus_positive, total),
                 names_to = "category",
                 values_to = "count")


  # Plot the transformed data
  plot <- ggplot(long_result, aes(x = scientific_name, y = count, fill = category)) +
    geom_bar(stat = 'identity', position = 'stack') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = 'Species', y = 'SRA Runs') +
    ggtitle('Incidence of vRNA Positive Runs by Source Species') +
    scale_fill_manual(values = c('virus_positive' = '#5e23f7', 'total' = '#d3d3db'),
                      labels = c("Total Runs Searched", "vRNA Positive Runs")) +
    theme_light()

  print(plot)

}

