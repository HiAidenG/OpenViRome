#palmPrevalence
#' @title palmPrevalence
#' @description For every sotu in a virome, plot the number of runs in which
#' it was found against the number of bioprojects
#' @param virome A virome object
#' @return A plotly scatterplot
#' @importFrom dplyr select group_by summarise n_distinct
#' @importFrom plotly plot_ly layout add_trace
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales log
#' @importFrom stats cut
#' @export
palmPrevalence <- function(virome) {

  # Collect unique palm_ids, and calculate the mean coverage across all runs
  # Then for each palm_id, sum up the run and bio_project counts
  # Add the gb_pid as a column
  palms <- virome %>%
    dplyr::select(palm_id, run, bio_project, node_coverage, gb_pid, tax_family) %>%
    dplyr::group_by(palm_id) %>%
    dplyr::summarise(runs = n_distinct(run),
                     bio_projects = n_distinct(bio_project),
                     mean_coverage = mean(node_coverage),
                     gb_pid = gb_pid[1],
                     tax_family = tax_family[1])

  # Set a buffer size for the axes
  buffer_size <- 0.5

  # Adjust axis ranges to include dynamic buffer
  x_range <- c(min(palms$runs) - buffer_size, max(palms$runs) + buffer_size)
  y_range <- c(min(palms$bio_projects) - buffer_size, max(palms$bio_projects) + buffer_size)

  # Scale the mean coverage
  palms$mean_coverage <- log(palms$mean_coverage)

  # Add jitter to the data to avoid overplotting
  set.seed(123) # Setting seed for reproducibility
  jitter_amount <- 0.7
  palms$runs_jitter <- jitter(palms$runs, amount=jitter_amount)
  palms$bio_projects_jitter <- jitter(palms$bio_projects, amount=jitter_amount)

  bins <- c(0, 40, 60, 80, 90, 100)
  color_labels <- c("0-40", "40-60", "60-80", "80-90", "90-100")
  colors <- c("#ed2939", "#ff7900", '#ffc40c', '#fada5e', '#708090')

  # Create a named vector for colors
  named_colors <- setNames(colors, color_labels)

  # Use cut to create a factor for color assignment
  palms$color <- cut(palms$gb_pid, bins, labels=color_labels)

  # Plot
  plot <- plotly::plot_ly(palms, x=~runs_jitter, y=~bio_projects_jitter, color=~color,
                          size=~mean_coverage, colors = named_colors, alpha=0.9,
                          type = 'scatter', mode = 'markers',
                          text = ~paste("Runs: ", runs,
                                        "<br>Bioprojects: ", bio_projects,
                                        "<br>palm_id: ", palm_id,
                                        '<br>log(Mean coverage): ',
                                        round(mean_coverage, digits=2),
                                        '<br>Genbank identity: ', gb_pid, '%'),
                          marker=list(sizemode = 'diameter',
                                      sizeref = 5,
                                      sizemin = 2)) %>%
    plotly::layout(title='Palm Prevalence',
                   xaxis = list(title = 'Number of Runs',
                                tickmode = 'logarithmic',
                                automargin = TRUE),
                   yaxis = list(title = 'Number of BioProjects',
                                tickmode = 'linear',
                                automargin = TRUE),
                   legend = list(title = list(text = 'GenBank ID (%)')))

  # Return the plot object
  return(plot)
}
