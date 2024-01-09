# palmPrevalence.R

#' @title palmPrevalence
#'
#' @description For every sotu in a virome, plot the number of runs in which
#' it was found against the number of bioprojects. Points are sized by
#' normalized mean coverage. Colors are genbank divergence.
#'
#' @param virome A virome object
#'
#' @return A plotly scatterplot
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome("Meloidogyne", con = con)
#' palmPrevalence(virome)
#'
#'
#' @importFrom dplyr select group_by summarise n_distinct
#' @importFrom plotly plot_ly layout add_trace
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
palmPrevalence <- function(virome) {
  virome <- virome[[1]]

  # Collect unique sotus, and calculate the mean coverage across all runs
  # Then for each sotu, sum up the run and bio_project counts
  # Add the gb_pid as a column
  palms <- virome %>%
    dplyr::select(sotu, run, bio_project, node_coverage_norm, gb_pid,
                  tax_family) %>%
    dplyr::group_by(sotu) %>%
    dplyr::summarise(runs = n_distinct(run),
                     bio_projects = n_distinct(bio_project),
                     mean_coverage = mean(node_coverage_norm),
                     gb_pid = gb_pid[1],
                     tax_family = tax_family[1])

  # Scale the mean coverage
  palms$mean_coverage <- log(palms$mean_coverage)

  # Add jitter to the data to avoid overplotting
  set.seed(123) # Setting seed for reproducibility
  jitter_amount <- 0.4

  # Custom jitter function that ensures values do not go below 1
  safe_jitter <- function(x, amount) {
    jittered <- jitter(x, amount = amount)
    pmax(jittered, 1)  # Ensuring that the jittered values do not fall below 1
  }

  palms$runs_jitter <- safe_jitter(palms$runs, amount = jitter_amount)
  palms$bio_projects_jitter <- safe_jitter(palms$bio_projects, amount = jitter_amount)

  bins <- c(0, 40, 60, 80, 90, 100)
  color_labels <- c("0-40", "40-60", "60-80", "80-90", "90-100")
  colors <- c("#ed2939", "#ff7900", '#ffc40c', '#fada5e', '#708090')

  # Create a named vector for colors
  named_colors <- setNames(colors, color_labels)

  # Use cut to create a factor for color assignment
  palms$color <- cut(palms$gb_pid, bins, labels=color_labels)

  # Plot
  plot <- plotly::plot_ly(palms, x=~runs_jitter, y=~bio_projects_jitter, color=~color,
                          size=~mean_coverage, colors = named_colors, alpha=0.5,
                          type = 'scatter', mode = 'markers',
                          text = ~paste("Runs: ", runs,
                                        "<br>Bioprojects: ", bio_projects,
                                        "<br>sotu: ", sotu,
                                        '<br>log(Normalized Mean Coverage): ',
                                        round(mean_coverage, digits=2),
                                        '<br>Genbank identity: ', gb_pid, '%'),
                          marker=list(sizemode = 'diameter',
                                      sizeref = 5,
                                      sizemin = 2)) %>%
    plotly::layout(title='sOTU Virome Prevalence',
                   xaxis = list(title = 'Runs',
                                type = 'log',  # Set x-axis to logarithmic scale
                                automargin = FALSE,
                                range = c(log10(1-0.1), log10(max(palms$runs_jitter)+10))),
                   yaxis = list(title = 'BioProjects',
                                automargin = FALSE,
                                type = 'log',
                                range = c(log10(1-0.1), log10(max(palms$bio_projects_jitter)+1))),
                   legend = list(title = list(text = 'GenBank ID (%)')))
  # Return the plot object
  return(plot)
}

#' @title palmPrevalenceKernels
#'
#' @description For each of the four distributions present in the palmPrevalence
#' plot (SRA runs, bio_projects, mean coverage, and genbank identity), plot the
#' kernel density estimate.
#'
#' @param virome A virome object
#'
#' @return A plotly plot
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome("Meloidogyne", con = con)
#' palmPrevalenceKernels(virome)
#'
#' @import dplyr ggplot2
#'
#' @export
palmPrevalenceKernels <- function(virome) {
  virome <- virome[[1]]

  runs <- virome %>% select(sotu, run) %>% group_by(sotu) %>% summarise(runs = n_distinct(run))
  bio_projects <- virome %>% select(sotu, bio_project) %>% group_by(sotu) %>% summarise(bio_projects = n_distinct(bio_project))
  mean_coverage <- virome %>% select(sotu, node_coverage_norm) %>% group_by(sotu) %>% summarise(mean_coverage = mean(node_coverage_norm))
  gb_pid <- virome %>% select(sotu, gb_pid) %>% group_by(sotu) %>% summarise(gb_pid = gb_pid[1])

  sotuTable <- runs %>% left_join(bio_projects, by = 'sotu') %>% left_join(mean_coverage, by = 'sotu') %>% left_join(gb_pid, by = 'sotu')

  # Plot each variable on a separate facet
  plot <- sotuTable %>%
    gather(variable, value, -sotu) %>%
    ggplot(aes(x = value)) +
    geom_density() +
    geom_rug() +
    facet_wrap(~variable, scales = 'free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, title = 'Kernel Density Estimates for sOTU Prevalence Variables')

  return(plot)
}


# [END]
