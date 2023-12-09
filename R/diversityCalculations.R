# diversityCalculations.R

#' @title getAlphaDiveristy
#'
#' @description A wrapper function for calculating various metrics of alpha
#' diversity for a virome object. Options to calculate Shannon diversity,
#' Simpson diversity, species richness, and Pielou evenness. Each bioSample
#' is treated as an independent sample for calculating these metrics. A more
#' detailed explanation in the context of SRA virome data is available in the
#' vignette.
#'
#' @param virome A virome object
#' @param bioSample Optional argument for filtering the virome by BioSample
#' @param mode A character vector of diversity metrics to calculate. Options:
#' "shannon", "simpson", "richness", "evenness". Default is "shannon".
#'
#' @return A tibble with columns for each diversity metric calculated.
#'
#' @import dplyr
#'
#' @examples
#' # Calculate Shannon diversity for all virome data
#' data("TylenchoideaVirome")
#' getAlphaDiversity(TylenchoideaVirome) # Shannon diversity by default
#'
#' # Calculate Shannon diversity for a specific biosample
#' getAlphaDiversity(TylenchoideaVirome, bioSample = "SAMN25258043")
#'
#' # Calculate Simpson diversity for all virome data
#' getAlphaDiversity(TylenchoideaVirome, mode = "simpson")
#'
#' # Calculate all metrics for all virome data
#' getAlphaDiversity(TylenchoideaVirome, mode = c("shannon", "simpson",
#'                   "richness", "evenness"))
#'
#' @export
getAlphaDiversity <- function(virome = NULL, mode = "shannon",
                              bioSample = NULL) {
  if (is.null(virome)) {
    stop("Please provide a virome object")
  }

  # Check if a BioSample was provided
  if (!is.null(bioSample)) {

    # Check that the BioSample is valid
    if (!bioSample %in% virome$BioSample) {
      stop("Error: BioSample not found in virome object")
    }

    virome <- virome %>%
      filter(bio_sample == bioSample)

  }

  # Check if mode is a character vector
  if (!is.character(mode)) {
    stop("Error: mode must be a character vector")
  }

  # Check if mode is a valid option
  if (!all(mode %in% c("shannon", "simpson", "richness", "evenness"))) {
    stop("Error: mode must be one of 'shannon', 'simpson', 'richness', or
         'evenness'")
  }


  returnTable <- tibble(bio_sample = unique(virome$bio_sample))
  for (i in 1:length(mode)) {
    if (mode[i] == "shannon") {
      returnTable <- returnTable %>%
        full_join(getDiversity(virome, mode="shannon"),
                  by = "bio_sample")
    }
    else if (mode[i] == "simpson") {
      returnTable <- returnTable %>%
        full_join(getDiversity(virome, mode="simpson"),
                  by = "bio_sample")
    }
    else if (mode[i] == "richness") {
      returnTable <- returnTable %>%
        full_join(getRichness(virome), by = "bio_sample")
    }
    else if (mode[i] == "evenness") {
      returnTable <- returnTable %>%
        full_join(getEvenness(virome), by = "bio_sample")
    }
  }

  return(returnTable)

}

#' @title getRichness
#'
#' @description Calculate species richness (i.e. number of unique viral
#' families) of a virome object.
#' NOT EXPORTED
#'
#' @param virome A virome object
#'
#' @keywords internal
#'
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for richness.
#'
#' @import dplyr
getRichness <- function(virome) {

  virome <- virome[[1]]

  richness <- virome %>%
    group_by(bio_sample) %>%
    summarise(richness = n_distinct(sotu)) %>%
    left_join(virome, by = "bio_sample") %>%
    select(bio_sample, scientific_name, richness) %>%
    distinct()

  return(richness)

}

#' @title getDiversity
#'
#' @description Calculate either Shannon or Simpson diversity of a virome object
#' Note: There are multiple metrics referred to in the literature as "Simpson
#' diveristy". This function calculates the D_0, which is the probability that
#' two randomly selected reads will be of the same sotu.
#' NOT EXPORTED
#'
#' @param virome A virome object
#' @param mode Either "shannon" or "simpson"
#' @param con A database connection object (to Serratus SQL).
#'
#' @keywords internal
#'
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for type of diversity calculated. Each bio_sample is treated as a
#' sample, and each sotu is treated as a species. Reads mapping to each sotu
#' are treated as the abundance of that species.
#'
#' @import dplyr
getDiversity <- function(virome, mode = "shannon") {
  virome <- virome[[1]]

  # Calculate proportion of coverage in each biosample
  virome <- virome %>%
    group_by(bio_sample) %>%
    mutate(total_coverage = sum(node_coverage_norm)) %>%
    mutate(prop = node_coverage_norm / total_coverage) %>%
    select(bio_sample, sotu, prop, scientific_name)

  if (mode == "shannon") {
    diversity <- virome %>%
      group_by(bio_sample) %>%
      summarise(shannon = -sum(prop * log(prop)), .groups = 'drop')
  } else if (mode == "simpson") {
    diversity <- virome %>%
      group_by(bio_sample) %>%
      summarise(simpson = sum(prop^2), .groups = 'drop')

  } else {
    stop("Invalid mode specified. Use 'shannon' or 'simpson'.")
  }

  # Add scientific name
  diversity <- diversity %>%
    left_join(virome %>% select(bio_sample, scientific_name) %>% distinct(),
              by = "bio_sample")

  return(diversity)
}

#' @title getEvenness
#'
#' @description Calculate Pielou evenness of a virome object.
#' Note: Pielou evenness is defined as the Shannon diversity divided by the
#' natural log of the richness. Evenness metrics have known issues as they
#' are largely dependent on species (in this case reads) abundance, which
#' can be variable between samples (due to sequencing depth). Some values
#' may be NaN if the richness is 1.
#' NOT EXPORTED
#'
#' @param virome A virome object
#'
#' @keywords internal
#'
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for evenness.
#'
#' @import dplyr
#'
#' @export
getEvenness <- function(virome) {


  # Calculate Shannon diversity
  shannon <- getDiversity(virome, mode = "shannon")

  # Calculate richness
  richness <- getRichness(virome)
  richness <- richness %>% select(-scientific_name)
  # Join the two tables
  evenness <- shannon %>%
    full_join(richness, by = "bio_sample")

  # Calculate evenness
  evenness <- evenness %>%
    mutate(evenness = shannon / log(richness))

  # Make NA and Infinite values 0
  evenness[is.na(evenness)] <- 0

  evenness <- evenness %>% select(bio_sample, scientific_name, evenness)


  return(evenness)

}

#' @title viromeToOTUTable
#'
#' @description Convert a virome object to the ampvis2 otutable format where
#' each row is a viral sOTU and each column is a biosample. The last two
#' columns provide the tax_species and tax_phylum for each sOTU.
#'
#' @param virome A virome object
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Salidae", con = con)
#' otuTable <- viromeToOTUTable(virome)
#' head(otuTable)
#'
#' @import dplyr tidyr
#'
#' @export
viromeToOTUTable <- function(virome) {
  data <- virome[[1]]

  # Get the sOTU table
  sotuTable <- data %>%
    group_by(bio_sample, sotu) %>%
    summarise(count = sum(node_coverage_norm)) %>%
    pivot_wider(names_from = bio_sample, values_from = count, values_fill = 0)

  # Add the taxonomic information
  sotuTable <- sotuTable %>%
    left_join(data %>% select(sotu, tax_phylum) %>% distinct(), by = "sotu",
              multiple = "first") %>%
    left_join(data %>% select(sotu, tax_species) %>% distinct(), by = "sotu",
              multiple = "first")

  # Rename so ampvis will recognize
  sotuTable <- sotuTable %>%
    rename(phylum = tax_phylum, species = tax_species, otu=sotu)

  return(sotuTable)
}

#' @title viromeToMetadataTable
#'
#' @description Convert a virome object to the metadata table format specified
#' in ampvis2. First column specifies biosample id's, second column specifies
#' bioproject id's, third column is the scientific_name of the sample.
#'
#' @import dplyr
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Salidae", con = con)
#' metadata <- viromeToMetadataTable(virome)
#' head(metadata)
#'
#' @export
viromeToMetadataTable <- function(virome = NULL) {
  if (is.null(virome)) {
    stop("Must provide a virome object.")
  }

  data <- virome[[1]]

  # Get the metadata table
  metadata <- data %>%
    select(bio_sample, scientific_name, bio_project) %>%
    distinct()

  return(metadata)
}


# TODO: This needs to be fixed, breaks when plotting many genera.
#' @title plotAlphaDiversity
#'
#' @description Automatically plots Simpson, Shannon, and Pielou evenness for
#' a virome object.
#'
#' @param virome A virome object
#'
#' @param mode Whether to calculate Shannon or Simpson diversity.
#' Default = 'shannon'.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Salidae", con = con)
#' plotAlphaDiversity(virome, mode = "simpson")
#'
#' @export
plotAlphaDiversity <- function(virome, mode = "shannon") {
  # Get alpha diversity
  alphaDiversity <- getDiversity(virome, mode = mode)

  # Get evenness
  evenness <- getEvenness(virome)
  evenness <- evenness %>% select(-scientific_name)

  # Join the tables
  alphaDiversity <- alphaDiversity %>%
    left_join(evenness, by = "bio_sample")

  # Get distinct genera
  alphaDiversity <- alphaDiversity %>%
    mutate(genus = gsub(" .*", "", scientific_name))
  genera <- alphaDiversity %>% select(genus) %>% distinct()

  # Set colors for each genus
  colors <- grDevices::palette.colors(n = nrow(genera), alpha=0.5)

  # Assign colors to each genus
  alphaDiversity <- alphaDiversity %>%
    left_join(genera, by = "genus") %>%
    mutate(color = colors)

  alphaDiversity$hover_info <- paste("Biosample: ", alphaDiversity$bio_sample)

  # Plot
    if (mode == "shannon") {
        p1 <- plot_ly(alphaDiversity, x = ~scientific_name, y = ~shannon,
                      color = ~color, type = "box",
                      boxpoints = "all", jitter = 0.3, pointpos = 0.0,
                      line = list(width = 0),
                      fillcolor = "rgba(0,0,0,0)",
                      name = mode,
                      text = ~hover_info,  # Add hover text
                      hoverinfo = "text+y") %>%
            layout(yaxis = list(title = "Shannon Diversity"),
                   xaxis = list(title = "Source Species"))
    } else if (mode == "simpson") {
        p1 <- plot_ly(alphaDiversity, x = ~scientific_name, y = ~simpson,
                      color = ~color, type = "box",
                      boxpoints = "all", jitter = 0.3, pointpos = 0.0,
                      line = list(width = 0),
                      fillcolor = "rgba(0,0,0,0)",
                      name = mode,
                      text = ~hover_info,  # Add hover text
                      hoverinfo = "text+y") %>%
            layout(yaxis = list(title = "Simpson Diversity"),
                   xaxis = list(title = "Source Species"))
    }

    return(p1)
}

#' @title plotBetaDiversity
#'
#' @description Automatically plots beta diversity for a virome object using
#' ampvis2. This is the plot used for the shiny frontend, and uses Bray-Curtis
#' dissimilarity with PCOA. There are other options for plotting beta diversity,
#' see https://kasperskytte.github.io/ampvis2/articles/ampvis2.html#ordination.
#'
#' @param virome A virome object
#'
#' @import ampvis2
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax = "Salidae", con = con)
#' plotBetaDiversity(virome)
#'
#' @export
plotBetaDiversity <- function(virome = NULL) {

  if (is.null(virome)) {
    stop("Must provide a virome object.")
  }

  # Get the sOTU table
  sotuTable <- viromeToOTUTable(virome)

  # Get the metadata table
  metadata <- viromeToMetadataTable(virome)

  # Plot
  d <- ampvis2::amp_load(otutable = sotuTable, metadata = metadata)
  p1 <- ampvis2::amp_ordinate(d, type = "PCOA", distmeasure = "bray",
                               transform = "none",
                               sample_color_by = "scientific_name",
                               sample_plotly = "all",
                               filter_species = 0.0
                               )

  return(p1)
}

# [END]
