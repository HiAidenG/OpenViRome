#' @title getShannonDiversity
#' @description Calculate Shannon diversity of a virome object
#' @param virome A virome object
#' @param sourceSpecies Optional argument for filtering by source species
#' @return A numeric value (Shannon Diversity Index)
#' @importFrom dplyr filter group_by %>% summarise mutate
#' @export
getShannonDiversity <- function(virome = NULL, sourceSpecies = NULL) {
  # Check if a source species was provided
  if (!is.null(sourceSpecies)) {

    #Check that the source species is valid
    if (!sourceSpecies %in% virome$scientific_name) {
      stop("Error: source species not found in virome object")
    }

    virome <- virome %>%
      filter(scientific_name == sourceSpecies)

  }
    # Calculate frequencies
    freq <- table(virome$tax_family) / length(virome$tax_family)

    # Calculate Shannon diversity
    shannon <- -sum(freq * log(freq))

    return(shannon)
  }


#' @title getSimpsonDiversity
#' @description Calculate Simpson diversity of a virome object
#' @param virome A virome object
#' @param sourceSpecies Optional argument for filtering by source species
#' @return A numeric value (Simpson Diversity Index)
#' @importFrom dplyr filter group_by %>% summarise mutate
#' @export
getSimpsonDiversity <- function(virome = NULL, sourceSpecies = NULL) {
  # Check if a source species was provided
  if (!is.null(sourceSpecies)) {

    # Check that the source species is valid
    if (!sourceSpecies %in% virome$scientific_name) {
     stop("Error: source species not found in virome object")
   }

    virome <- virome %>%
      filter(scientific_name == sourceSpecies)

  }
  # Calculate Simpson diversity
  simpson <- virome %>%
    group_by(tax_family) %>%
    summarise(count = n()) %>%
    mutate(freq = count / sum(count)) %>%
    mutate(simpson = freq^2) %>%
    summarise(simpson = sum(simpson))

  return(simpson[[1]])
}

#' @title getSpeciesRichness
#' @description Calculate species richness of a virome object
#' @param virome A virome object
#' @param sourceSpecies Optional argument for filtering by source species
#' @return A numeric value (Species Richness)
#' @importFrom dplyr filter group_by %>% summarise mutate
#' @export
getSpeciesRichness <- function(virome = NULL, sourceSpecies = NULL) {
  # Check if a source species was provided
  if (!is.null(sourceSpecies)) {

    #Check that the source species is valid
    if (!sourceSpecies %in% virome$scientific_name) {
      stop("Error: source species not found in virome object")
    }

    virome <- virome %>%
      filter(scientific_name == sourceSpecies)

  }
  # Calculate species richness
  speciesRichness <- virome %>%
    group_by(tax_family) %>%
    summarise(count = n()) %>%
    summarise(speciesRichness = n())

  return(speciesRichness[[1]])
}

