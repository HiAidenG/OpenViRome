#' @title getAlphaDiveristy
#' @description A wrapper function for calculating various metrics of alpha
#' diversity for a virome object. Options to calculate Shannon diversity,
#' Simpson diversity, species richness, and Pielou evenness. Each bioSample
#' is treated as an independent sample for calculating these metrics. A more
#' detailed explanation in the context of SRA virome data is available in the
#' vignette.
#' @param virome A virome object
#' @param bioSample Optional argument for filtering the virome by BioSample
#' @param mode A character vector of diversity metrics to calculate. Options:
#' "shannon", "simpson", "richness", "evenness". Default is "shannon".
#' @return A tibble with columns for each diversity metric calculated.
#' @importfrom dplyr filter group_by %>% summarise mutate
#' @export
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

  # connect to Serratus
  con <- palmid::SerratusConnect()

  returnTable <- tibble(bio_sample = unique(virome$bio_sample))
  for (i in 1:length(mode)) {
    if (mode[i] == "shannon") {
      returnTable <- returnTable %>%
        full_join(getDiversity(virome, mode="shannon", con=con),
                  by = "bio_sample")
    }
    else if (mode[i] == "simpson") {
      returnTable <- returnTable %>%
        full_join(getDiversity(virome, mode="simpson", con=con),
                  by = "bio_sample")
    }
    else if (mode[i] == "richness") {
      returnTable <- returnTable %>%
        full_join(getRichness(virome), by = "bio_sample")
    }
    else if (mode[i] == "evenness") {
      returnTable <- returnTable %>%
        full_join(getEvenness(virome), by = "bio_sample", con=con)
    }
  }

  return(returnTable)

}

#' @title getRichness
#' @description Calculate species richness (i.e. number of unique viral
#' families) of a virome object.
#' @param virome A virome object
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for richness.
#' @importFrom dplyr filter group_by %>% summarise mutate
#' NOT EXPORTED
getRichness <- function(virome) {

  richness <- virome %>%
    group_by(bio_sample) %>%
    summarise(richness = n_distinct(sotu))

  return(richness)

}

#' @title getDiversity
#' @description Calculate either Shannon or Simpson diversity of a virome object
#' Note: There are multiple metrics referred to in the literature as "Simpson
#' diveristy". This function calculates the D_0, which is the probability that
#' two randomly selected reads will be of the same sotu.
#' @param virome A virome object
#' @param mode Either "shannon" or "simpson"
#' @param con A database connection object (to Serratus SQL).
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for type of diversity calculated. Each bio_sample is treated as a
#' sample, and each sotu is treated as a species. Reads mapping to each sotu
#' are treated as the abundance of that species.
#' @importFrom dplyr filter group_by %>% summarise mutate
#' NOT EXPORTED
getDiversity <- function(virome, mode = "shannon", con) {
  # Calculate the sum of node_coverage for each sotu within each bio_sample
  virome <- virome %>%
    group_by(bio_sample, sotu) %>%
    summarise(node_coverage_sum = sum(node_coverage), .groups = 'drop')

  # Calculate the total coverage and proportion for each sotu in each bio_sample
  virome <- virome %>%
    group_by(bio_sample) %>%
    mutate(total_coverage = sum(node_coverage_sum)) %>%
    mutate(prop = node_coverage_sum / total_coverage) %>%
    select(bio_sample, sotu, prop)

  # Get the library size for normalization
  bioSamples <- unique(virome$bio_sample)
  librarySize <- tbl(con, "srarun") %>%
    filter(bio_sample %in% bioSamples) %>%
    group_by(bio_sample) %>%
    summarise(library_size = sum(spots), .groups = 'drop') %>%
    select(bio_sample, library_size) %>%
    collect()

  # Make the numbers a bit nicer, scaling by 1e8
  librarySize$library_size <- librarySize$library_size / 1e8

  # Scale proportion by library size
  virome <- virome %>%
    left_join(librarySize, by = "bio_sample") %>%
    mutate(prop = prop / library_size)

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

  return(diversity)
}

#' @title getEvenness
#' @description Calculate Pielou evenness of a virome object.
#' Note: Pielou evenness is defined as the Shannon diversity divided by the
#' natural log of the richness. Evenness metrics have known issues as they
#' are largely dependent on species (in this case reads) abundance, which
#' can be variable between samples (due to sequencing depth). Some values
#' may be NaN if the richness is 1.
#' @param virome A virome object
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for evenness.
#' @importFrom dplyr filter group_by %>% summarise mutate
#' NOT EXPORTED
getEvenness <- function(virome, con) {

  # Calculate Shannon diversity
  shannon <- getDiversity(virome, mode = "shannon", con=con)

  # Calculate richness
  richness <- getRichness(virome)

  # Join the two tables
  evenness <- shannon %>%
    full_join(richness, by = "bio_sample")

  # Calculate evenness
  evenness <- evenness %>%
    mutate(evenness = shannon / log(richness))

  # only want to return the evenness column
  evenness <- evenness %>%
    select(bio_sample, evenness)

  return(evenness)

}

# TODO: Beta-diversity (Bray-curtis)

# [END]






