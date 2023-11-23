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
#' @import dplyr
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
#' @description Calculate species richness (i.e. number of unique viral
#' families) of a virome object.
#' NOT EXPORTED
#' @param virome A virome object
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for richness.
#' @import dplyr
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
#' NOT EXPORTED
#' @param virome A virome object
#' @param mode Either "shannon" or "simpson"
#' @param con A database connection object (to Serratus SQL).
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for type of diversity calculated. Each bio_sample is treated as a
#' sample, and each sotu is treated as a species. Reads mapping to each sotu
#' are treated as the abundance of that species.
#' @importFrom dplyr filter group_by %>% summarise mutate
getDiversity <- function(virome, mode = "shannon") {

  con <- palmid::SerratusConnect()

  virome <- getSpeciesCounts(virome, con) # get normalized counts for each sotu

  # Calculate proportion of coverage in each biosample
  virome <- virome %>%
    group_by(bio_sample) %>%
    mutate(total_coverage = sum(node_coverage_sum)) %>%
    mutate(prop = node_coverage_sum / total_coverage) %>%
    select(bio_sample, sotu, prop)

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
#' NOT EXPORTED
#' @param virome A virome object
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for evenness.
#' @importFrom dplyr filter group_by %>% summarise mutate
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

#' @title getSpeciesCounts
#' @description Calculate the number of reads mapping to each sotu in a virome
#' object.
#' NOT EXPORTED
#' @param virome A virome object
#' @param con A database connection object (to Serratus SQL).
#' @keywords internal
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for each sotu, with the normalized number of reads mapping to that
#' sotu.
getSpeciesCounts <- function(virome = NULL, con = NULL) {

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

  # Calculate the sum of node_coverage for each sotu within each bio_sample
  virome <- virome %>%
    group_by(bio_sample, sotu) %>%
    summarise(node_coverage_sum = sum(node_coverage), .groups = 'drop')

  # Scale node_coverage_sum by library size
  virome <- virome %>%
    left_join(librarySize, by = "bio_sample") %>%
    mutate(node_coverage_sum = node_coverage_sum / library_size) %>%
    select(bio_sample, sotu, node_coverage_sum)

  return(virome)

}

#' @title getBetaDiversity
#' @description Calculate pairwise Bray-curtis dissimilarity of a virome object
#' for all bioSamples present in the virome. Bray-curtis dissimilarity is
#' defined as sum(|x_ij - x_ik|) / sum(x_ij + x_ik), where where x_ij and
#' x_ik are the number of reads mapping to sotu i and bioSample j and k,
#' respectively.
#' @param virome A virome object
#' @return A tibble with rows for each bioSample present in the virome and a
#' column for each bioSample, with the Bray-curtis dissimilarity between the
#' two bioSamples.
#' @export
getBetaDiversity <- function(virome) {

  con <- palmid::SerratusConnect()

  virome <- getSpeciesCounts(virome, con)

  # Calculate the pairwise Bray-curtis dissimilarity
  betaDiversity <- virome %>%
    pivot_wider(names_from = sotu, values_from = prop,
                values_fill = 0) %>%
    column_to_rownames("bio_sample") %>%
    as.matrix() %>%
    vegdist(method = "bray", na.rm=TRUE)

  # Convert to a tibble
  betaDiversity <- as_tibble(as.matrix(betaDiversity))
  rownames(betaDiversity) <- colnames(betaDiversity)
  return(betaDiversity)
}

#' @title getAmpvisData
#' @description Return an OTU-table and metadata table in ampvis format
#' (OTUs as rows, samples as columns), values are the normalized number of
#' reads mapping to each OTU in each sample. Metadata table is a subset of the
#' palm_virome table. Columns are:
#' - bioSample: BioSample accession
#' - bioProject: BioProject accession
#' - source: Metadata annotation of the BioSample source species
#' NOT EXPORTED
#' @param virome A virome object
#' @param con A database connection object (to Serratus SQL).
#' @keywords internal
#' @return A data.frame with OTUs as rows and samples as columns, with values
#' being the normalized number of reads mapping to each OTU in each sample.
getAmpvisCounts <- function(virome, con) {

  counts <- getSpeciesCounts(virome, con)
  families <- virome %>%
    select(sotu, tax_family) %>%
    distinct()
  counts <- counts %>%
    left_join(families, by = c("sotu" = "sotu"))

  counts <- counts %>%
    pivot_wider(names_from = bio_sample, values_from = prop,
                values_fill = 0)

  # rename columns so ampvis recognizes them
  colnames(counts)[1] <- "otu"
  colnames(counts)[2] <- "family"

  # Make all other columns numeric
  counts[,3:ncol(counts)] <- apply(counts[,3:ncol(counts)], 2, as.numeric)

  counts$family[is.na(counts$family)] <- "null"
  # make missing numeric values 0 (some runs have 0 spots?)
  counts[is.na(counts)] <- 0

  # get metadata table
  metadata <- virome %>%
    select(bio_sample, bio_project, scientific_name) %>%
    distinct()

  return(list(counts, metadata))
}


# [END]






