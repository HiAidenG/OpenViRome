# viromeQueries.R
# This file contains functionality relating to querying the Serratus SQL and
# NCBI taxonomy for assembling a virome object.

#' @title getVirome
#'
#' @description Returns a virome object given either:
#' 1) a taxonomic term and (optional) taxonomic level
#' 2) a character vector of SRA accessions
#' Note: the broader the taxonomic search term, the longer the query will
#' take.
#'
#' @param tax A taxonomic term
#' @param sra A character vector of SRA accessions
#' @param con A connection to the Serratus database
#' @return A virome object
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' getVirome(tax = "Meloidogyne", con = con)
#' getVirome(sra = c("SRR17756040", "SRR5942326"), con = con)
#'
#' @importFrom taxizedb classification downstream
#' @importFrom dplyr filter select collect left_join pull
#'
#' @references
#' Wickham H, François R, Henry L, Müller K, Vaughan D (2023).
#' dplyr: A Grammar of Data Manipulation.
#' https://dplyr.tidyverse.org, https://github.com/tidyverse/dplyr.
#'
#' A. Babaian and R. C. Edgar (2022),
#'Ribovirus classification by a polymerase barcode sequence. PeerJ
#'
#' R. C. Edgar et al. (2021),
#' Petabase-scale sequence alignment catalyses viral discovery. Nature.
#' @export
getVirome <- function(tax = NULL, sra = NULL, con = NULL) {
  if (is.null(con)) {
    stop("Please provide a connection to the Serratus database
         (see palmid::SerratusConnect)")
  }
  if (is.null(tax) & is.null(sra)) {
    stop("Must provide either a taxonomic term or a character vector of
         SRA accessions")
  }
  if (!is.null(tax) & !is.null(sra)) {
    stop("Must provide either a taxonomic term or a character vector of
         SRA accessions")
  }
  if (!is.null(tax) & is.null(sra)) {
    runDF <- taxLookup(tax, con)

    # Check there wasn't an error
    if (is.null(runDF)) {
      stop("Error: could not find taxonomic term in NCBI taxonomy database")
    }

    runDF <- taxLookup(tax = tax, con = con)
    runDF <- data.frame(runDF)
    runs <- runDF %>% dplyr::pull(run)
    virome <- tbl(con, "palm_virome") %>%
      dplyr::filter(run %in% runs) %>%
      dplyr::collect()

    # Mutate runDF to include a column for whether that run was virus positive
    runDF <- runDF %>%
      dplyr::mutate(virus_positive = ifelse(run %in% virome$run, TRUE, FALSE))

    # Add normalized coverage
    virome <- virome %>%
      dplyr::left_join(runDF %>% dplyr::select(run, spots), by = 'run') %>%
      dplyr::mutate(node_coverage_norm = (node_coverage / spots) * 1e8) %>%
      dplyr::select(-spots)

    # Add taxonomic information
    virusSpecies <- virome %>% select(sotu) %>% unique() %>% pull(sotu)
    phyla <- getVirusTaxonomy(otu = virusSpecies, con = con)

    # Join the taxonomic information to the virome
    virome <- virome %>% dplyr::left_join(phyla, by = 'sotu',
                                          relationship = "many-to-many")

    # Rename NA to "Unclassified"
    virome <- virome %>% dplyr::mutate(tax_phylum = ifelse(is.na(tax_phylum),
                                                       "Unclassified",
                                                       tax_phylum))
  }


  return(list(virome, runDF))
  # TODO: add support for sra accession input

}

#' @title taxLookup
#'
#' @description Return a list of all runs processed by Serratus that match
#' tax. Used internally by getVirome.
#  NOT EXPORTED

#' @param tax A taxon defined in NCBI taxonomy. Must be type char.
#' @param con A connection to the Serratus database
#'
#' @return A character vector of SRA accessions.
#'
#' @import dplyr
taxLookup <- function(tax = NULL, con = NULL) {

  if (is.null(con)) {
    stop("Please provide a connection to the Serratus database
         (see palmid::SerratusConnect)")
  }

  if (is.null(tax)) {
    stop("Please provide a taxonomic term")
  }

  # Get ranking of taxonomic term
  class <- taxizedb::classification(tax, db = 'ncbi')[[1]]

  if (length(class) == 1) { # Should be length > 1 if valid term
    stop("Error: could not find taxonomic term in NCBI taxonomy database")
  }

  # Check if a species was provided
  rank <- class[class$name == tax, 'rank']
  if (rank == 'species') {
    searchTerms <- class[class$name == tax, 'id']
  }
  else {
    # Collect all child taxa
    taxid <- as.character(class[class$name == tax, 'id'])
    searchTerms <- taxizedb::downstream(taxid, db ='ncbi', downto = 'species')
    searchTerms <- searchTerms[[taxid]][["childtaxa_id"]]
  }
  # Get all SRA accessions
  query <- tbl(con, "srarun") %>%
    dplyr::filter(tax_id %in% searchTerms) %>%
    dplyr::select(tax_id, scientific_name, run, spots) %>%
    dplyr::collect()
  return(query)
}

#' @title getVirusTaxonomy
#'
#' @description Return the taxonomic information at a specified rank for
#' a given virus.
#'
#' @param virus A character vector of virus names
#' @param rank A character vector of taxonomic ranks
#' @param con A connection to the Serratus database
#'
#' @return A character vector of taxons
#'
#' @import dplyr
#'
#' @export
getVirusTaxonomy <- function(otu = NULL, con = NULL) {
  if (is.null(otu)) {
    stop("Must provide a character vector of virus names")
  }
  if (is.null(con)) {
    stop("Please provide a connection to the Serratus database
         (see palmid::SerratusConnect)")
  }

  virusClass <- tbl(con, "palm_tax") %>%
    dplyr::filter(sotu %in% otu) %>%
    dplyr::select(tax_phylum, sotu) %>%
    distinct() %>%
    dplyr::collect()


  return(virusClass)
}

#' @title nameVecToRank
#'
#' @description Convert a vector of scientific names to the specified taxonomic
#' ranks for those names.
#'
#' @param names A character vector of scientific names
#' @param taxRank A character vector of taxonomic ranks
#'
#' @return A character vector of taxons
#'
#' @import taxizedb
nameVecToRank <- function(names = NULL, taxRank = NULL) {
  if (is.null(names) | is.null(taxRank)) {
    stop("Must provide both a vector of names and a vector of taxonomic ranks")
  }
  if (length(names) != length(taxRank)) {
    stop("Length of names and taxRank must be equal")
  }

  taxons <- rep(NA, length(names))

  for (i in 1:length(names)) {
    if (!is.null(names[i])) {
      class <- tryCatch({
        taxizedb::classification(names[i], db = 'ncbi')[[1]]
      }, error = function(e) NULL)

      # Check if class is a dataframe or list and contains the rank column
      if (!is.null(class) && "data.frame" %in% class(class) && "rank" %in%
          names(class)) {
        # Filter the classification for the desired rank
        rankRow <- class[class$rank == taxRank[i], 'name', drop = FALSE]
        if (nrow(rankRow) > 0) {
          taxons[i] <- rankRow$name
        }
      }
    }
  }

  return(taxons)
}

#' @title getViromeSummary
#'
#' @description Return a number of summary statistics for a given virome.
#'
#' @param virome A virome object
#'
#' @return A list of summary statistics (# of unique sOTUs, mean normalized
#' coverage, median normalized coverage, max normalized coverage, virus
#' positive runs, total runs processed).
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax="Canis", con=con)
#' summary <- getViromeSummary(virome = virome)
#'
#' @import dplyr
#'
#' @export
getViromeSummary <- function(virome = NULL) {

  if (is.null(virome)) {
    stop("Please provide a virome object (see getVirome)")
  }

  runData <- virome[[2]]
  virome <- virome[[1]]

  # Get unique sOTU count
  sotuCount <- virome %>% dplyr::select(sotu) %>% unique() %>% nrow()

  # Get median normalized coverage
  medianNormCov <- virome %>% dplyr::select(node_coverage_norm) %>%
    dplyr::summarise(median = median(node_coverage_norm)) %>%
    dplyr::pull(median)

  # Get mean normalized coverage
  meanNormCov <- virome %>% dplyr::select(node_coverage_norm) %>%
    dplyr::summarise(mean = mean(node_coverage_norm)) %>%
    dplyr::pull(mean)

  # Get max normalized coverage
  maxNormCov <- virome %>% dplyr::select(node_coverage_norm) %>%
    dplyr::summarise(max = max(node_coverage_norm)) %>%
    dplyr::pull(max)

  # Get virus positive runs
  virusPos <- virome %>% dplyr::select(run) %>% unique() %>% nrow()

  # Get total runs processed
  totalRuns <- runData %>% dplyr::select(run) %>% unique() %>% nrow()

  returnVec <- c(sotuCount, meanNormCov, medianNormCov, maxNormCov, virusPos,
                 totalRuns)
  names(returnVec) <- c("Unique sOTUs", "Mean Coverage", "Median Coverage",
                        "Max Coverage", "Virus Positive Runs", "Runs Processed")

  return(returnVec)
}

#' @title getAvailablePhyla
#'
#' @description Returns a list of phyla that are present in the virome.
#' Mostly for use by the shiny frontend.
#'
#' @param virome A virome object
#'
#' @return A character vector of phyla
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax="Canis", con=con)
#' phyla <- getAvailablePhyla(virome = virome)
#'
#' @import dplyr
#'
#' @export
getAvailablePhyla <- function(virome = NULL) {
  if (is.null(virome)) {
    stop("Please provide a virome object (see getVirome)")
  }

  virome <- virome[[1]]

  phyla <- virome %>% dplyr::select(tax_phylum) %>% unique() %>%
    dplyr::pull(tax_phylum)

  return(phyla)
}

#' @title viromeFormatCheck
#'
#' @description Checks if the virome object is in the correct format. Mostly for
#' use with the shiny frontend.
#'
#' @param virome A virome object
#'
#' @return A boolean value
#'
#' @examples
#' con <- palmid::SerratusConnect()
#' virome <- getVirome(tax="Canis", con=con)
#' if (viromeFormatCheck(virome = virome)) {
#'  print("Virome is in the correct format.")
#' }
#'
#' @import dplyr
viromeFormatCheck <- function(virome = NULL) {
  if (is.null(virome)) {
    stop("Please provide a virome object (see getVirome)")
  }

  if (!is.data.frame(virome)) {
    return(FALSE)
  }

  if (!all(c("sotu", "tax_phylum", "tax_species", "node_coverage_norm",
             "run", "bio_sample", "bio_project", "scientific_name",
             "gb_pid") %in%
             colnames(virome))) {
    return(FALSE)
  }

  return(TRUE)
}

# [END]
