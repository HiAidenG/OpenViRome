# viromeQueries.R
# This file contains functionality relating to querying the Serratus SQL and
# NCBI taxonomy for assembling a virome object.

#' @title getVirome
#' @description Returns a virome object given either:
#' 1) a taxonomic term and (optional) taxonomic level
#' 2) a character vector of SRA accessions
#' Note: the broader the taxonomic search term, the longer the query will
#' take.
#' @param tax A taxonomic term
#' @param sra A character vector of SRA accessions
#' @param con A connection to the Serratus database
#' @return A virome object
#' @export
#' @examples
#' con <- palmid::SerratusConnect()
#' getVirome(tax = "Meloidogyne", con = con)
#' getVirome(sra = c("SRR17756040", "SRR5942326"), con = con)
#' @importFrom dplyr tbl collect filter mutate group_by summarise n
#' @importFrom taxize classification downstream
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
      dplyr::mutate(node_coverage_norm = (node_coverage / spots) * 1e6) %>%
      dplyr::select(-spots)

    # Add taxonomic information
    virusSpecies <- virome %>% dplyr::pull(tax_species) %>% unique()
    ranks <- rep('phylum', length(virusSpecies))
    virusPhyla <- nameVecToRank(names=virusSpecies, taxRank=ranks)
    taxInfo <- tibble::tibble(tax_species = virusSpecies,
                              tax_phylum = virusPhyla)
    virome <- virome %>% dplyr::left_join(taxInfo, by = 'tax_species')
  }


  return(list(virome, runDF))
  # TODO: add support for sra accession input

}


#' @title taxLookup
#' @description Return a list of all runs processed by Serratus that match
#' tax. Used internally by getVirome.
#  NOT EXPORTED
#' @param tax A taxon defined in NCBI taxonomy. Must be type char.
#' @param con A connection to the Serratus database
#' @return A character vector of SRA accessions.
#' @import dplyr
taxLookup <- function(tax, con) {
  # Get ranking of taxonomic term
  class <- taxize::classification(tax, db = 'ncbi')[[1]]
  # Check if a species was provided
  rank <- class[class$name == tax, 'rank']
  if (is.null(rank)) {
    stop("Error: could not find taxonomic term in NCBI taxonomy database")
  }
  else if (rank == 'species') {
    searchTerms <- class[class$name == tax, 'id']
  }
  else {
    # Collect all child taxa
    taxid <- as.character(class[class$name == tax, 'id'])
    searchTerms <- taxize::downstream(taxid, db ='ncbi', downto = 'species')
    searchTerms <- searchTerms[[1]][,'childtaxa_id']
  }
  # Get all SRA accessions
  query <- tbl(con, "srarun") %>%
    dplyr::filter(tax_id %in% searchTerms) %>%
    dplyr::select(tax_id, scientific_name, run, spots) %>%
    dplyr::collect()
  return(query)
}

#' @title scientificNametoTaxID
#' @description Convert a scientific name to a taxonomic ID
#' @param name A character vector of scientific names
#' @return A character vector of taxonomic IDs
#' @import taxize
#' NOT EXPORTED. Used internally as a helper.
scientificNametoTaxID <- function(name) {
  class <- taxize::classification(name, db = 'ncbi')[[1]]
  taxid <- class[class$name == name, 'id']
  return(taxid)
}

#' @title nameVecToRank
#' @description Convert a vector of scientific names to the specified taxonomic
#' ranks for those names.
#' @param names A character vector of scientific names
#' @param taxRank A character vector of taxonomic ranks
#' @return A character vector of taxons
#' @import taxize
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
        taxize::classification(names[i], db = 'ncbi')[[1]]
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



# [END]
